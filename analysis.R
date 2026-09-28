# ============================================================
# Analysis: ww_national_mean
# Author: Moritz Wagner (moritz.wagner@bag.admin.ch)
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(gratia)
library(INLA)
library(MASS)
library(mgcv)

# -------------------------------------------------------------------------
# Config
# -------------------------------------------------------------------------

PATHOGEN <- "sars-cov-2"
# "influenza_a", "influenza_b", "respiratory_syncytial_virus", "sars-cov-2"

N_SAMPLES <- 2000
SEED <- 123
INLA_TIME_STEP <- 3L
USE_POPULATION_WEIGHTS <- TRUE
TREND_DAYS <- 14L
TREND_STABLE_THRESHOLD <- 0.10
RAW_SMOOTH_WINDOW <- 7L

# Shift applied to viral_load before the Gamma fit, so exact zeros
# (non-detects) are valid under a strictly-positive family. NOTE: 1 is a
# nominal placeholder given viral_load ~ 1e10-1e15 -- if this is meant to
# encode "indistinguishable from below-LOQ", it should probably sit near
# the assay's actual LOQ on the raw scale, not a flat +1.
SHIFT <- 1

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

get_site_info <- function(dat) {
  last_day <- max(dat$day)
  last_date <- max(dat$date)  # forced global, not site-specific, so every
  # site's grid extends to the same end point
  dat |>
    group_by(wwtp) |>
    summarise(
      population = median(population, na.rm = TRUE),
      first_day = min(day), last_day = last_day,
      first_date = min(date), last_date = last_date,
      .groups = "drop"
    )
}

get_prediction_days <- function(dat, inla_time_step = INLA_TIME_STEP) {
  tibble(time = sort(unique(dat$time))) |>
    mutate(
      day = (time - 1L) * inla_time_step + floor(inla_time_step / 2) + 1L,
      day = pmin(day, max(dat$day)),
      date = min(dat$date) + day - 1L
    )
}

build_site_time_grid <- function(dat, prediction_days, site_info) {
  crossing(prediction_days, wwtp = levels(dat$wwtp)) |>
    mutate(
      wwtp = factor(wwtp, levels = levels(dat$wwtp)),
      wwtp_id = as.integer(wwtp),
      time_ch = time, time_wwtp = time
    ) |>
    left_join(site_info, by = "wwtp") |>
    # don't extrapolate a site's trajectory beyond when it was surveilled
    filter(day >= first_day, day <= last_day) |>
    group_by(time) |>
    mutate(weight = population / sum(population)) |>
    ungroup()
}

# Aggregates a set of already-predicted site-level draws into a
# population-weighted national series. Unchanged from before -- now called
# with the `draws`/`grid` produced by predict_*_site() below.
aggregate_national_draws <- function(response_draws, grid, model_name, weighted = TRUE) {
  days <- sort(unique(grid$day))
  index <- grid |> distinct(day, date) |> filter(day %in% days) |> arrange(day)

  national_draws <- do.call(cbind, lapply(days, function(d) {
    ii <- which(grid$day == d)
    w <- if (weighted) grid$weight[ii] else rep(1 / length(ii), length(ii))
    drop(response_draws[, ii, drop = FALSE] %*% w)
  }))

  summary <- index |>
    mutate(
      estimate = apply(national_draws, 2, median),
      lower = apply(national_draws, 2, quantile, probs = 0.025),
      upper = apply(national_draws, 2, quantile, probs = 0.975),
      model = model_name
    )

  list(summary = summary, draws = national_draws, index = index, model = model_name)
}

# Site-level (WWTP x time) summary directly from the draws -- for comparing
# model fit against raw per-site observations.
summarise_site_draws <- function(response_draws, grid, model_name) {
  grid |>
    ungroup() |>
    mutate(
      estimate = apply(response_draws, 2, median),
      lower = apply(response_draws, 2, quantile, probs = 0.025),
      upper = apply(response_draws, 2, quantile, probs = 0.975),
      model = model_name
    ) |>
    dplyr::select(wwtp, day, date, estimate, lower, upper, model)
}

# Change over exactly `trend_days` by interpolating each posterior draw on
# the (3-day) model time grid, worked out on the log1p scale so differences
# become proportional changes.
summarise_trend_draws <- function(draws, index, model_name, trend_days = TREND_DAYS,
                                  stable_threshold = TREND_STABLE_THRESHOLD) {

  x <- index$day
  target_x <- x - trend_days

  log_draws <- log1p(pmax(draws, 0))
  lagged_log_draws <- t(apply(log_draws, 1, function(z) {
    approx(x = x, y = z, xout = target_x, rule = 1)$y
  }))

  log_change <- log_draws - lagged_log_draws
  pct_change <- 100 * expm1(log_change)
  increase_cutoff <- log1p(stable_threshold)
  decrease_cutoff <- log1p(-stable_threshold)

  index |>
    mutate(
      estimate = apply(pct_change, 2, median, na.rm = TRUE),
      lower = apply(pct_change, 2, quantile, probs = 0.025, na.rm = TRUE),
      upper = apply(pct_change, 2, quantile, probs = 0.975, na.rm = TRUE),
      p_increasing = colMeans(log_change > increase_cutoff, na.rm = TRUE),
      p_stable = colMeans(log_change >= decrease_cutoff & log_change <= increase_cutoff, na.rm = TRUE),
      p_decreasing = colMeans(log_change < decrease_cutoff, na.rm = TRUE),
      direction = case_when(
        p_increasing >= 0.8 ~ "increasing",
        p_decreasing >= 0.8 ~ "decreasing",
        p_stable >= 0.8 ~ "stable",
        TRUE ~ "uncertain"
      ),
      model = model_name, trend_days = trend_days, stable_threshold = stable_threshold
    )
}

summarise_trend_draws_ols <- function(draws, index, model_name,
                                      trend_days = TREND_DAYS,
                                      stable_threshold = TREND_STABLE_THRESHOLD) {
  x_all <- index$day
  log_draws <- log1p(pmax(draws, 0))  # n_samples x n_days

  # per-day equivalent of the existing +/-stable_threshold-over-trend_days band
  increase_cutoff <- log1p(stable_threshold) / trend_days
  decrease_cutoff <- log1p(-stable_threshold) / trend_days

  map_dfr(seq_along(x_all), function(i) {
    target_day <- x_all[i]
    win <- which(x_all >= target_day - trend_days & x_all <= target_day)

    if (length(win) < 3) {  # too few grid points for a meaningful slope
      return(tibble(day = target_day, date = index$date[i],
                    estimate = NA_real_, lower = NA_real_, upper = NA_real_,
                    p_increasing = NA_real_, p_stable = NA_real_, p_decreasing = NA_real_,
                    direction = NA_character_))
    }

    xw <- x_all[win]
    yw <- log_draws[, win, drop = FALSE]        # n_samples x length(win)
    xc <- xw - mean(xw)                          # centre x -> intercept drops out
    slopes <- as.numeric(yw %*% xc) / sum(xc^2)   # one slope per draw, fully vectorized

    pct <- 100 * expm1(slopes * trend_days)       # express as "% over trend_days", same units as before

    tibble(
      day = target_day, date = index$date[i],
      estimate = median(pct), lower = quantile(pct, 0.025), upper = quantile(pct, 0.975),
      p_increasing = mean(slopes > increase_cutoff),
      p_stable = mean(slopes >= decrease_cutoff & slopes <= increase_cutoff),
      p_decreasing = mean(slopes < decrease_cutoff),
      direction = case_when(
        mean(slopes > increase_cutoff) >= 0.8 ~ "increasing",
        mean(slopes < decrease_cutoff) >= 0.8 ~ "decreasing",
        mean(slopes >= decrease_cutoff & slopes <= increase_cutoff) >= 0.8 ~ "stable",
        TRUE ~ "uncertain"
      )
    )
  }) |>
    mutate(model = model_name, trend_days = trend_days, stable_threshold = stable_threshold)
}

summarise_raw_data <- function(raw_data, weighted = TRUE, window = RAW_SMOOTH_WINDOW) {
  first_date <- min(raw_data$date)
  last_date <- max(raw_data$date)

  raw_data |>
    group_by(wwtp) |>
    # daily grid so every active WWTP contributes even on days it wasn't sampled
    complete(date = seq.Date(min(date), last_date, by = "day")) |>
    arrange(date) |>
    mutate(
      population = median(population, na.rm = TRUE),
      viral_load_smoothed = slider::slide_index_dbl(
        viral_load, date, ~ median(.x, na.rm = TRUE),
        .before = window - 1, .after = 0, .complete = TRUE
      )
    ) |>
    ungroup() |>
    mutate(day = as.integer(date - first_date) + 1L) |>
    group_by(day, date) |>
    summarise(
      estimate = if (weighted) {
        bag.epi.stats::calculate_summary_statistic(viral_load_smoothed, weights = population, method = "mean")
      } else {
        bag.epi.stats::calculate_summary_statistic(viral_load_smoothed, method = "mean")
      },
      model = "WWTP 14-day rolling median + national mean",
      .groups = "drop"
    )
}

fit_gam <- function(dat) {
  gam(
    y ~
      s(wwtp, bs = "re") +               # persistent WWTP-specific baseline
      s(day, bs = "tp", k = 40) +        # shared temporal component
      s(day, wwtp, bs = "sz", k = 10),   # WWTP-specific deviations from it
    data = dat, method = "REML", select = TRUE
  )
}

# -------------------------------------------------------------------------
# Site-level draws (WWTP x time) -- the shared building block for both
# the national aggregate and the per-site fit-vs-raw comparison.
# -------------------------------------------------------------------------

predict_gam_site <- function(m_gam, gam_grid, n_samples = N_SAMPLES, seed = SEED) {
  Xp <- predict(m_gam, newdata = gam_grid, type = "lpmatrix")

  set.seed(seed)
  beta_draws <- MASS::mvrnorm(n = n_samples, mu = coef(m_gam), Sigma = vcov(m_gam, unconditional = TRUE))
  response_draws <- expm1(beta_draws %*% t(Xp))

  list(draws = response_draws, grid = gam_grid)
}

predict_gam_national <- function(m_gam, gam_grid, n_samples = N_SAMPLES,
                                 seed = SEED, weighted = TRUE) {
  site <- predict_gam_site(m_gam, gam_grid, n_samples, seed)
  aggregate_national_draws(site$draws, site$grid, model_name = "GAM", weighted = weighted)
}

build_inla_data <- function(dat, inla_grid) {
  bind_rows(
    dat |> mutate(is_prediction = FALSE),
    inla_grid |>
      transmute(
        date, wwtp, viral_load = NA_real_, population, y = NA_real_,
        day, time, time_ch, time_wwtp, wwtp_id, is_prediction = TRUE
      )
  )
}

fit_inla <- function(inla_data) {
  inla_data <-
    inla_data |>
    dplyr::mutate(
      wwtp_id2 = wwtp_id,   # only used by the commented RW1+group spec below
      viral_load = viral_load + SHIFT
      # Not needed with transformation in raw data
      # + SHIFT
    )

  # formula_inla <- y ~ 1 +
  formula_inla <- viral_load ~ 1 +
    f(wwtp_id, model = "iid") +
    f(time_ch, model = "rw2", scale.model = TRUE, constr = TRUE)+
    f(time_wwtp, model = "ar1", replicate = wwtp_id)
  # f(time_ch, model="rw2", scale.model=TRUE, constr=TRUE, hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01))))+
  # f(time_wwtp, model="rw1", scale.model=TRUE, constr=TRUE, group=wwtp_id2, control.group=list(model="iid"),
  #   hyper=list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

  inla(
    # formula_inla, family = "gaussian", data = inla_data,
    formula_inla, family = "gamma", data = inla_data,
    control.predictor = list(compute = TRUE),
    control.compute = list(config = TRUE, waic = TRUE)
  )
}

predict_inla_site <- function(m_inla, inla_data, inla_grid, n_samples = N_SAMPLES,
                              seed = SEED, shift = SHIFT) {
  set.seed(seed)
  samples <- inla.posterior.sample(n = n_samples, result = m_inla)

  latent_names <- rownames(samples[[1]]$latent)
  idx_predictor <- grep("^Predictor:", latent_names)
  predictor_rows <- as.integer(sub("^Predictor:", "", latent_names[idx_predictor]))

  pred_rows <- which(inla_data$is_prediction)
  use_idx <- match(pred_rows, predictor_rows)
  if (anyNA(use_idx)) stop("Could not match all prediction rows to INLA predictor latent indices.")

  eta_draws <- t(vapply(samples, function(s) s$latent[idx_predictor[use_idx], 1], numeric(length(pred_rows))))
  # log link -> exp(eta) recovers the mean of what was fit, i.e. viral_load + shift.
  # Undo the shift here so response_draws is back on the viral_load scale.
  response_draws <- exp(eta_draws) - shift

  list(draws = response_draws, grid = inla_grid)
}

predict_inla_national <- function(m_inla, inla_data, inla_grid, n_samples = N_SAMPLES,
                                  seed = SEED, weighted = TRUE, shift = SHIFT) {
  site <- predict_inla_site(m_inla, inla_data, inla_grid, n_samples, seed, shift)
  aggregate_national_draws(site$draws, site$grid, model_name = "INLA", weighted = weighted)
}

# -------------------------------------------------------------------------
# Main pipeline: load, prepare, fit
# -------------------------------------------------------------------------

wastewater <- data.table::fread(
  "https://idd.bag.admin.ch/api/v1/export/current/RESPVIRUSES_wastewater/csv"
) |> tibble::as_tibble()

dat <- wastewater |>
  filter(!is.na(value), georegion2 != "liechtenstein", valueCategory == PATHOGEN) |>
  # Sketchy handling of zeros
  # dplyr::group_by(georegion2) |>
  # dplyr::mutate(conc = dplyr::case_when(conc>0 ~ conc,
  #                                       # TODO Very arbitrary!! Think of something better here
  #                                       TRUE ~ min(conc[conc != 0])/100),
  #               # ww_concentration * 1000 * ww_flow / (ww_population / inz_ref)
  #               value = 1000*conc*flow/(pop/1e5)) |>
  # dplyr::ungroup() |>
  transmute(date = as.Date(temporal), wwtp = factor(georegion2), viral_load = value, population = pop) |>
  mutate(
    y = log1p(viral_load),
    day = as.integer(date - min(date)) + 1L,
    time = floor((day - 1L) / INLA_TIME_STEP) + 1L,
    time_ch = time, time_wwtp = time,
    wwtp_id = as.integer(wwtp),
  )

site_info <- get_site_info(dat)
site_time_grid <- build_site_time_grid(dat, get_prediction_days(dat), site_info)

m_gam <- fit_gam(dat)
gam_site <- predict_gam_site(m_gam, site_time_grid)
gam_result <- aggregate_national_draws(gam_site$draws, gam_site$grid, model_name = "GAM", weighted = USE_POPULATION_WEIGHTS)
gam_wwtp <- summarise_site_draws(gam_site$draws, gam_site$grid, model_name = "GAM")

inla_data <- build_inla_data(dat, site_time_grid)
m_inla <- fit_inla(inla_data)
inla_site <- predict_inla_site(m_inla, inla_data, site_time_grid)
inla_result <- aggregate_national_draws(inla_site$draws, inla_site$grid, model_name = "INLA", weighted = USE_POPULATION_WEIGHTS)
inla_wwtp <- summarise_site_draws(inla_site$draws, inla_site$grid, model_name = "INLA")

raw_ch <- summarise_raw_data(dat, weighted = USE_POPULATION_WEIGHTS)

models_ch <- bind_rows(gam_result$summary, inla_result$summary, raw_ch)
model_trends <- bind_rows(
  summarise_trend_draws(gam_result$draws, gam_result$index, first(gam_result$summary$model)),
  summarise_trend_draws(inla_result$draws, inla_result$index, first(inla_result$summary$model))
)
model_trends_ols <- bind_rows(
  summarise_trend_draws_ols(gam_result$draws, gam_result$index, first(gam_result$summary$model)),
  summarise_trend_draws_ols(inla_result$draws, inla_result$index, first(inla_result$summary$model))
)

comparison_subtitle <- if (USE_POPULATION_WEIGHTS) {
  "Population-weighted national signal: GAM vs INLA"
} else {
  "Unweighted national signal: GAM vs INLA"
}

latest_trend_summary <- model_trends |>
  filter(date == max(date, na.rm = TRUE)) |>
  dplyr::select(model, date, estimate, lower, upper, p_increasing, p_stable, p_decreasing, direction)

# -------------------------------------------------------------------------
# Diagnostics
# -------------------------------------------------------------------------

summary(m_gam)
gratia::draw(m_gam, residuals = TRUE)
gratia::appraise(m_gam)
gam.check(m_gam)

summary(m_inla)
m_inla$summary.hyperpar

# GAM terms, split by smooth -- se.fit lets us build a CI for the
# (until now missing) WWTP baseline term.
gam_terms <- predict(m_gam, newdata = site_time_grid, type = "terms", se.fit = TRUE)
colnames(gam_terms$fit)  # confirm exact term labels before indexing below --
# should include "s(wwtp)", "s(day)", "s(day,wwtp)"
# alternative using gratia
gratia::smooth_estimates(m_gam)

grid_terms <- site_time_grid |>
  mutate(
    site_baseline = gam_terms$fit[, "s(wwtp)"],
    site_baseline_se = gam_terms$se.fit[, "s(wwtp)"],
    national_trend = gam_terms$fit[, "s(day)"],
    site_deviation = gam_terms$fit[, "s(day,wwtp)"]
  )

p1 <-
  dat |>
  ggplot(aes(time, viral_load, colour = wwtp))+
  geom_line(alpha = .5)+
  scale_y_log10()+
  theme_classic() +
  labs(title = "Raw viral loads by WWTPs")

p2 <-
  ggplot() +
  geom_ribbon(data = m_inla$summary.random$time_ch,
              aes(ID, ymin = `0.025quant`, ymax = `0.975quant`), alpha  = .2)+
  geom_line(data = m_inla$summary.random$time_ch,
            aes(ID, mean)) +
  # geom_ribbon(data = m_inla$summary.random$time_wwtp |>
  #             group_by(ID) |>
  #             mutate(wwtp = levels(dat$wwtp)),
  #           aes(ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = wwtp),
  #           alpha = 0.2)+
  geom_line(data = m_inla$summary.random$time_wwtp |>
              group_by(ID) |>
              mutate(wwtp = levels(dat$wwtp)),
            aes(ID, mean, colour = wwtp),
            alpha = 0.5)+
  theme_classic() +
  labs(title = "INLA: National RW2 vs WWTP-specific AR1 components")

p3 <-
  ggplot() +
  geom_line(data =  grid_terms |>
              distinct(day, date, national_trend),
            aes(date, national_trend)
  ) +
  geom_line(data = grid_terms,
            aes(date, site_deviation, colour = wwtp, group = wwtp), alpha = 0.5) +
  theme_classic() +
  labs(title = "GAM: National spline vs WWTP-specific splines")

gam_baseline <- grid_terms |>
  distinct(wwtp, site_baseline, site_baseline_se) |>
  transmute(
    wwtp, estimate = site_baseline,
    lower = site_baseline - 1.96 * site_baseline_se,
    upper = site_baseline + 1.96 * site_baseline_se,
    model = "GAM"
  )

inla_baseline <- m_inla$summary.random$wwtp_id |>
  mutate(wwtp = levels(dat$wwtp)) |>
  transmute(wwtp, estimate = mean, lower = `0.025quant`, upper = `0.975quant`, model = "INLA")

p4 <-
  list(inla_baseline, gam_baseline) |>
  purrr::map(~.x |>
  ggplot(aes(wwtp, estimate, colour = wwtp)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
  # coord_flip() +
  # facet_wrap(~model, scales = "free_y", nrow = 2)+
  theme_classic() +
  labs(title = paste0(first(.x$model), ": Persistent WWTP baseline"),
       x = NULL, y = "Offset from overall level")
  )

((p1/p2/p3)|(plot_spacer()/p4))+
  plot_layout(widths = c(2, 1), guides = "collect")

# p1 <-
#   ggplot() +
#   geom_line(data = m_inla$summary.random$time_ch,
#             aes(ID, mean)) +
#   theme_classic() +
#   labs(title = "INLA: National RW2")
#
# p2 <-
#   ggplot() +
#   geom_hline(yintercept = 0, linetype = 2)+
#   geom_line(data = m_inla$summary.random$time_wwtp |>
#               group_by(ID) |>
#               mutate(wwtp = levels(dat$wwtp)),
#             aes(ID, mean, colour = wwtp))+
#   theme_classic() +
#   labs(title = "WWTP-specific RW1 components")
# p1/p2

# ---- NEW: persistent WWTP baseline, GAM (re) vs INLA (iid) -----------
# NOTE: not on the same numeric scale -- GAM's is on the log1p(viral_load)
# / identity-link scale, INLA's is on the log(mean viral_load) / gamma
# log-link scale. Comparable in spirit (both ~ log-fold offset from the
# overall level) but not directly overlayable as identical units.


# ---- NEW: model fit vs raw data, per WWTP -----------------------------
# Uses the site-level draws directly (gam_wwtp / inla_wwtp), on the
# viral_load scale, against the raw points for the same site.

wwtp_selection <- c("aire", "werdhoelzli", "chur", "laupen")
bind_rows(gam_wwtp, inla_wwtp) |>
  dplyr::filter(wwtp %in% wwtp_selection) |>
  ggplot() +
  geom_point(data = dat |>
               dplyr::filter(wwtp %in% wwtp_selection) ,
             aes(date, viral_load), size = 0.4, alpha = .8, colour = "grey40") +
  geom_ribbon(aes(date, ymin = lower, ymax = upper, fill = model), alpha = .2, colour = NA) +
  geom_line(aes(date, estimate, colour = model), linewidth = 0.6) +
  scale_y_log10() +
  facet_wrap(~ wwtp, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(title = "Model fit vs raw data, by WWTP", x = NULL, y = "Viral load")

# -------------------------------------------------------------------------
# Plots
# -------------------------------------------------------------------------

p1 <- ggplot(models_ch
             # |>
               # dplyr::filter(date >="2023-01-01", model != "GAM")
             ) +
  geom_ribbon(aes(date, ymin = lower, ymax = upper, fill = model), alpha = 0.12, colour = NA) +
  geom_line(aes(date, estimate, colour = model), linewidth = 1, alpha = 0.6) +
  # scale_y_log10() +
  theme_classic() +
  # scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_colour_manual(values = c("#D32F2F", "#388E3C", "#1976D2")) +
  scale_fill_manual(values = c("#D32F2F", "#388E3C", "#1976D2")) +
  theme(legend.position = "right") +
  labs(x = NULL, y = "Viral load", colour = "Model", fill = "Model",
       title = "Swiss wastewater signal", subtitle = comparison_subtitle)

p2 <- ggplot(model_trends
             # |>
               # dplyr::filter(date >="2023-01-01", model != "GAM")
               ,
             aes(date, estimate, colour = model, fill = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 100 * c(-TREND_STABLE_THRESHOLD, TREND_STABLE_THRESHOLD), linetype = "dotted") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("#D32F2F", "#388E3C")) +
  scale_fill_manual(values = c("#D32F2F", "#388E3C")) +
  theme_classic() +
  labs(x = NULL, y = glue::glue("Change over {TREND_DAYS} days (%)"),
       colour = "Model", fill = "Model",
       title = "Change in Swiss wastewater signal",
       subtitle = glue::glue("Dotted lines show the \u00b1{100 * TREND_STABLE_THRESHOLD}% stable range"))

p3 <- ggplot(model_trends_ols
             # |>
             # dplyr::filter(date >="2023-01-01", model != "GAM")
             ,
             aes(date, estimate, colour = model, fill = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 100 * c(-TREND_STABLE_THRESHOLD, TREND_STABLE_THRESHOLD), linetype = "dotted") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("#D32F2F", "#388E3C")) +
  scale_fill_manual(values = c("#D32F2F", "#388E3C")) +
  theme_classic() +
  labs(x = NULL, y = glue::glue("Change over {TREND_DAYS} days (%)"),
       colour = "Model", fill = "Model",
       title = "Change in Swiss wastewater signal",
       subtitle = glue::glue("Dotted lines show the \u00b1{100 * TREND_STABLE_THRESHOLD}% stable range"))

p1 / p2

# ============================================================
# Nowcast backtest: GAM vs INLA
# ------------------------------------------------------------
# For a handful of scenario cutoff dates (wave start, peak,
# declining, low season), truncate the data to what would have
# been available at that point, refit both models using the
# functions above, and compare their nowcast for the preceding
# weeks against the "final" estimate from the full-data fit.
# ============================================================

# Eyeballed off the plot -- swap in the exact dates once you can see the
# underlying daily values. Pick each cutoff a few days *into* the phase so
# there's real signal for the model to react to, not right at a turning point.
# BACKTEST_CUTOFFS <- tribble(
#   ~phase,        ~cutoff_date,
#   "wave_start",  as.Date("2025-10-15"),
#   "peak",        as.Date("2025-12-01"),
#   "declining",   as.Date("2026-02-15"),
#   "low_season",  as.Date("2026-06-15")
# )
#
# NOWCAST_WINDOW_DAYS <- 21L   # how many days back from each cutoff to score
# PLOT_WINDOW_DAYS <- 100L      # how many days back from each cutoff to plot
#
# refit_at_cutoff <- function(cutoff_date, full_dat, inla_time_step = INLA_TIME_STEP,
#                             weighted = USE_POPULATION_WEIGHTS) {
#   dat_trunc <- full_dat |>
#     filter(date <= cutoff_date) |>
#     mutate(
#       day = as.integer(date - min(date)) + 1L,
#       time = floor((day - 1L) / inla_time_step) + 1L,
#       time_ch = time, time_wwtp = time,
#       wwtp_id = as.integer(wwtp)
#     )
#
#   grid_trunc <- build_site_time_grid(dat_trunc, get_prediction_days(dat_trunc, inla_time_step), get_site_info(dat_trunc))
#
#   gam_trunc <- predict_gam_national(fit_gam(dat_trunc), grid_trunc, weighted = weighted)$summary |>
#     mutate(model = "GAM")
#
#   inla_data_trunc <- build_inla_data(dat_trunc, grid_trunc)
#   inla_trunc <- predict_inla_national(fit_inla(inla_data_trunc), inla_data_trunc, grid_trunc, weighted = weighted)$summary |>
#     mutate(model = "INLA")
#
#   bind_rows(gam_trunc, inla_trunc) |> mutate(cutoff_date = cutoff_date)
# }
#
# backtest_results <- map(BACKTEST_CUTOFFS$cutoff_date, refit_at_cutoff, full_dat = dat) |>
#   set_names(BACKTEST_CUTOFFS$phase) |>
#   bind_rows(.id = "phase")
#
# final_estimates <- models_ch |>
#   filter(model %in% c("GAM", "INLA")) |>
#   dplyr::select(date, model, final_estimate = estimate, final_lower = lower, final_upper = upper)
#
# backtest_comparison <- backtest_results |>
#   inner_join(final_estimates, by = c("date", "model")) |>
#   group_by(phase, cutoff_date) |>
#   filter(date >= cutoff_date - NOWCAST_WINDOW_DAYS, date <= cutoff_date) |>
#   ungroup() |>
#   mutate(
#     days_before_cutoff = as.integer(cutoff_date - date),
#     pct_error = 100 * (estimate - final_estimate) / final_estimate,
#     covered = final_estimate >= lower & final_estimate <= upper
#   )
#
# backtest_summary <- backtest_comparison |>
#   group_by(phase, model) |>
#   summarise(
#     mean_pct_error = mean(pct_error, na.rm = TRUE),  # bias: + = nowcast overshoots what was later confirmed
#     mae_pct = mean(abs(pct_error), na.rm = TRUE),     # typical size of the miss
#     coverage = mean(covered, na.rm = TRUE),           # should be ~0.95 if uncertainty is well-calibrated
#     .groups = "drop"
#   ) |>
#   arrange(phase, model)
#
# print(backtest_summary)
#
# plot_backtest <- function(backtest_results, final_estimates, cutoffs, plot_window_days = PLOT_WINDOW_DAYS) {
#   final_windowed <- pmap_dfr(cutoffs, function(phase, cutoff_date) {
#     final_estimates |>
#       filter(date >= cutoff_date - plot_window_days, date <= cutoff_date + 14) |>
#       mutate(phase = phase, cutoff_date = cutoff_date, source = "final (all data)") |>
#       rename(estimate = final_estimate, lower = final_lower, upper = final_upper)
#   })
#
#   combined <- bind_rows(
#     final_windowed |> dplyr::select(phase, cutoff_date, date, model, estimate, lower, upper, source),
#     backtest_results |> mutate(source = "nowcast (at cutoff)") |>
#       dplyr::filter(date >= cutoff_date - 365) |>
#       dplyr::select(phase, cutoff_date, date, model, estimate, lower, upper, source)
#   ) |>
#     mutate(phase = factor(phase, levels = cutoffs$phase))
#
#   ggplot(combined, aes(date, estimate, colour = source, fill = source)) +
#     geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
#     geom_line(linewidth = 1) +
#     geom_vline(aes(xintercept = cutoff_date), linetype = "dashed", colour = "grey40") +
#     facet_grid(model ~ phase, scales = "free_x") +
#     scale_y_log10() +
#     theme_classic() +
#     theme(legend.position = "bottom") +
#     labs(x = NULL, y = "Viral load", colour = NULL, fill = NULL,
#          title = "Nowcast at cutoff vs. final (fully revised) estimate",
#          subtitle = "Dashed vertical line = cutoff date")
# }
#
# p_backtest <- plot_backtest(backtest_results, final_estimates, BACKTEST_CUTOFFS)
# p_backtest
