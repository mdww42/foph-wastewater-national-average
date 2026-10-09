# ============================================================
# Analysis: Swiss-wide wastewater average -- GAM vs INLA models
# Ticket:   MTDATA-5048
# Author:   Moritz Wagner (moritz.wagner@bag.admin.ch)
#
# Run from this folder. Functions live in R/.
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(mgcv)
library(INLA)

current_dir <- fs::path_dir(rstudioapi::documentPath())
walk(list.files(fs::path(current_dir, "R"), pattern = "\\.R$", full.names = TRUE), source)

# -------------------------------------------------------------------------
# Config
# -------------------------------------------------------------------------

PATHOGEN <- "sars-cov-2"
# "influenza_a", "influenza_b", "respiratory_syncytial_virus", "sars-cov-2"

N_SAMPLES <- 2000                 # posterior draws per model
SEED <- 123
TIME_STEP <- 7L                   # days per INLA time bin (also the prediction grid); 3-day
                                  # bins fit better (WAIC ~1800 lower) but take ~5x longer
USE_POPULATION_WEIGHTS <- TRUE
TREND_DAYS <- 21L
TREND_STABLE_THRESHOLD <- 0.10    # +/- 10% over TREND_DAYS counts as "stable"
RAW_SMOOTH_WINDOW <- 7L           # days, for the model-free reference

# INLA models to fit and compare, label = time structure of the WWTP
# deviations (see fit_inla()):
#   "ar1": pulled back to 0 (fitted weekly rho ~0.8 sars-cov-2, ~0.7
#          influenza A); best WAIC for both
#   "rw1": random walk as in Riou et al. (multilevel_wbe); deviations persist
#          after a plant's last sample. Same deviations in the backtest;
#          add `INLA_RW1 = "rw1"` as a sensitivity check.
INLA_MODELS <- c(INLA_AR1 = "ar1")

# National series per INLA model, all from the same posterior draws:
#   "<model>":               population-weighted arithmetic mean of plant curves
#   "<model> geometric":     population-weighted geometric mean (log scale)
#   "<model> national term": the model's national RW2, equal-weight geometric

# Positive samples at or below the LOQ (sars-cov-2 ~16%, influenza A ~48%).
# The lab considers them usable concentrations.
#   "censored":       INLA treats them as left-censored at the LOQ in load
#                     units (see add_loq()). For influenza A that left ~88% of
#                     samples censored and the low season without a lower
#                     bound (intervals spanning 1e5-1e13).
#   "observed":       used as measured
#   "observed_extra": used as measured, plus an estimated extra spread
BLQ_HANDLING <- "observed_extra"

# Zero loads are non-detects (sars-cov-2 45, influenza A ~41%). On the log
# scale a zero sits ~30 units below typical values: for sars-cov-2, 45 zeros
# caused 82% of the GAM residual variance (sig2 3.25 vs robust 0.34),
# inflating its mean correction ~5x and pulling the fitted curves down.
#   "censor": INLA treats them as left-censored at the plant's detection
#             floor (see add_loq()): "below what the lab detects"
#   "drop":   INLA ignores them; biases low periods upwards
# The GAM always drops zeros; the raw reference always uses all samples.
ZERO_HANDLING <- "censor"

# Samples more than OUTLIER_FACTOR x away from the median of the plant's
# samples within +/- OUTLIER_WINDOW_DAYS are excluded from the models and the
# raw reference, but shown (flagged) in the raw-data plots.
OUTLIER_WINDOW_DAYS <- 7L
OUTLIER_FACTOR <- 10

WWTP_SELECTION <- c("aire", "werdhoelzli", "chur", "laupen")  # for the fit-vs-data plot

RUN_BACKTEST <- FALSE             # refits all models per cutoff -- slow
EXPORT_FIGURES <- FALSE           # PNGs for presentation/ (see end of script)

# -------------------------------------------------------------------------
# Data
# -------------------------------------------------------------------------

# all samples, with `outlier` flag -- for the raw-data plots
dat <- load_wastewater() |>
  prepare_wastewater(pathogen = PATHOGEN, time_step = TIME_STEP) |>
  add_outlier_flag(window_days = OUTLIER_WINDOW_DAYS, factor = OUTLIER_FACTOR)

# without outliers -- for everything that is estimated
dat_clean <- dat |>
  filter(!outlier) |>
  add_time_index(TIME_STEP)

# -------------------------------------------------------------------------
# Fit models and summarise
# -------------------------------------------------------------------------

models <- fit_models(
  dat_clean,
  time_step = TIME_STEP, n_samples = N_SAMPLES, seed = SEED,
  weighted = USE_POPULATION_WEIGHTS, inla_models = INLA_MODELS,
  blq_handling = BLQ_HANDLING, zero_handling = ZERO_HANDLING
)
m_gam <- models$fits$GAM
inla_fits <- models$fits[names(INLA_MODELS)]

national_summary <- bind_rows(
  bind_rows(map(models$national, \(x) x$summary)),
  summarise_raw_data(dat_clean, window = RAW_SMOOTH_WINDOW, weighted = USE_POPULATION_WEIGHTS)
)

site_summary <- imap(models$site_draws, \(draws, model) summarise_site_draws(draws, models$grid, model)) |>
  bind_rows()

trends <- map(models$national, summarise_trend_draws,
              trend_days = TREND_DAYS, stable_threshold = TREND_STABLE_THRESHOLD) |>
  bind_rows()

site_trends <- imap(models$site_draws, \(draws, model) {
  summarise_site_trends(draws, models$grid, model, TREND_DAYS, TREND_STABLE_THRESHOLD)
}) |>
  bind_rows()

latest_trend_summary <- trends |>
  filter(date == max(date, na.rm = TRUE)) |>
  dplyr::select(model, date, estimate, lower, upper, p_increasing, p_stable, p_decreasing, direction)
latest_trend_summary

latest_site_trends(site_trends, dat_clean) |>
  filter(model %in% names(INLA_MODELS)) |>
  dplyr::select(model, wwtp, date, estimate, lower, upper, direction, days_since_last_obs) |>
  arrange(model, desc(estimate))

# -------------------------------------------------------------------------
# Diagnostics
# -------------------------------------------------------------------------

summary(m_gam)
gratia::draw(m_gam, residuals = TRUE)
gratia::appraise(m_gam)
gam.check(m_gam)

# zeros / extreme low values: sig2 should be close to the robust mad(r)^2
c(sig2 = m_gam$sig2, robust = mad(residuals(m_gam))^2)
dat |>
  group_by(wwtp) |>
  summarise(
    n = n(),
    n_zero = sum(viral_load == 0),
    n_blq = sum(blq),
    n_outlier = sum(outlier),
    n_loq_proxy = sum(loq_source == "proxy"),
    median_loq = median(loq),
    median_detection_floor = median(detection_floor)
  ) |>
  filter(n_zero > 0 | n_blq > 0 | n_outlier > 0)

# sanity check of column names / units: `value` recomputed from
# concentration x flow should match (ratio ~1)
dat |>
  filter(viral_load > 0) |>
  mutate(ratio = concentration_to_load(concentration, flow, population) / viral_load) |>
  summarise(n = sum(!is.na(ratio)), min = min(ratio, na.rm = TRUE),
            median = median(ratio, na.rm = TRUE), max = max(ratio, na.rm = TRUE))

map(inla_fits, summary)
map(inla_fits, \(m) m$summary.hyperpar)

# extra SD (log scale) of below-LOQ samples on top of the Gamma noise; ~0 =
# they behave like quantified samples
if (BLQ_HANDLING == "observed_extra") {
  imap(inla_fits, \(m, model) {
    prec <- m$summary.hyperpar["Precision for blq_obs", ]
    tibble(model = model, sd = 1 / sqrt(prec$`0.5quant`),
           lower = 1 / sqrt(prec$`0.975quant`), upper = 1 / sqrt(prec$`0.025quant`))
  }) |>
    bind_rows()
}

# lower WAIC = better fit (same data and likelihood for all INLA models)
imap(inla_fits, \(m, model) tibble(model = model, waic = m$waic$waic, p_eff = m$waic$p.eff)) |>
  bind_rows()

# site deviations must satisfy their constraints (range ~0):
#   AR1: sum over all plants is 0 in every time bin -- if not, constr = TRUE
#        did not apply per group
#   RW1: the explicit extraconstr rows (per-bin sums and per-plant means)
if ("INLA_AR1" %in% names(inla_fits)) {
  inla_site_deviations(inla_fits$INLA_AR1, levels(dat$wwtp)) |>
    group_by(time) |>
    summarise(total = sum(mean)) |>
    pull(total) |>
    range()
}
if ("INLA_RW1" %in% names(inla_fits)) {
  range(rw1_site_constraints(models$inla_data)$A %*% inla_fits$INLA_RW1$summary.random$wwtp_dev$mean)
}

plot_model_components(dat, inla_fits, gam_term_grid(m_gam, models$grid))
plot_site_fit(site_summary, dat, WWTP_SELECTION)
plot_site_deviations(inla_fits, models$grid, WWTP_SELECTION)

# -------------------------------------------------------------------------
# National signal and trend
# -------------------------------------------------------------------------

(plot_national(national_summary, weighted = USE_POPULATION_WEIGHTS)+scale_y_log10()) /
  plot_trend(trends, "interpolated difference")

# -------------------------------------------------------------------------
# WWTP trends
# -------------------------------------------------------------------------

plot_site_trend_directions(site_trends, since = max(dat$date) - 180)

plot_trend(filter(site_trends, wwtp %in% WWTP_SELECTION), "by WWTP") +
  facet_wrap(~ wwtp)

# -------------------------------------------------------------------------
# Nowcast backtest
# -------------------------------------------------------------------------

if (RUN_BACKTEST) {
  # Eyeballed off the plot -- swap in exact dates. Pick each cutoff a few
  # days *into* the phase so there's real signal for the model to react to.
  backtest_cutoffs <- tribble(
    ~phase,        ~cutoff_date,
    "wave_start",  as.Date("2025-10-15"),
    "peak",        as.Date("2025-12-01"),
    "declining",   as.Date("2026-02-15"),
    "low_season",  as.Date("2026-06-15")
  )
  nowcast_window_days <- 21L   # days back from each cutoff to score
  plot_window_days <- 100L     # days back from each cutoff to plot

  backtest_results <- run_backtest(
    dat_clean, backtest_cutoffs,
    time_step = TIME_STEP, n_samples = N_SAMPLES, seed = SEED,
    weighted = USE_POPULATION_WEIGHTS, inla_models = INLA_MODELS,
    blq_handling = BLQ_HANDLING, zero_handling = ZERO_HANDLING
  )

  final_estimates <- national_summary |>
    filter(model %in% names(models$national)) |>
    dplyr::select(date, model, final_estimate = estimate, final_lower = lower, final_upper = upper)

  backtest_results |>
    compare_backtest(final_estimates, nowcast_window_days) |>
    summarise_backtest() |>
    print()

  print(plot_backtest(backtest_results, final_estimates, backtest_cutoffs, plot_window_days))
}

# -------------------------------------------------------------------------
# Figures for the team presentation
# -------------------------------------------------------------------------

# Writes presentation/figures/<PATHOGEN>_<name>.png. Run once per pathogen
# (the page uses sars-cov-2 and influenza_a); with RUN_BACKTEST also the
# backtest figure.
if (EXPORT_FIGURES) {
  fig_dir <- fs::path(current_dir, "presentation", "figures")
  fs::dir_create(fig_dir)
  save_figure <- function(plot, name, width = 11, height = 6) {
    ggsave(fs::path(fig_dir, paste0(PATHOGEN, "_", name, ".png")), plot,
           width = width, height = height, dpi = 150, bg = "white")
  }

  save_figure(
    (plot_national(national_summary, weighted = USE_POPULATION_WEIGHTS) + scale_y_log10()) /
      plot_trend(trends, "interpolated difference"),
    "national", height = 8
  )
  save_figure(plot_site_fit(site_summary, dat, WWTP_SELECTION), "site_fit")
  save_figure(plot_site_deviations(inla_fits, models$grid, WWTP_SELECTION), "deviations")
  save_figure(plot_model_components(dat, inla_fits, gam_term_grid(m_gam, models$grid)),
              "components", width = 12, height = 9)
  if (RUN_BACKTEST) {
    save_figure(plot_backtest(backtest_results, final_estimates, backtest_cutoffs, plot_window_days),
                "backtest", width = 12, height = 6)
  }
}
