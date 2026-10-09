# -------------------------------------------------------------------------
# Summaries of posterior draws: national aggregate, per-site fit, trends
# -------------------------------------------------------------------------

# Column-wise median and 95% interval of a draws matrix (samples x points).
summarise_draw_columns <- function(draws) {
  tibble(
    estimate = apply(draws, 2, median, na.rm = TRUE),
    lower = apply(draws, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(draws, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
}

# Population-weighted (or unweighted) national series from site-level draws,
# over the plants active on each day. Returns the per-day summary plus the
# national draws, which the trend summaries need.
#   "arithmetic": weighted mean of the loads (dominated by high-level plants)
#   "geometric":  weighted mean on the log scale (the typical plant)
aggregate_national_draws <- function(site_draws, grid, model_name, weighted = TRUE,
                                     mean = c("arithmetic", "geometric")) {
  mean <- match.arg(mean)
  index <- grid |> distinct(day, date) |> arrange(day)

  national_draws <- vapply(index$day, function(d) {
    ii <- which(grid$day == d)
    w <- if (weighted) grid$weight[ii] else rep(1 / length(ii), length(ii))
    x <- site_draws[, ii, drop = FALSE]
    if (mean == "geometric") drop(exp(log(x) %*% w)) else drop(x %*% w)
  }, numeric(nrow(site_draws)))

  list(
    summary = bind_cols(index, summarise_draw_columns(national_draws)) |> mutate(model = model_name),
    draws = national_draws,
    index = index,
    model = model_name
  )
}

# Site-level (WWTP x time) summary, for comparing fit against raw observations.
summarise_site_draws <- function(site_draws, grid, model_name) {
  grid |>
    dplyr::select(wwtp, day, date) |>
    bind_cols(summarise_draw_columns(site_draws)) |>
    mutate(model = model_name)
}

# ---- Trends -------------------------------------------------------------

classify_direction <- function(p_increasing, p_stable, p_decreasing, prob = 0.8) {
  case_when(
    p_increasing >= prob ~ "increasing",
    p_decreasing >= prob ~ "decreasing",
    p_stable >= prob ~ "stable",
    TRUE ~ "uncertain"
  )
}

# Matrix L with L %*% y = y (observed at `days`) linearly interpolated at
# `days - lag`; rows are NA where that falls before the first day.
lag_interpolation_matrix <- function(days, lag) {
  n <- length(days)
  if (n < 2) return(matrix(NA_real_, n, n))
  matrix(
    vapply(seq_len(n), function(k) {
      approx(days, as.numeric(seq_len(n) == k), xout = days - lag, rule = 1)$y
    }, numeric(n)),
    nrow = n
  )
}

# Change over exactly `trend_days` for each column of `draws` (samples x
# days), interpolating each draw on the model time grid. Worked out on the
# log1p scale so differences become proportional changes.
summarise_trend <- function(draws, days, trend_days, stable_threshold) {
  log_draws <- log1p(pmax(draws, 0))
  log_change <- log_draws - log_draws %*% t(lag_interpolation_matrix(days, trend_days))

  increase_cutoff <- log1p(stable_threshold)
  decrease_cutoff <- log1p(-stable_threshold)

  summarise_draw_columns(100 * expm1(log_change)) |>
    mutate(
      p_increasing = colMeans(log_change > increase_cutoff, na.rm = TRUE),
      p_stable = colMeans(log_change >= decrease_cutoff & log_change <= increase_cutoff, na.rm = TRUE),
      p_decreasing = colMeans(log_change < decrease_cutoff, na.rm = TRUE),
      direction = classify_direction(p_increasing, p_stable, p_decreasing),
      trend_days = trend_days, stable_threshold = stable_threshold
    )
}

summarise_trend_draws <- function(national, trend_days, stable_threshold) {
  national$index |>
    bind_cols(summarise_trend(national$draws, national$index$day, trend_days, stable_threshold)) |>
    mutate(model = national$model)
}

# Same trend summary per WWTP, from the site-level draws.
summarise_site_trends <- function(site_draws, grid, model_name, trend_days, stable_threshold) {
  grid |>
    mutate(draw_col = row_number()) |>
    group_by(wwtp) |>
    arrange(day, .by_group = TRUE) |>
    group_modify(\(site, key) {
      site |>
        dplyr::select(day, date) |>
        bind_cols(summarise_trend(site_draws[, site$draw_col, drop = FALSE], site$day,
                                  trend_days, stable_threshold))
    }) |>
    ungroup() |>
    mutate(model = model_name)
}

# Latest trend per WWTP, with how stale that site's data is: a site without
# recent samples mostly shows the national trend (an AR1 deviation reverts
# to zero, an RW1 deviation stays flat), not its own.
latest_site_trends <- function(site_trends, dat) {
  last_obs <- dat |>
    group_by(wwtp) |>
    summarise(last_obs_date = max(date), .groups = "drop")

  site_trends |>
    group_by(model, wwtp) |>
    slice_max(date, n = 1) |>
    ungroup() |>
    left_join(last_obs, by = "wwtp") |>
    mutate(days_since_last_obs = as.integer(max(dat$date) - last_obs_date))
}

# ---- Model-free reference -------------------------------------------------

# Per-WWTP rolling median over `window` days, then a (weighted) national mean.
summarise_raw_data <- function(dat, window, weighted = TRUE) {
  first_date <- min(dat$date)
  last_date <- max(dat$date)

  dat |>
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
      .groups = "drop"
    ) |>
    mutate(model = paste0("WWTP ", window, "-day rolling median + national mean"))
}
