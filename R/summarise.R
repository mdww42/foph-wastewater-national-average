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

# Population-weighted (or unweighted) national series from site-level draws.
# Returns the per-day summary plus the national draws, which the trend
# summaries need.
aggregate_national_draws <- function(site_draws, grid, model_name, weighted = TRUE) {
  index <- grid |> distinct(day, date) |> arrange(day)

  national_draws <- vapply(index$day, function(d) {
    ii <- which(grid$day == d)
    w <- if (weighted) grid$weight[ii] else rep(1 / length(ii), length(ii))
    drop(site_draws[, ii, drop = FALSE] %*% w)
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

# Change over exactly `trend_days`, by interpolating each national draw on
# the model time grid. Worked out on the log1p scale so differences become
# proportional changes.
summarise_trend_draws <- function(national, trend_days, stable_threshold) {
  x <- national$index$day

  log_draws <- log1p(pmax(national$draws, 0))
  lagged_log_draws <- t(apply(log_draws, 1, function(z) {
    approx(x = x, y = z, xout = x - trend_days, rule = 1)$y
  }))
  log_change <- log_draws - lagged_log_draws

  increase_cutoff <- log1p(stable_threshold)
  decrease_cutoff <- log1p(-stable_threshold)

  national$index |>
    bind_cols(summarise_draw_columns(100 * expm1(log_change))) |>
    mutate(
      p_increasing = colMeans(log_change > increase_cutoff, na.rm = TRUE),
      p_stable = colMeans(log_change >= decrease_cutoff & log_change <= increase_cutoff, na.rm = TRUE),
      p_decreasing = colMeans(log_change < decrease_cutoff, na.rm = TRUE),
      direction = classify_direction(p_increasing, p_stable, p_decreasing),
      model = national$model, trend_days = trend_days, stable_threshold = stable_threshold
    )
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
