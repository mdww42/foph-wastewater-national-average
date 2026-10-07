# -------------------------------------------------------------------------
# Nowcast backtest: GAM vs INLA
#
# For a handful of cutoff dates (wave start, peak, declining, low season),
# truncate the data to what would have been available at that point, refit
# both models and compare their nowcast for the preceding weeks against the
# "final" estimate from the full-data fit.
# -------------------------------------------------------------------------

# `...` is passed on to fit_models() (n_samples, seed, weighted, zero_handling).
refit_at_cutoff <- function(cutoff_date, dat, time_step, ...) {
  dat_trunc <- dat |>
    filter(date <= cutoff_date) |>
    add_time_index(time_step)

  fit_models(dat_trunc, time_step = time_step, ...)$national |>
    map(\(x) x$summary) |>
    bind_rows() |>
    mutate(cutoff_date = cutoff_date)
}

run_backtest <- function(dat, cutoffs, time_step, ...) {
  map(cutoffs$cutoff_date, refit_at_cutoff, dat = dat, time_step = time_step, ...) |>
    set_names(cutoffs$phase) |>
    bind_rows(.id = "phase")
}

compare_backtest <- function(backtest_results, final_estimates, nowcast_window_days) {
  backtest_results |>
    inner_join(final_estimates, by = c("date", "model")) |>
    filter(date >= cutoff_date - nowcast_window_days, date <= cutoff_date) |>
    mutate(
      days_before_cutoff = as.integer(cutoff_date - date),
      pct_error = 100 * (estimate - final_estimate) / final_estimate,
      covered = final_estimate >= lower & final_estimate <= upper
    )
}

summarise_backtest <- function(backtest_comparison) {
  backtest_comparison |>
    group_by(phase, model) |>
    summarise(
      mean_pct_error = mean(pct_error, na.rm = TRUE),  # bias: + = nowcast overshoots what was later confirmed
      mae_pct = mean(abs(pct_error), na.rm = TRUE),     # typical size of the miss
      coverage = mean(covered, na.rm = TRUE),           # should be ~0.95 if uncertainty is well-calibrated
      .groups = "drop"
    ) |>
    arrange(phase, model)
}

plot_backtest <- function(backtest_results, final_estimates, cutoffs, plot_window_days) {
  final_windowed <- pmap_dfr(cutoffs, function(phase, cutoff_date) {
    final_estimates |>
      filter(date >= cutoff_date - plot_window_days, date <= cutoff_date + 14) |>
      transmute(phase = phase, cutoff_date = cutoff_date, date, model,
                estimate = final_estimate, lower = final_lower, upper = final_upper,
                source = "final (all data)")
  })

  nowcasts <- backtest_results |>
    filter(date >= cutoff_date - plot_window_days) |>
    transmute(phase, cutoff_date, date, model, estimate, lower, upper,
              source = "nowcast (at cutoff)")

  bind_rows(final_windowed, nowcasts) |>
    mutate(phase = factor(phase, levels = cutoffs$phase)) |>
    ggplot(aes(date, estimate, colour = source, fill = source)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    geom_vline(aes(xintercept = cutoff_date), linetype = "dashed", colour = "grey40") +
    facet_grid(model ~ phase, scales = "free_x") +
    scale_y_log10() +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(x = NULL, y = "Viral load", colour = NULL, fill = NULL,
         title = "Nowcast at cutoff vs. final (fully revised) estimate",
         subtitle = "Dashed vertical line = cutoff date")
}
