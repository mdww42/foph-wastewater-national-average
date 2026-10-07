# ============================================================
# Analysis: Swiss-wide wastewater average -- GAM vs INLA
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
TIME_STEP <- 3L                   # days per INLA time bin (also the prediction grid)
USE_POPULATION_WEIGHTS <- TRUE
TREND_DAYS <- 14L
TREND_STABLE_THRESHOLD <- 0.10    # +/- 10% over TREND_DAYS counts as "stable"
RAW_SMOOTH_WINDOW <- 7L           # days, for the model-free reference

# Zero loads are non-detects (below the LOQ). On the log scale a zero sits
# ~30 units below typical values: for sars-cov-2, 45 zeros caused 82% of
# the GAM residual variance (sig2 3.25 vs robust 0.34), inflating its mean
# correction ~5x and pulling the fitted curves down.
#   "censor": INLA treats them, and positive samples at or below the LOQ,
#             as left-censored at the LOQ in load units (see add_loq()).
#             Needed for influenza / RSV.
#   "drop":   INLA ignores them; biases low periods upwards.
# The GAM always drops zeros; the raw reference always uses all samples.
ZERO_HANDLING <- "censor"

WWTP_SELECTION <- c("aire", "werdhoelzli", "chur", "laupen")  # for the fit-vs-data plot

RUN_BACKTEST <- FALSE             # refits both models per cutoff -- slow

# -------------------------------------------------------------------------
# Data
# -------------------------------------------------------------------------

dat <- load_wastewater() |>
  prepare_wastewater(pathogen = PATHOGEN, time_step = TIME_STEP)

# -------------------------------------------------------------------------
# Fit models and summarise
# -------------------------------------------------------------------------

models <- fit_models(
  dat,
  time_step = TIME_STEP, n_samples = N_SAMPLES, seed = SEED,
  weighted = USE_POPULATION_WEIGHTS, zero_handling = ZERO_HANDLING
)
m_gam <- models$fits$GAM
m_inla <- models$fits$INLA

national_summary <- bind_rows(
  bind_rows(map(models$national, \(x) x$summary)),
  summarise_raw_data(dat, window = RAW_SMOOTH_WINDOW, weighted = USE_POPULATION_WEIGHTS)
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

latest_site_trends(site_trends, dat) |>
  filter(model == "INLA") |>
  dplyr::select(wwtp, date, estimate, lower, upper, direction, days_since_last_obs) |>
  arrange(desc(estimate))

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
    n_censored = sum(censored),
    n_loq_proxy = sum(loq_source == "proxy"),
    median_loq = median(loq)
  ) |>
  filter(n_censored > 0)

# sanity check of column names / units: `value` recomputed from
# concentration x flow should match (ratio ~1)
dat |>
  filter(viral_load > 0) |>
  mutate(ratio = concentration_to_load(concentration, flow, population) / viral_load) |>
  summarise(n = sum(!is.na(ratio)), min = min(ratio, na.rm = TRUE),
            median = median(ratio, na.rm = TRUE), max = max(ratio, na.rm = TRUE))

summary(m_inla)
m_inla$summary.hyperpar

plot_model_components(dat, m_inla, gam_term_grid(m_gam, models$grid))
plot_site_fit(site_summary, dat, WWTP_SELECTION)

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
    dat, backtest_cutoffs,
    time_step = TIME_STEP, n_samples = N_SAMPLES, seed = SEED,
    weighted = USE_POPULATION_WEIGHTS, zero_handling = ZERO_HANDLING
  )

  final_estimates <- national_summary |>
    filter(model %in% c("GAM", "INLA")) |>
    dplyr::select(date, model, final_estimate = estimate, final_lower = lower, final_upper = upper)

  backtest_results |>
    compare_backtest(final_estimates, nowcast_window_days) |>
    summarise_backtest() |>
    print()

  print(plot_backtest(backtest_results, final_estimates, backtest_cutoffs, plot_window_days))
}
