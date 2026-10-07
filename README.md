# MTDATA-5048: Approaches for a Swiss-wide wastewater average

Author: Moritz Wagner (moritz.wagner@bag.admin.ch)  
Date: 2026-09-28

## Question

How should per-WWTP wastewater viral loads be combined into a single Swiss-wide
signal (and trend), given that plants differ in baseline level, sampling frequency and
start of surveillance?

## Approaches compared

All three produce a population-weighted (configurable) national series.

| Approach | Model | Site structure |
|---|---|---|
| **GAM** (`mgcv`) | Gaussian on `log1p(viral_load)`, back-transformed to the mean via `expm1(mu + sig2/2)` | WWTP random intercept + shared smooth of time + WWTP-specific deviation smooths (`bs = "sz"`) |
| **INLA** | Gamma (log link) on `viral_load`; non-detects left-censored (`gammasurv`) | WWTP iid intercept + national RW2 + WWTP deviations (AR1 over 3-day bins, sum to zero per bin); weak PC priors |
| **Raw reference** | none | per-WWTP rolling median, then weighted national mean |

For both models, 2000 posterior draws are taken at each WWTP × time-bin point, then
aggregated with population weights into national draws. Trends (% change over 14 days,
with P(increasing / stable / decreasing)) come from those draws in two ways:
interpolated point-to-point difference and an OLS slope over the window.

A nowcast backtest (refit on data truncated at a few cutoff dates, compared with the
full-data estimate) is included but off by default (`RUN_BACKTEST <- FALSE`).

## Structure

```
analysis.R      main script: config, data, fits, diagnostics, plots, backtest
R/data.R        load and prepare the IDD wastewater export
R/grid.R        WWTP x time prediction grid and population weights
R/models.R      GAM / INLA fits, site-level posterior draws, fit_models() wrapper
R/summarise.R   national aggregation, site summaries, trend summaries, raw reference
R/plots.R       component, fit-vs-data, national and trend plots
R/backtest.R    nowcast backtest
```

## Running

Open `analysis.R` with this folder as working directory and run it top to bottom.
Pathogen and other settings are in the Config block. Data is fetched from the public
IDD export (`RESPVIRUSES_wastewater`).

Required packages: dplyr, tidyr, purrr, ggplot2, patchwork, mgcv, INLA, MASS, gratia,
slider, data.table, tibble and bag.epi.stats.

## Open points

- Non-detects and samples at or below the LOQ are left-censored in INLA (`ZERO_HANDLING`)
  at the lab LOQ converted to load units with that day's flow and population. Where LOQ
  or flow is missing, a per-plant proxy (lowest positive load) is used. The GAM drops
  zeros; for sars-cov-2 the 45 zeros drove 82% of its residual variance.
- The GAM and INLA WWTP baselines are on different scales (log1p identity vs. gamma log
  link), so they compare only qualitatively.
- Backtest cutoff dates are eyeballed and still need to be pinned down.
