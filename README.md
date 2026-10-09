# MTDATA-5048: Approaches for a Swiss-wide wastewater average

Author: Moritz Wagner (moritz.wagner@bag.admin.ch)  
Date: 2026-09-28

## Question

How should per-WWTP wastewater viral loads be combined into a single Swiss-wide
signal (and trend), given that plants differ in baseline level, sampling frequency and
start of surveillance?

## Approaches compared

All produce a population-weighted (configurable) national series.

| Approach | Model | Site structure |
|---|---|---|
| **GAM** (`mgcv`) | Gaussian on `log1p(viral_load)`, back-transformed to the mean via `expm1(mu + sig2/2)` | shared smooth of time + WWTP-specific deviation smooths (`bs = "sz"`, which carry the WWTP levels) |
| **INLA_AR1** | Gamma (log link) on `viral_load`; non-detects left-censored at the LOQ (`gammasurv`), below-LOQ samples per `BLQ_HANDLING` | WWTP iid intercept + national RW2 + WWTP deviations (AR1 over weekly bins, sum to zero per bin); INLA default priors |
| **INLA_RW1** (optional, off by default) | as INLA_AR1 | as INLA_AR1, but the WWTP deviations are a RW1 over time bins (as in Riou et al., [multilevel_wbe](https://github.com/foph-modelling/multilevel_wbe)), summing to zero per bin and centred per WWTP |
| **Raw reference** | none | per-WWTP rolling median, then weighted national mean |

The INLA models to fit are set in `INLA_MODELS` in `analysis.R`; WAIC, the
constraint checks and a side-by-side plot of the WWTP deviations (log scale)
compare them.

For all models, 2000 posterior draws are taken at each WWTP × time-bin point, then
aggregated with population weights into national draws. Trends (% change over 21 days,
with P(increasing / stable / decreasing)) are computed from those draws, nationally and
per WWTP.

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

- Non-detects are left-censored in INLA (`ZERO_HANDLING`) at an empirical detection floor:
  no LODs are available, so per plant the 5th percentile of its positive concentrations
  (influenza A: ~3% of the LOQ), converted to load units with that day's flow and
  population. Positive samples below
  the LOQ (`BLQ_HANDLING`) are censored, used as measured, or used as measured with an
  estimated extra spread; the lab considers them usable. Censoring them left ~88% of
  influenza A samples censored and the low season unidentified. The GAM drops zeros;
  for sars-cov-2 the 45 zeros drove 82% of its residual variance.
- Outliers: quantified samples more than 10x away from the median of the plant's positive
  samples (quantified and below-LOQ) within ±7 days are excluded from all estimates (shown
  flagged in raw-data plots). To confirm
  with the lab that such values are artefacts.
- Backtest cutoff dates are eyeballed and still need to be pinned down.
