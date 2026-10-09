# -------------------------------------------------------------------------
# Data loading and preparation
# -------------------------------------------------------------------------

WASTEWATER_URL <- "https://idd.bag.admin.ch/api/v1/export/current/RESPVIRUSES_wastewater/csv"

load_wastewater <- function(url = WASTEWATER_URL) {
  data.table::fread(url) |> tibble::as_tibble()
}

# One row per WWTP x sampling date for a single pathogen, with the
# modelling columns (log1p response, day and INLA time indices).
prepare_wastewater <- function(wastewater, pathogen, time_step) {
  wastewater |>
    filter(!is.na(value), georegion2 != "liechtenstein", valueCategory == pathogen) |>
    transmute(
      date = as.Date(temporal),
      wwtp = factor(georegion2),
      viral_load = value,
      population = pop,
      y = log1p(viral_load),
      # inputs for the LOQ in load units -- check the column names against the export
      concentration = conc,          # gene copies / L
      flow,                          # m3 / day
      loq_concentration = loq,       # gene copies / L
      below_loq = belowLoq           # concentration <= LOQ (FALSE where either is NA)
    ) |>
    add_loq() |>
    add_time_index(time_step)
}

# gene copies / L x m3 / day -> gene copies / day per 100,000 inhabitants,
# the unit of `value` in the export.
concentration_to_load <- function(concentration, flow, population) {
  concentration * flow * 1000 / (population / 1e5)
}

# The sample's LOQ in load units (lab LOQ converted with that day's flow and
# population), and `blq`: a positive load at or below the LOQ (zeros are
# non-detects, not `blq`). Where the lab LOQ or flow is missing, `loq` falls
# back to a per-plant proxy: the lowest positive load the plant has reported.
#
# `detection_floor`: the bound for non-detects, in load units. The lab
# reports positive concentrations well below the LOQ (influenza A down to
# ~2-3% of it at every plant), so a zero means "below what the lab detects",
# not "below the LOQ". Per plant, the `detection_quantile` quantile of its
# positive concentrations, converted with that day's flow and population.
#
# Which samples are censored is decided at fit time (see fit_models()).
add_loq <- function(dat, detection_quantile = 0.05) {
  dat |>
    group_by(wwtp) |>
    mutate(
      loq_proxy = min(viral_load[viral_load > 0]),
      loq_lab = concentration_to_load(loq_concentration, flow, population),
      loq_source = if_else(is.na(loq_lab), "proxy", "lab"),
      loq = coalesce(loq_lab, loq_proxy),
      detection_floor = concentration_to_load(
        quantile(concentration[concentration > 0], detection_quantile, na.rm = TRUE, names = FALSE),
        flow, population
      ),
      blq = below_loq %in% TRUE & viral_load > 0
    ) |>
    ungroup() |>
    dplyr::select(-loq_proxy, -loq_lab)
}

# Single samples far from their plant's neighbouring samples (e.g. chur at
# the end of 2023, ~30x its neighbours) dominate a time bin and with it the
# plant's trend. Flags quantified samples more than `factor` x away from the
# median of the plant's positive samples (quantified and below-LOQ) within
# +/- `window_days` (the sample itself included), given at least 3 such
# samples. Below-LOQ samples count as neighbours, so the rule also works in
# low seasons (influenza: mostly below-LOQ samples and zeros), but are not
# flagged themselves: they are noisier by definition. Zeros are neither.
add_outlier_flag <- function(dat, window_days, factor) {
  dat |>
    group_by(wwtp) |>
    arrange(date, .by_group = TRUE) |>
    mutate(
      log_load = if_else(viral_load > 0, log10(viral_load), NA_real_),
      neighbour_median = slider::slide_index_dbl(
        log_load, date, \(x) median(x, na.rm = TRUE),
        .before = window_days, .after = window_days
      ),
      n_neighbours = slider::slide_index_int(
        log_load, date, \(x) sum(!is.na(x)),
        .before = window_days, .after = window_days
      ),
      outlier = !blq & !is.na(log_load) & n_neighbours >= 3 &
        abs(log_load - neighbour_median) > log10(factor)
    ) |>
    ungroup() |>
    dplyr::select(-log_load, -neighbour_median, -n_neighbours)
}

# (Re)computes the day / time indices relative to the first date in `dat`.
# Separate from prepare_wastewater() so truncated data (backtest) can be
# re-indexed.
add_time_index <- function(dat, time_step) {
  dat |>
    mutate(
      day = as.integer(date - min(date)) + 1L,
      time = floor((day - 1L) / time_step) + 1L,
      time_ch = time,
      time_wwtp = time,
      wwtp_id = as.integer(wwtp)
    )
}
