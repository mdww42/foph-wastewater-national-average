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

# Non-detects and samples at or below the LOQ are censored: the true load is
# somewhere below the sample's LOQ, converted to load units with that day's
# flow and population. Where the lab LOQ or flow is missing, `loq` falls back
# to a per-plant proxy: the lowest positive load the plant has reported.
add_loq <- function(dat) {
  dat |>
    group_by(wwtp) |>
    mutate(
      loq_proxy = min(viral_load[viral_load > 0]),
      loq_lab = concentration_to_load(loq_concentration, flow, population),
      loq_source = if_else(is.na(loq_lab), "proxy", "lab"),
      loq = coalesce(loq_lab, loq_proxy),
      censored = below_loq %in% TRUE | viral_load == 0
    ) |>
    ungroup() |>
    dplyr::select(-loq_proxy, -loq_lab)
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
