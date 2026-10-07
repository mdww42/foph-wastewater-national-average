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
      y = log1p(viral_load)
    ) |>
    add_loq() |>
    add_time_index(time_step)
}

# Zero loads are non-detects: the true load is somewhere below the sample's
# LOQ. The LOQ in load units is not in the export, so `loq` is a per-plant
# proxy: the lowest positive load that plant has reported. Replace with the
# per-sample lab LOQ (LOQ_conc x flow x 1000 / (pop / 1e5)) once available.
add_loq <- function(dat) {
  dat |>
    group_by(wwtp) |>
    mutate(
      censored = viral_load == 0,
      loq = min(viral_load[viral_load > 0])
    ) |>
    ungroup()
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
