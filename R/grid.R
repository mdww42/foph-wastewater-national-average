# -------------------------------------------------------------------------
# Prediction grid: WWTP x time points at which both models are evaluated
# -------------------------------------------------------------------------

get_site_info <- function(dat) {
  # last day/date are global, not site-specific, so every site's grid
  # extends to the same end point
  last_day <- max(dat$day)
  last_date <- max(dat$date)

  dat |>
    group_by(wwtp) |>
    summarise(
      population = median(population, na.rm = TRUE),
      first_day = min(day), last_day = last_day,
      first_date = min(date), last_date = last_date,
      .groups = "drop"
    )
}

# One prediction day per INLA time bin (the bin centre, capped at the last day).
get_prediction_days <- function(dat, time_step) {
  tibble(time = sort(unique(dat$time))) |>
    mutate(
      day = (time - 1L) * time_step + floor(time_step / 2) + 1L,
      day = pmin(day, max(dat$day)),
      date = min(dat$date) + day - 1L
    )
}

build_site_time_grid <- function(dat, time_step) {
  prediction_days <- get_prediction_days(dat, time_step)
  site_info <- get_site_info(dat)

  crossing(prediction_days, wwtp = levels(dat$wwtp)) |>
    mutate(
      wwtp = factor(wwtp, levels = levels(dat$wwtp)),
      wwtp_id = as.integer(wwtp),
      time_ch = time,
      time_wwtp = time
    ) |>
    left_join(site_info, by = "wwtp") |>
    # don't extrapolate a site's trajectory before it was surveilled
    filter(day >= first_day, day <= last_day) |>
    group_by(time) |>
    mutate(weight = population / sum(population)) |>
    ungroup()
}
