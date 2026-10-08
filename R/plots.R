# -------------------------------------------------------------------------
# Plots
# -------------------------------------------------------------------------

# GAM red, raw reference blue, INLA models in INLA_COLOURS (in label order).
INLA_COLOURS <- c("#388E3C", "#F57C00", "#7B1FA2", "#00897B")

model_colours <- function(models) {
  models <- sort(unique(models))
  inla <- grep("^INLA", models, value = TRUE)
  colours <- case_when(
    models == "GAM" ~ "#D32F2F",
    models %in% inla ~ INLA_COLOURS[match(models, inla)],
    TRUE ~ "#1976D2"
  )
  set_names(colours, models)
}

# ---- Model components -----------------------------------------------------

# All samples; outliers (excluded from the models) marked as red crosses.
plot_raw_loads <- function(dat) {
  ggplot(dat, aes(time, viral_load)) +
    geom_line(aes(colour = wwtp), alpha = 0.5) +
    geom_point(data = filter(dat, outlier), shape = 4, colour = "red", size = 2) +
    scale_y_log10() +
    theme_classic() +
    labs(title = "Raw viral loads by WWTPs", subtitle = "Red crosses: outliers, excluded from the models")
}

plot_inla_components <- function(m_inla, wwtp_levels, model_name) {
  national <- m_inla$summary.random$time_ch
  site <- inla_site_deviations(m_inla, wwtp_levels)

  ggplot() +
    geom_ribbon(data = national, aes(ID, ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.2) +
    geom_line(data = national, aes(ID, mean)) +
    geom_line(data = site, aes(time, mean, colour = wwtp), alpha = 0.5) +
    theme_classic() +
    labs(title = paste0(model_name, ": National RW2 vs WWTP-specific deviations (sum to zero per time bin)"),
         x = "Time bin")
}

plot_gam_components <- function(gam_terms) {
  ggplot() +
    geom_line(data = distinct(gam_terms, day, date, national_trend), aes(date, national_trend)) +
    geom_line(data = gam_terms, aes(date, site_deviation, colour = wwtp, group = wwtp), alpha = 0.5) +
    theme_classic() +
    labs(title = "GAM: National spline vs WWTP-specific splines")
}

# Persistent WWTP baseline (INLA iid term), on the log scale. The GAM has
# no separate baseline: its WWTP levels sit in the sz smooths.
plot_inla_baseline <- function(m_inla, wwtp_levels, model_name) {
  m_inla$summary.random$wwtp_id |>
    mutate(wwtp = wwtp_levels) |>
    ggplot(aes(wwtp, mean, colour = wwtp)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`)) +
    theme_classic() +
    labs(title = paste0(model_name, ": Persistent WWTP baseline"),
         x = NULL, y = "Offset from overall level (log)")
}

# inla_fits: named list of INLA fits, one row of panels each.
plot_model_components <- function(dat, inla_fits, gam_terms) {
  wwtp_levels <- levels(dat$wwtp)

  left <- c(
    list(plot_raw_loads(dat)),
    imap(inla_fits, \(m, model) plot_inla_components(m, wwtp_levels, model)),
    list(plot_gam_components(gam_terms))
  )
  right <- c(
    list(plot_spacer()),
    imap(inla_fits, \(m, model) plot_inla_baseline(m, wwtp_levels, model)),
    list(plot_spacer())
  )

  (wrap_plots(left, ncol = 1) | wrap_plots(right, ncol = 1)) +
    plot_layout(widths = c(2, 1), guides = "collect")
}

# WWTP deviations from the national trend (log scale, 0 = in line), median
# and 95% interval per time bin, for comparing INLA models. Only bins in the
# prediction grid (from a plant's first sample on) are shown.
plot_site_deviations <- function(inla_fits, grid, wwtps) {
  wwtp_levels <- levels(grid$wwtp)
  deviations <- imap(inla_fits, \(m, model) inla_site_deviations(m, wwtp_levels) |> mutate(model = model)) |>
    bind_rows() |>
    inner_join(distinct(grid, wwtp, time, date), by = c("wwtp", "time")) |>
    filter(wwtp %in% wwtps)

  ggplot(deviations, aes(date, `0.5quant`, colour = model, fill = model)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = model_colours(deviations$model)) +
    scale_fill_manual(values = model_colours(deviations$model)) +
    facet_wrap(~ wwtp) +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(x = NULL, y = "Deviation from national trend (log)", colour = "Model", fill = "Model",
         title = "WWTP deviations from the national trend",
         subtitle = "Posterior median and 95% interval; 0 = in line with the national trend")
}

# ---- Fit vs data ------------------------------------------------------------

plot_site_fit <- function(site_summary, dat, wwtps) {
  ggplot(filter(site_summary, wwtp %in% wwtps)) +
    geom_point(data = filter(dat, wwtp %in% wwtps, !outlier), aes(date, viral_load),
               size = 0.4, alpha = 0.8, colour = "grey40") +
    geom_point(data = filter(dat, wwtp %in% wwtps, outlier), aes(date, viral_load),
               shape = 4, colour = "red", size = 2) +
    geom_ribbon(aes(date, ymin = lower, ymax = upper, fill = model), alpha = 0.2, colour = NA) +
    geom_line(aes(date, estimate, colour = model), linewidth = 0.6) +
    scale_y_log10() +
    facet_wrap(~ wwtp, scales = "free_y") +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(title = "Model fit vs raw data, by WWTP", subtitle = "Red crosses: outliers, excluded from the models",
         x = NULL, y = "Viral load")
}

# ---- National signal and trend ----------------------------------------------

plot_national <- function(national_summary, weighted) {
  ggplot(national_summary) +
    geom_ribbon(aes(date, ymin = lower, ymax = upper, fill = model), alpha = 0.12, colour = NA) +
    geom_line(aes(date, estimate, colour = model), linewidth = 1, alpha = 0.6) +
    scale_colour_manual(values = model_colours(national_summary$model)) +
    scale_fill_manual(values = model_colours(national_summary$model)) +
    theme_classic() +
    labs(x = NULL, y = "Viral load", colour = "Model", fill = "Model",
         title = "Swiss wastewater signal",
         subtitle = paste(if (weighted) "Population-weighted" else "Unweighted",
                          "national signal by model"))
}

plot_trend <- function(trends, method_label) {
  trend_days <- first(trends$trend_days)
  stable_threshold <- first(trends$stable_threshold)

  ggplot(trends, aes(date, estimate, colour = model, fill = model)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 100 * c(-stable_threshold, stable_threshold), linetype = "dotted") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = model_colours(trends$model)) +
    scale_fill_manual(values = model_colours(trends$model)) +
    theme_classic() +
    labs(x = NULL, y = paste0("Change over ", trend_days, " days (%)"),
         colour = "Model", fill = "Model",
         title = paste0("Change in Swiss wastewater signal (", method_label, ")"),
         subtitle = paste0("Dotted lines show the ±", 100 * stable_threshold, "% stable range"))
}

# ---- WWTP trends --------------------------------------------------------------

DIRECTION_COLOURS <- c(increasing = "#D32F2F", stable = "#9E9E9E",
                       decreasing = "#1976D2", uncertain = "#EEEEEE")

# Trend direction per WWTP over time, one tile per model time bin.
plot_site_trend_directions <- function(site_trends, since) {
  site_trends |>
    filter(date >= since) |>
    ggplot(aes(date, wwtp, fill = direction)) +
    geom_tile() +
    scale_fill_manual(values = DIRECTION_COLOURS, na.value = "white") +
    facet_wrap(~ model) +
    theme_classic() +
    labs(x = NULL, y = NULL, fill = NULL,
         title = "Trend direction by WWTP",
         subtitle = paste0("Change over ", first(site_trends$trend_days), " days, ",
                           "±", 100 * first(site_trends$stable_threshold), "% stable range"))
}
