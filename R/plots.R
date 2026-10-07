# -------------------------------------------------------------------------
# Plots
# -------------------------------------------------------------------------

# GAM, INLA, raw reference (in alphabetical order of the model labels)
MODEL_COLOURS <- c("#D32F2F", "#388E3C", "#1976D2")

# ---- Model components -----------------------------------------------------

plot_raw_loads <- function(dat) {
  ggplot(dat, aes(time, viral_load, colour = wwtp)) +
    geom_line(alpha = 0.5) +
    scale_y_log10() +
    theme_classic() +
    labs(title = "Raw viral loads by WWTPs")
}

plot_inla_components <- function(m_inla, wwtp_levels) {
  national <- m_inla$summary.random$time_ch
  site <- inla_site_deviations(m_inla, wwtp_levels)

  ggplot() +
    geom_ribbon(data = national, aes(ID, ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.2) +
    geom_line(data = national, aes(ID, mean)) +
    geom_line(data = site, aes(time, mean, colour = wwtp), alpha = 0.5) +
    theme_classic() +
    labs(title = "INLA: National RW2 vs WWTP-specific deviations (sum to zero per time bin)")
}

plot_gam_components <- function(gam_terms) {
  ggplot() +
    geom_line(data = distinct(gam_terms, day, date, national_trend), aes(date, national_trend)) +
    geom_line(data = gam_terms, aes(date, site_deviation, colour = wwtp, group = wwtp), alpha = 0.5) +
    theme_classic() +
    labs(title = "GAM: National spline vs WWTP-specific splines")
}

# Persistent WWTP baseline, GAM (re) vs INLA (iid).
# NOTE: not on the same numeric scale -- GAM's is on the log1p(viral_load)
# identity-link scale, INLA's on the gamma log-link scale. Comparable in
# spirit (both ~ log-fold offset from the overall level), not in units.
site_baselines <- function(gam_terms, m_inla, wwtp_levels) {
  gam_baseline <- gam_terms |>
    distinct(wwtp, site_baseline, site_baseline_se) |>
    transmute(
      wwtp, estimate = site_baseline,
      lower = site_baseline - 1.96 * site_baseline_se,
      upper = site_baseline + 1.96 * site_baseline_se,
      model = "GAM"
    )

  inla_baseline <- m_inla$summary.random$wwtp_id |>
    mutate(wwtp = wwtp_levels) |>
    transmute(wwtp, estimate = mean, lower = `0.025quant`, upper = `0.975quant`, model = "INLA")

  list(INLA = inla_baseline, GAM = gam_baseline)
}

plot_site_baseline <- function(baseline) {
  ggplot(baseline, aes(wwtp, estimate, colour = wwtp)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_pointrange(aes(ymin = lower, ymax = upper)) +
    theme_classic() +
    labs(title = paste0(first(baseline$model), ": Persistent WWTP baseline"),
         x = NULL, y = "Offset from overall level")
}

plot_model_components <- function(dat, m_inla, gam_terms) {
  wwtp_levels <- levels(dat$wwtp)
  baselines <- site_baselines(gam_terms, m_inla, wwtp_levels) |> map(plot_site_baseline)

  left <- plot_raw_loads(dat) / plot_inla_components(m_inla, wwtp_levels) / plot_gam_components(gam_terms)
  right <- plot_spacer() / baselines$INLA / baselines$GAM

  (left | right) + plot_layout(widths = c(2, 1), guides = "collect")
}

# ---- Fit vs data ------------------------------------------------------------

plot_site_fit <- function(site_summary, dat, wwtps) {
  ggplot(filter(site_summary, wwtp %in% wwtps)) +
    geom_point(data = filter(dat, wwtp %in% wwtps), aes(date, viral_load),
               size = 0.4, alpha = 0.8, colour = "grey40") +
    geom_ribbon(aes(date, ymin = lower, ymax = upper, fill = model), alpha = 0.2, colour = NA) +
    geom_line(aes(date, estimate, colour = model), linewidth = 0.6) +
    scale_y_log10() +
    facet_wrap(~ wwtp, scales = "free_y") +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(title = "Model fit vs raw data, by WWTP", x = NULL, y = "Viral load")
}

# ---- National signal and trend ----------------------------------------------

plot_national <- function(national_summary, weighted) {
  ggplot(national_summary) +
    geom_ribbon(aes(date, ymin = lower, ymax = upper, fill = model), alpha = 0.12, colour = NA) +
    geom_line(aes(date, estimate, colour = model), linewidth = 1, alpha = 0.6) +
    scale_colour_manual(values = MODEL_COLOURS) +
    scale_fill_manual(values = MODEL_COLOURS) +
    theme_classic() +
    labs(x = NULL, y = "Viral load", colour = "Model", fill = "Model",
         title = "Swiss wastewater signal",
         subtitle = paste(if (weighted) "Population-weighted" else "Unweighted",
                          "national signal: GAM vs INLA"))
}

plot_trend <- function(trends, method_label) {
  trend_days <- first(trends$trend_days)
  stable_threshold <- first(trends$stable_threshold)

  ggplot(trends, aes(date, estimate, colour = model, fill = model)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 100 * c(-stable_threshold, stable_threshold), linetype = "dotted") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = MODEL_COLOURS) +
    scale_fill_manual(values = MODEL_COLOURS) +
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
