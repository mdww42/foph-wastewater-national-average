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

# Fold-change axis (x0.5 ... x4) on a log scale, so rises and falls of the
# same size look the same and influenza's x10+ season starts don't flatten
# everything else.
FOLD_BREAKS <- c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16)
fold_labels <- function(x) paste0("×", format(x, drop0trailing = TRUE))

plot_trend <- function(trends, method_label) {
  trend_days <- first(trends$trend_days)
  stable_threshold <- first(trends$stable_threshold)

  trends |>
    mutate(across(c(estimate, lower, upper), \(x) 1 + x / 100)) |>
    ggplot(aes(date, estimate, colour = model, fill = model)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 1 + c(-stable_threshold, stable_threshold), linetype = "dotted") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 1) +
    scale_y_log10(breaks = FOLD_BREAKS, labels = fold_labels) +
    scale_colour_manual(values = model_colours(trends$model)) +
    scale_fill_manual(values = model_colours(trends$model)) +
    theme_classic() +
    labs(x = NULL, y = paste0("Change over ", trend_days, " days"),
         colour = "Model", fill = "Model",
         title = paste0("Change in Swiss wastewater signal (", method_label, ")"),
         subtitle = paste0("Dotted lines show the ±", 100 * stable_threshold, "% stable range"))
}

# What the dashboard would show: one model's national signal over the last
# `months` months, with the 21-day trend label per time bin underneath.
plot_national_dashboard <- function(national_summary, trends, model_name, months = 12) {
  since <- max(national_summary$date) - round(months * 30.44)
  signal <- filter(national_summary, model == model_name, date >= since)
  labels <- filter(trends, model == model_name, date >= since, !is.na(direction))
  bin_days <- as.numeric(median(diff(sort(unique(labels$date)))))
  x_scale <- scale_x_date(limits = range(signal$date) + c(-bin_days, bin_days) / 2, expand = c(0, 0))

  top <- ggplot(signal, aes(date, estimate)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = model_colours(model_name), alpha = 0.2) +
    geom_line(colour = model_colours(model_name), linewidth = 1) +
    scale_y_log10() +
    x_scale +
    theme_classic() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = NULL, y = "Viral load (gc/day/100k)",
         title = "Swiss wastewater signal",
         subtitle = paste0("Population-weighted national signal (", model_name, "), median and 95% interval"))

  strip <- ggplot(labels, aes(date, 1, fill = direction)) +
    geom_tile(width = bin_days, height = 1) +
    scale_fill_manual(values = DIRECTION_COLOURS, drop = FALSE) +
    x_scale +
    theme_classic() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          legend.position = "bottom") +
    labs(x = NULL, y = NULL, fill = paste0(first(labels$trend_days), "-day trend"))

  top / strip + plot_layout(heights = c(5, 1))
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

# ---- WWTPs compared with the national trend ---------------------------------

STATUS_COLOURS <- c(above = "#D95F02", "in line" = "#BDBDBD", below = "#7570B3",
                    uncertain = "#F0F0F0", "no recent data" = "#FFFFFF")

status_subtitle <- function(site_status) {
  threshold <- first(site_status$threshold)
  paste0("Deviation from the population-weighted national trend, after accounting\n",
         "for each plant's usual level; mean over the last ", first(site_status$window_bins), " time bins.\n",
         "Above / below: >= 80% probability beyond ×", format(threshold), " / ×",
         format(round(1 / threshold, 2)), "; in line: >= 80% within.")
}

# Every plant over time: the status label (`fill = "status"`) or the median
# deviation itself (`fill = "deviation"`, x0.5 to x2).
plot_site_status_heatmap <- function(site_status, since = min(site_status$date), fill = c("status", "deviation")) {
  fill <- match.arg(fill)
  d <- filter(site_status, date >= since, !is.na(status))
  bin_days <- as.numeric(median(diff(sort(unique(d$date)))))

  p <- ggplot(d, aes(date, wwtp)) +
    facet_wrap(~ model) +
    theme_classic() +
    labs(x = NULL, y = NULL, title = "WWTPs compared with the national trend", subtitle = status_subtitle(d))
  if (fill == "status") {
    p + geom_tile(aes(fill = factor(status, levels = names(STATUS_COLOURS))), width = bin_days) +
      scale_fill_manual(values = STATUS_COLOURS, drop = FALSE) +
      guides(fill = guide_legend(override.aes = list(colour = "grey70"))) +
      labs(fill = NULL)
  } else {
    p + geom_tile(aes(fill = log(estimate)), width = bin_days) +
      scale_fill_gradient2(low = STATUS_COLOURS[["below"]], mid = "#F5F5F5", high = STATUS_COLOURS[["above"]],
                           midpoint = 0, limits = log(c(0.5, 2)), oob = scales::squish,
                           breaks = log(c(0.5, 0.8, 1, 1.25, 2)), labels = fold_labels(c(0.5, 0.8, 1, 1.25, 2))) +
      labs(fill = "vs national")
  }
}

# Latest status per plant: median and 95% interval of the deviation, sorted,
# with the "in line" band shaded. Plants without recent samples are labelled.
plot_site_status_latest <- function(site_status) {
  latest <- site_status |>
    group_by(model) |>
    filter(date == max(date)) |>
    ungroup() |>
    mutate(wwtp = reorder(wwtp, estimate),
           note = if_else(status == "no recent data", paste0("no sample for ", days_since_last_obs, " days"), NA_character_))
  threshold <- first(latest$threshold)

  ggplot(latest, aes(estimate, wwtp)) +
    annotate("rect", xmin = 1 / threshold, xmax = threshold, ymin = -Inf, ymax = Inf, fill = "grey92") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_errorbar(aes(xmin = lower, xmax = upper, colour = status), width = 0, orientation = "y", linewidth = 0.8) +
    geom_point(aes(colour = status), size = 2.5) +
    geom_text(aes(x = upper, label = note), hjust = -0.1, size = 3, colour = "grey40", na.rm = TRUE) +
    scale_x_log10(breaks = sort(c(FOLD_BREAKS, 1 / threshold, threshold)), labels = fold_labels) +
    # darker greys than the heatmap, so points stay visible
    scale_colour_manual(values = c(STATUS_COLOURS[c("above", "below")], "in line" = "#616161",
                                   uncertain = "#9E9E9E", "no recent data" = "#BDBDBD"), drop = FALSE) +
    facet_wrap(~ model) +
    theme_classic() +
    labs(x = "Compared with the national trend", y = NULL, colour = NULL,
         title = paste("WWTPs compared with the national trend,", format(max(latest$date), "%d %b %Y")),
         subtitle = status_subtitle(latest))
}
