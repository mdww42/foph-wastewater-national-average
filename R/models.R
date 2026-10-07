# -------------------------------------------------------------------------
# Model fitting and site-level (WWTP x time) posterior draws
#
# Both predict_*_site() functions return a n_samples x nrow(grid) matrix of
# draws on the viral_load scale -- the shared building block for the
# national aggregate and the per-site fit-vs-raw comparison.
# -------------------------------------------------------------------------

# ---- GAM ----------------------------------------------------------------

# Gaussian GAM on log1p(viral_load).
fit_gam <- function(dat) {
  gam(
    y ~
      s(wwtp, bs = "re") +               # persistent WWTP-specific baseline
      s(day, bs = "tp", k = 80) +        # shared temporal component (k = 40: edf 38.6/39, k-index 0.93)
      s(day, wwtp, bs = "sz", k = 10),   # WWTP-specific deviations from it
    data = dat, method = "REML", select = TRUE
  )
}

predict_gam_site <- function(m_gam, grid, n_samples, seed) {
  Xp <- predict(m_gam, newdata = grid, type = "lpmatrix")

  set.seed(seed)
  beta_draws <- MASS::mvrnorm(n = n_samples, mu = coef(m_gam), Sigma = vcov(m_gam, unconditional = TRUE))
  # log1p(viral_load) ~ N(mu, sig2), so E[viral_load] = expm1(mu + sig2 / 2).
  # expm1(mu) alone would be the median, not the mean that INLA estimates.
  # sig2 is a point estimate, so its uncertainty is not propagated.
  expm1(beta_draws %*% t(Xp) + m_gam$sig2 / 2)
}

# GAM smooth terms evaluated on the grid, for the component / baseline plots.
gam_term_grid <- function(m_gam, grid) {
  terms <- predict(m_gam, newdata = grid, type = "terms", se.fit = TRUE)
  stopifnot(all(c("s(wwtp)", "s(day)", "s(day,wwtp)") %in% colnames(terms$fit)))

  grid |>
    mutate(
      site_baseline = terms$fit[, "s(wwtp)"],
      site_baseline_se = terms$se.fit[, "s(wwtp)"],
      national_trend = terms$fit[, "s(day)"],
      site_deviation = terms$fit[, "s(day,wwtp)"]
    )
}

# ---- INLA ---------------------------------------------------------------

# Observations plus the prediction grid (response NA) in one data set, so
# INLA returns the linear predictor at the grid points.
build_inla_data <- function(dat, grid) {
  bind_rows(
    dat |> mutate(is_prediction = FALSE),
    grid |>
      transmute(
        date, wwtp, viral_load = NA_real_, population, y = NA_real_,
        day, time, time_ch, time_wwtp, wwtp_id, is_prediction = TRUE
      )
  )
}

# Gamma (log link) model on viral_load:
#   WWTP iid baseline + national RW2 + WWTP-specific deviations.
#
# The site deviations are iid across WWTPs within each time bin, AR1 over
# time bins (`group`), and sum to zero within each bin (`constr = TRUE` on a
# grouped model applies per group). Without that constraint a wave common to
# all plants can sit in either the national term or in every plant's
# deviation; the fit put it in the deviations (site SD ~1, rho ~0.97, i.e.
# ~3-month memory) and left the national RW2 nearly linear.
#
# Weak PC priors, P(SD > u) = alpha, on the log scale:
#   baseline and site deviations: P(SD > 1) = 0.01
#   national RW2:                 P(SD > 1) = 0.5 (weak, lets it carry waves)
#
# censor = TRUE: non-detects enter as left-censored at `loq` (gammasurv);
# otherwise `inla_data` must not contain zeros.
fit_inla <- function(inla_data, censor) {
  inla_data <- inla_data |> mutate(wwtp_dev = wwtp_id)

  hyper_baseline <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  hyper_national <- list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
  hyper_site_dev <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

  latent <- ~ 1 +
    f(wwtp_id, model = "iid", constr = TRUE, hyper = hyper_baseline) +
    f(time_ch, model = "rw2", scale.model = TRUE, constr = TRUE, hyper = hyper_national) +
    f(wwtp_dev, model = "iid", constr = TRUE, hyper = hyper_site_dev,
      group = time_wwtp, control.group = list(model = "ar1"))

  if (censor) {
    # inla.surv event: 1 = observed, 2 = left-censored (true value < time).
    # Prediction rows have NA time, so they are predicted, not fitted.
    inla_data <- inla_data |>
      mutate(
        surv_time = if_else(censored %in% TRUE, loq, viral_load),
        surv_event = if_else(censored %in% TRUE, 2L, 1L)
      )
    formula_inla <- update(latent, inla.surv(surv_time, surv_event) ~ .)
    family <- "gammasurv"
  } else {
    formula_inla <- update(latent, viral_load ~ .)
    family <- "gamma"
  }

  inla(
    formula_inla, family = family, data = inla_data,
    control.predictor = list(compute = TRUE),
    control.compute = list(config = TRUE, waic = TRUE)
  )
}

# Posterior summary of the site deviations, one row per WWTP x time bin.
# The grouped term's latent vector runs over WWTPs within each time bin.
inla_site_deviations <- function(m_inla, wwtp_levels) {
  m_inla$summary.random$wwtp_dev |>
    as_tibble() |>
    mutate(
      wwtp = factor(wwtp_levels[ID], levels = wwtp_levels),
      time = (row_number() - 1L) %/% n_distinct(ID) + 1L
    )
}

predict_inla_site <- function(m_inla, inla_data, n_samples, seed) {
  set.seed(seed)
  samples <- inla.posterior.sample(n = n_samples, result = m_inla)

  latent_names <- rownames(samples[[1]]$latent)
  idx_predictor <- grep("^Predictor:", latent_names)
  predictor_rows <- as.integer(sub("^Predictor:", "", latent_names[idx_predictor]))

  pred_rows <- which(inla_data$is_prediction)
  use_idx <- match(pred_rows, predictor_rows)
  if (anyNA(use_idx)) stop("Could not match all prediction rows to INLA predictor latent indices.")

  eta_draws <- t(vapply(samples, function(s) s$latent[idx_predictor[use_idx], 1], numeric(length(pred_rows))))
  # log link -> exp(eta) is the mean viral_load
  exp(eta_draws)
}

# ---- Both models --------------------------------------------------------

# Fits GAM and INLA on `dat` and returns, per model, the fit, the site-level
# draws and the national aggregate. Used for the main analysis and for each
# backtest cutoff.
#
# zero_handling: "censor" fits INLA with non-detects left-censored at `loq`;
# "drop" fits INLA on positive loads only. The GAM always drops zeros (a
# censored GAM would need mgcv::cnorm()).
fit_models <- function(dat, time_step, n_samples, seed, weighted,
                       zero_handling = c("censor", "drop")) {
  zero_handling <- match.arg(zero_handling)
  grid <- build_site_time_grid(dat, time_step)
  dat_positive <- dat |> filter(viral_load > 0)

  m_gam <- fit_gam(dat_positive)

  censor <- zero_handling == "censor"
  inla_data <- build_inla_data(if (censor) dat else dat_positive, grid)
  m_inla <- fit_inla(inla_data, censor)

  site_draws <- list(
    GAM = predict_gam_site(m_gam, grid, n_samples, seed),
    INLA = predict_inla_site(m_inla, inla_data, n_samples, seed)
  )

  list(
    grid = grid,
    fits = list(GAM = m_gam, INLA = m_inla),
    site_draws = site_draws,
    national = imap(site_draws, \(draws, model) aggregate_national_draws(draws, grid, model, weighted))
  )
}
