# -------------------------------------------------------------------------
# Model fitting and site-level (WWTP x time) posterior draws
#
# Both predict_*_site() functions return a n_samples x nrow(grid) matrix of
# draws on the viral_load scale -- the shared building block for the
# national aggregate and the per-site fit-vs-raw comparison.
# -------------------------------------------------------------------------

# ---- GAM ----------------------------------------------------------------

# Gaussian GAM on log1p(viral_load). The sz smooths include each WWTP's
# level, so a separate s(wwtp, bs = "re") baseline was redundant (edf ~0).
fit_gam <- function(dat) {
  gam(
    y ~
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

# GAM smooth terms evaluated on the grid, for the component plot.
gam_term_grid <- function(m_gam, grid) {
  terms <- predict(m_gam, newdata = grid, type = "terms")
  stopifnot(all(c("s(day)", "s(day,wwtp)") %in% colnames(terms)))

  grid |>
    mutate(
      national_trend = terms[, "s(day)"],
      site_deviation = terms[, "s(day,wwtp)"]
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
# The site deviations are iid across WWTPs within each time bin and sum to
# zero within each bin. Without that constraint a wave common to all plants
# can sit in either the national term or in every plant's deviation; the fit
# put it in the deviations (site SD ~1, rho ~0.97, i.e. ~3-month memory) and
# left the national RW2 nearly linear. With it, the national RW2 is the
# plants' average on the log scale.
#
# site_dev sets how a plant's deviation evolves over time bins (`group`):
#   "ar1": stationary, pulled back to 0 (~3-4 week memory whatever the prior:
#          PC priors pushing it slower were overridden by the data, ~1000
#          samples per plant). Reverts to 0 after a plant's last sample.
#          The per-bin sum to zero is `constr = TRUE` (applies per group).
#   "rw1": random walk, as in Riou et al. (multilevel_wbe). A deviation can
#          persist for weeks and is carried forward after the last sample.
#          The RW1 leaves each plant's level free, so the constraints are set
#          explicitly (see rw1_site_constraints()). They span plants and time
#          bins, which `extraconstr` on a grouped term cannot express (it is
#          applied within each group), so the term is a generic0 with the
#          precision of a scaled RW1 over bins x iid over plants.
# The paper's own term (RW1 per plant, centred over all bins, no sum to zero
# across plants) was also tried: same deviations, but much wider intervals,
# as the national RW2 and the plants' common shift are not identified.
# INLA's default priors are used. The baseline is estimated from only 10
# plants and is the one term where the prior can matter.
#
# censor = TRUE: rows with `censored` enter as left-censored at `censor_at`
# (gammasurv); otherwise `inla_data` must not contain zeros.
# blq_extra = TRUE: observed below-LOQ samples get an extra per-sample iid
# effect, i.e. more spread than the Gamma noise; its PC prior shrinks the
# extra SD towards 0, so ~0 means they behave like quantified samples.
#
# The sum to zero gives every plant the same weight. Population weights were
# tried and dropped: the national RW2 then follows the large plants, so what
# the small plants share (e.g. a slow common drift) has to sit in their
# independent deviations; for sars-cov-2 rho went to ~0.99 and the plant
# levels moved into the deviations. Population weighting is done when
# aggregating the plant curves instead (aggregate_national_draws()).
fit_inla <- function(inla_data, censor, site_dev = c("ar1", "rw1"), blq_extra = FALSE) {
  site_dev <- match.arg(site_dev)
  inla_data <- inla_data |> mutate(wwtp_dev = wwtp_id)

  if (site_dev == "ar1") {
    latent <- ~ 1 +
      f(wwtp_id, model = "iid", constr = TRUE) +
      f(time_ch, model = "rw2", scale.model = TRUE, constr = TRUE) +
      f(wwtp_dev, model = "iid", constr = TRUE,
        group = time_wwtp, control.group = list(model = "ar1"))
  } else {
    n_wwtp <- max(inla_data$wwtp_id)
    n_time <- max(inla_data$time_wwtp)
    inla_data <- inla_data |> mutate(wwtp_dev = (time_wwtp - 1L) * n_wwtp + wwtp_id)
    site_precision <- rw1_site_precision(n_wwtp, n_time)
    site_constr <- rw1_site_constraints(inla_data)
    latent <- ~ 1 +
      f(wwtp_id, model = "iid", constr = TRUE) +
      f(time_ch, model = "rw2", scale.model = TRUE, constr = TRUE) +
      # rankdef = number of constraints: INLA does not correct the prior's
      # normalising constant for `extraconstr`, and with rankdef = n_wwtp the
      # missing ~(n_constr - n_wwtp) / 2 * log(tau) pushed the precision to
      # ~2e6, i.e. switched the deviations off. The constraints include the
      # RW1's null space (each plant's level).
      f(wwtp_dev, model = "generic0", Cmatrix = site_precision, rankdef = nrow(site_constr$A),
        constr = FALSE, extraconstr = site_constr)
  }

  if (blq_extra) {
    inla_data <- inla_data |>
      mutate(blq_obs = if_else(blq %in% TRUE & !(censored %in% TRUE), row_number(), NA_integer_))
    latent <- update(latent, ~ . + f(blq_obs, model = "iid",
                                     hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))))
  }

  if (censor) {
    # inla.surv event: 1 = observed, 2 = left-censored (true value < time).
    # Prediction rows have NA time, so they are predicted, not fitted.
    inla_data <- inla_data |>
      mutate(
        surv_time = if_else(censored %in% TRUE, censor_at, viral_load),
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
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(config = TRUE, waic = TRUE)
  )
}

# Precision of the RW1 site deviations: RW1 over time bins (scaled to
# generalised variance 1, like `scale.model = TRUE`) x iid over plants. The
# Kronecker order puts plants within time bins, as in the grouped AR1 term.
rw1_site_precision <- function(n_wwtp, n_time) {
  # first differences, (n_time - 1) x n_time
  D <- Matrix::sparseMatrix(
    i = rep(seq_len(n_time - 1L), 2), j = c(seq_len(n_time - 1L), seq_len(n_time - 1L) + 1L),
    x = rep(c(-1, 1), each = n_time - 1L), dims = c(n_time - 1L, n_time)
  )
  Q_time <- inla.scale.model(Matrix::crossprod(D), constr = list(A = matrix(1, 1, n_time), e = 0))
  Matrix::kronecker(Q_time, Matrix::Diagonal(n_wwtp))
}

# Linear constraints (A x = 0) on the RW1 site deviations. Their latent
# vector runs over WWTPs within each time bin, so plant i in bin t is entry
# (t - 1) * n_wwtp + i. Only a plant's active bins (rows in `inla_data`: its
# first sample to the end of the grid) are constrained:
#   - per bin, the active plants' deviations sum to zero (the AR1 version's
#     `constr = TRUE`, minus plants not yet surveilled);
#   - per plant, its deviations average to zero over its active bins, as in
#     Riou et al. The RW1 does not penalise a plant's overall level, which
#     would otherwise trade off against the iid baseline.
# Both sets add up every active cell once, so one plant row is redundant and
# is dropped.
rw1_site_constraints <- function(inla_data) {
  n_wwtp <- max(inla_data$wwtp_id)
  n_time <- max(inla_data$time_wwtp)
  active <- inla_data |> distinct(wwtp_id, time_wwtp)
  index <- (active$time_wwtp - 1L) * n_wwtp + active$wwtp_id

  by_wwtp <- split(index, active$wwtp_id)
  rows <- c(split(index, active$time_wwtp), by_wwtp[-length(by_wwtp)])

  A <- matrix(0, nrow = length(rows), ncol = n_wwtp * n_time)
  for (k in seq_along(rows)) A[k, rows[[k]]] <- 1
  if (qr(A)$rank < nrow(A)) stop("RW1 site constraints are linearly dependent.")

  list(A = A, e = rep(0, nrow(A)))
}

# Posterior summary of the site deviations, one row per WWTP x time bin.
# In both versions the latent vector runs over WWTPs within each time bin.
inla_site_deviations <- function(m_inla, wwtp_levels) {
  n_wwtp <- length(wwtp_levels)
  m_inla$summary.random$wwtp_dev |>
    as_tibble() |>
    mutate(
      wwtp = factor(wwtp_levels[(row_number() - 1L) %% n_wwtp + 1L], levels = wwtp_levels),
      time = (row_number() - 1L) %/% n_wwtp + 1L
    )
}

predict_inla_site <- function(samples, inla_data) {
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

# Posterior draws of the plant deviations d (log scale) at the grid's plant x
# time-bin points, n_samples x nrow(grid). In both versions the latent vector
# runs over plants within time bins; for AR1 the bins are the group index
# (position among the time bins in the data), for RW1 the time bin itself.
inla_deviation_draws <- function(samples, grid, inla_data, site_dev) {
  latent_names <- rownames(samples[[1]]$latent)
  n_wwtp <- max(inla_data$wwtp_id)
  time_pos <- if (site_dev == "rw1") grid$time else match(grid$time, sort(unique(inla_data$time_wwtp)))
  rows <- match(paste0("wwtp_dev:", (time_pos - 1L) * n_wwtp + grid$wwtp_id), latent_names)
  if (anyNA(rows)) stop("Could not match all grid points to the deviation term in the posterior samples.")
  t(vapply(samples, function(s) s$latent[rows, 1], numeric(nrow(grid))))
}

# ---- All models ---------------------------------------------------------

# Fits the GAM and one INLA model per entry of `inla_models` on `dat` and
# returns, per model, the fit, the site-level draws and the national
# aggregate. Used for the main analysis and for each backtest cutoff.
#
# inla_models: named vector of fit_inla() `site_dev` values; the names are
# the model labels, e.g. c(INLA_AR1 = "ar1", INLA_RW1 = "rw1").
# blq_handling, positive samples at or below the LOQ in INLA:
#   "censored":       left-censored at `loq`
#   "observed":       used as measured
#   "observed_extra": used as measured, with extra estimated spread
# zero_handling, non-detects in INLA:
#   "censor": left-censored at the plant's `detection_floor`
#   "drop":   removed
# The GAM always uses below-LOQ samples as measured and drops zeros (a
# censored GAM would need mgcv::cnorm()).
fit_models <- function(dat, time_step, n_samples, seed, weighted, inla_models,
                       blq_handling = c("censored", "observed", "observed_extra"),
                       zero_handling = c("censor", "drop")) {
  blq_handling <- match.arg(blq_handling)
  zero_handling <- match.arg(zero_handling)
  grid <- build_site_time_grid(dat, time_step)
  dat_positive <- dat |> filter(viral_load > 0)

  dat_inla <- (if (zero_handling == "drop") dat_positive else dat) |>
    mutate(
      censored = viral_load == 0 | (blq & blq_handling == "censored"),
      censor_at = if_else(viral_load == 0, detection_floor, loq)
    )
  censor <- any(dat_inla$censored)
  inla_data <- build_inla_data(dat_inla, grid)
  blq_extra <- blq_handling == "observed_extra"

  fits <- c(
    list(GAM = fit_gam(dat_positive)),
    map(inla_models, \(site_dev) fit_inla(inla_data, censor, site_dev, blq_extra))
  )

  inla_samples <- map(fits[names(inla_models)], \(fit) {
    set.seed(seed)
    inla.posterior.sample(n = n_samples, result = fit)
  })

  site_draws <- c(
    list(GAM = predict_gam_site(fits$GAM, grid, n_samples, seed)),
    map(inla_samples, predict_inla_site, inla_data = inla_data)
  )

  deviation_draws <- imap(inla_samples, \(samples, model) {
    inla_deviation_draws(samples, grid, inla_data, inla_models[[model]])
  })

  list(
    grid = grid,
    inla_data = inla_data,
    fits = fits,
    site_draws = site_draws,
    deviation_draws = deviation_draws,
    national = imap(site_draws, \(draws, model) aggregate_national_draws(draws, grid, model, weighted))
  )
}
