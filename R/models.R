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
# constraint_weights: one weight per plant (in wwtp_id order) for the per-bin
# sum to zero of the deviations. Equal weights make the national RW2 the
# plants' plain average on the log scale; population weights make it the
# population-weighted one (see national_term_draws()). For AR1 the weights
# are fixed over time, so plants not yet surveilled still count in early bins
# (with deviations ~0); RW1 constrains active plants only.
fit_inla <- function(inla_data, censor, site_dev = c("ar1", "rw1"),
                     blq_extra = FALSE, constraint_weights = NULL) {
  site_dev <- match.arg(site_dev)
  inla_data <- inla_data |> mutate(wwtp_dev = wwtp_id)
  n_wwtp <- max(inla_data$wwtp_id)
  if (is.null(constraint_weights)) constraint_weights <- rep(1, n_wwtp)

  if (site_dev == "ar1") {
    # `extraconstr` on a grouped term is applied within each group, i.e. in
    # every time bin: sum_i w_i * d_it = 0.
    ar1_constr <- list(A = matrix(constraint_weights, nrow = 1), e = 0)
    latent <- ~ 1 +
      f(wwtp_id, model = "iid", constr = TRUE) +
      f(time_ch, model = "rw2", scale.model = TRUE, constr = TRUE) +
      f(wwtp_dev, model = "iid", n = n_wwtp, constr = FALSE, extraconstr = ar1_constr,
        group = time_wwtp, control.group = list(model = "ar1"))
  } else {
    n_wwtp <- max(inla_data$wwtp_id)
    n_time <- max(inla_data$time_wwtp)
    inla_data <- inla_data |> mutate(wwtp_dev = (time_wwtp - 1L) * n_wwtp + wwtp_id)
    site_precision <- rw1_site_precision(n_wwtp, n_time)
    site_constr <- rw1_site_constraints(inla_data, constraint_weights)
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
#   - per bin, the active plants' weighted deviations sum to zero (as in the
#     AR1 version, minus plants not yet surveilled);
#   - per plant, its deviations average to zero over its active bins, as in
#     Riou et al. The RW1 does not penalise a plant's overall level, which
#     would otherwise trade off against the iid baseline.
# With equal weights both sets add up every active cell once, so one plant
# row is redundant and is dropped.
rw1_site_constraints <- function(inla_data, weights = NULL) {
  n_wwtp <- max(inla_data$wwtp_id)
  n_time <- max(inla_data$time_wwtp)
  active <- inla_data |> distinct(wwtp_id, time_wwtp)
  if (is.null(weights)) weights <- rep(1, n_wwtp)
  index <- (active$time_wwtp - 1L) * n_wwtp + active$wwtp_id

  by_bin <- split(seq_along(index), active$time_wwtp)
  by_wwtp <- split(index, active$wwtp_id)

  A <- matrix(0, nrow = length(by_bin) + length(by_wwtp), ncol = n_wwtp * n_time)
  for (k in seq_along(by_bin)) {
    cells <- by_bin[[k]]
    A[k, index[cells]] <- weights[active$wwtp_id[cells]]
  }
  for (k in seq_along(by_wwtp)) A[length(by_bin) + k, by_wwtp[[k]]] <- 1
  if (qr(A)$rank < nrow(A)) A <- A[-nrow(A), , drop = FALSE]
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

# The model's own national series: exp(alpha + sum_i w_i b_i + omega_t), per
# draw and prediction day. With the deviations' weighted sum fixed at zero,
# this is exactly the w-weighted average of the plants' curves on the log
# scale (a weighted geometric mean), so no aggregation of plant draws is
# needed. Same output shape as aggregate_national_draws().
national_term_draws <- function(samples, m_inla, grid, weights, model_name) {
  index <- grid |> distinct(day, date, time) |> arrange(day)
  latent_names <- rownames(samples[[1]]$latent)
  row_of <- function(term, ids) match(paste0(term, ":", ids), latent_names)

  i_alpha <- match("(Intercept):1", latent_names)
  i_omega <- row_of("time_ch", match(index$time, m_inla$summary.random$time_ch$ID))
  i_level <- row_of("wwtp_id", seq_along(weights))
  if (anyNA(c(i_alpha, i_omega, i_level))) stop("Could not find the national terms in the posterior samples.")
  w <- weights / sum(weights)

  draws <- t(vapply(samples, function(s) {
    x <- s$latent[, 1]
    exp(x[i_alpha] + sum(w * x[i_level]) + x[i_omega])
  }, numeric(nrow(index))))

  index <- dplyr::select(index, day, date)
  list(
    summary = bind_cols(index, summarise_draw_columns(draws)) |> mutate(model = model_name),
    draws = draws,
    index = index,
    model = model_name
  )
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
# constraint_weights: "equal" or "population" weights for the deviations'
# per-bin sum to zero (see fit_inla()). Each INLA model also returns its
# national term as an extra national series, labelled "<model> national term".
fit_models <- function(dat, time_step, n_samples, seed, weighted, inla_models,
                       blq_handling = c("censored", "observed", "observed_extra"),
                       zero_handling = c("censor", "drop"),
                       constraint_weights = c("equal", "population")) {
  blq_handling <- match.arg(blq_handling)
  zero_handling <- match.arg(zero_handling)
  constraint_weights <- match.arg(constraint_weights)
  grid <- build_site_time_grid(dat, time_step)
  plant_weights <- if (constraint_weights == "population") {
    site_info <- get_site_info(dat)
    site_info$population[match(levels(dat$wwtp), site_info$wwtp)]
  } else {
    rep(1, nlevels(dat$wwtp))
  }
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
    map(inla_models, \(site_dev) fit_inla(inla_data, censor, site_dev, blq_extra, plant_weights))
  )

  inla_samples <- map(fits[names(inla_models)], \(fit) {
    set.seed(seed)
    inla.posterior.sample(n = n_samples, result = fit)
  })

  site_draws <- c(
    list(GAM = predict_gam_site(fits$GAM, grid, n_samples, seed)),
    map(inla_samples, predict_inla_site, inla_data = inla_data)
  )

  national_terms <- imap(inla_samples, \(samples, model) {
    national_term_draws(samples, fits[[model]], grid, plant_weights, paste(model, "national term"))
  })
  names(national_terms) <- paste(names(national_terms), "national term")

  list(
    grid = grid,
    inla_data = inla_data,
    fits = fits,
    plant_weights = plant_weights,
    site_draws = site_draws,
    national = c(
      imap(site_draws, \(draws, model) aggregate_national_draws(draws, grid, model, weighted)),
      national_terms
    )
  )
}
