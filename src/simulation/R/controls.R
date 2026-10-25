# E3 power-control panel: the engineered aligned-error T2 control and U10 (the stage-2 weights
# omitted). The panel's host is R-a PSU-CF AIPW; each control is calibrated by a single constant
# (kappa = 0.30 for the fixed-shift control); U10 runs in its own oracle cell; the aligned-error
# construction must stay inside the fitted clip bound. G6 limb (a) reads each control against
# the boundary of its compliant region.
#
# C1's control is the fixed shift
# Psi_N + kappa sigma_m,wor, kappa = 0.30; the g-formula-keyed Z is reported without a verdict.
# G6 limb (a)'s
# 3 MCSE test at the planned count is evaluated by runs/run_prep.R.
#
# CONVENTIONS:
# C-AE-1 Aligned-error construction. A deterministic shift of the MARGINAL truth eta0 ("a
#        deterministic nuisance shift with a known nonzero P0-product remainder"), one constant
#        kappa >= 0 on the logit scale of all three nuisances:
#          logit Qhat_a = logit Q0_a + kappa (a = 0, 1),   logit ghat = logit g0 + kappa.
#        The errors are ALIGNED: ghat - g0 > 0 and Qhat_a - Q0_a > 0 everywhere, so under the
#        proofs' sign (Supplement S1, R/remainders.R) both arm terms of P0 f add and
#        T2(kappa) = P0 f > 0, quadratic in kappa near 0. The shift does not depend on the sample,
#        so f_v = f in every fold and T2 is a closed-form constant on the P0 grid.
# C-AE-2 Clip bound. The perturbed nuisances must stay strictly inside the fitted clip
#        (q0 = 0.005 for Q, delta = 0.01 for g), so that clipping never alters the
#        construction. eta0 is a monotone function of multilinear predictors on the covariate box,
#        so its extremes are at the box's vertices; kappa_max is the smallest logit distance from
#        those extremes to 1 - q0 (Q) and 1 - delta (g).
# C-AE-3 Calibration rule (fixed before any production population is
#        generated). kappa_AE is derived ONCE, by census algebra, from the complete prespecified
#        production panel population set -- populations AE_PANEL_POPS = 1..160 of master FS at the
#        panel cell (H = 5, m_h = 80, J = 800, m = 400) -- as the root of
#        Q_0.9(T2(kappa) / sigma_m,wor,b) = 0.30 over those populations; T2 is common to them, so
#        Q_0.9(T2/sigma_b) = T2 * Q_0.9(1/sigma_b) exactly (type-7 quantiles commute with positive
#        scaling). T2 is strictly increasing in kappa on [0, kappa_max] (checked), so the root is
#        unique; if T2(kappa_max) falls short, the control is NOT BUILDABLE. The value and its
#        provenance are saved (ae_frozen_write; runs/frozen/kappa_ae.json, written by
#        runs/run_prep.R --job panel) and reused unchanged by every dry and production job
#        and rerun (ae_frozen_load); it is never recalibrated from a job's subset of populations or
#        from fitted results. No runner solves kappa.
# C-AE-4 The aligned-error arm is PSU-CF AIPW on the host's sample and folds with this
#        deterministic "learner" (it ignores its training set), arm "AE-AIPW". No level-4 stream.
#        U3/U4/U5/U7, U8 and U9's brute-force limb apply (keeps U9 for it).
# C-U10-1 U10 (G8): the oracle Hajek mean with the stage-2 weights omitted, w1 = M_pop,h/m_h
#        per unit. Under the frozen fixed-n0 design N_hat1 = sum w1 = n0 * M_pop is constant, so
#        the estimator is linear: E[psi_U10] = Psi_U10 = (1/M_pop) sum_k phibar_k (the unweighted
#        mean of PSU means) and E[se2_U10] = sigma2_m,U10 exactly (exact enumeration). Its Z is
#        keyed to Psi_N, so it is biased by Psi_U10 - Psi_N whenever size is informative.

PANEL_KAPPA <- 0.30                    # the single shift constant
AE_PANEL_POPS <- 1:160                 # C-AE-3: the complete production panel population set (160 x 5)
AE_PANEL_CELL <- list(regime = "fixed", H = 5L, m_h = 80L, J = 800L, m = 400L, f = "1/10", e2c = "none")
AE_FROZEN_FILE <- file.path("runs", "frozen", "kappa_ae.json")   # relative to the engine directory
AE_RULE <- "C-AE-3 (R23 as amended, 2026-10-11; spec v9.8 Amendment A3 §A3.8)"
ARM_AE_AIPW <- "AE-AIPW"               # C-AE-4
ARM_ORACLE_U10 <- "oracle-U10"         # C-U10-1

# ---- aligned-error construction --------------------------------------------------------

ae_shift <- function(p, kappa) stats::plogis(stats::qlogis(p) + kappa)

# A deterministic nuisance "fit" for psu_cf_aipw(): the truth columns Q1, Q0, g are used when the
# data carry them (sample units), otherwise eta0 is computed at GH64 (P0 grids) -- the same
# function and nodes as the oracle truth, so both paths give identical values.
ae_fit <- function(kappa, co, q0 = 0.005, delta = 0.01, gh_n = 64L) {
  if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa) || kappa < 0)
    stop("kappa must be one finite number >= 0", call. = FALSE)
  gh <- gauss_hermite(gh_n)
  truth <- function(nd, which) {
    col <- switch(which, Q1 = "Q1", Q0 = "Q0", g = "g")
    if (!is.null(nd[[col]])) return(nd[[col]])
    switch(which, Q1 = Q0_marginal(1, nd, nd$KB, co, gh), Q0 = Q0_marginal(0, nd, nd$KB, co, gh),
           g = g0_marginal(nd, nd$KB, co, gh))
  }
  raw <- list(
    Q = function(a, newdata) {
      if (!length(a) %in% c(1L, nrow(newdata)) || !all(a %in% 0:1)) stop("a must be 0/1, scalar or one per row", call. = FALSE)
      a <- rep_len(a, nrow(newdata))
      out <- numeric(nrow(newdata))
      if (any(a == 1)) out[a == 1] <- ae_shift(truth(newdata[a == 1, , drop = FALSE], "Q1"), kappa)
      if (any(a == 0)) out[a == 0] <- ae_shift(truth(newdata[a == 0, , drop = FALSE], "Q0"), kappa)
      out
    },
    g = function(newdata) ae_shift(truth(newdata, "g"), kappa))
  c(clip_fit(raw, q0 = q0, delta = delta),
    list(raw = raw, diagnostics = NULL,
         conventions = list(status = "C-AE-1 (R23)", construction = "C-AE-1 (logit shift of eta0)",
                            kappa = kappa, q0 = q0, delta = delta)))
}

# T2(kappa) in closed form on a P0 grid with its truth: exactly p0_f_closed() for ae_fit(kappa),
# vectorised on the stored truth (used to solve kappa_AE; the runner checks it against the
# generic t2_closed_form() on the AE fits).
ae_t2_map <- function(kappa, truth, w) {
  g <- truth$g; gh <- ae_shift(g, kappa)
  q1 <- ae_shift(truth$Q1, kappa); q0 <- ae_shift(truth$Q0, kappa)
  sum(w * ((q1 - truth$Q1) * (gh - g) / gh - (q0 - truth$Q0) * ((1 - gh) - (1 - g)) / (1 - gh)))
}

# C-AE-2: extremes of eta0 on the covariate box (its 64 vertices x 2 classes) and kappa_max.
ae_clip_bound <- function(co, q0 = 0.005, delta = 0.01, gh_n = 64L) {
  v <- expand.grid(C1 = c(-1, 1), C2 = c(-1, 1), C3 = c(-1, 1), C4 = c(-1, 1), Cpsu = c(-1, 1), KB = 1:2,
                   KEEP.OUT.ATTRS = FALSE)
  tr <- truth_at(v, co, gh_n)
  Q <- c(tr$Q1, tr$Q0)
  list(kappa_max = min(stats::qlogis(1 - delta) - stats::qlogis(max(tr$g)),
                       stats::qlogis(1 - q0) - stats::qlogis(max(Q))),
       min_g = min(tr$g), max_g = max(tr$g), min_Q = min(Q), max_Q = max(Q),
       inside = min(tr$g) > delta && max(tr$g) < 1 - delta && min(Q) > q0 && max(Q) < 1 - q0,
       q0 = q0, delta = delta)
}

# C-AE-3: solve Q_0.9(T2(kappa) / sigma_b) = target on [0, kappa_max]. sigma_wor: the cell's
# populations' sigma_m,wor. Returns the root, the predicted statistic, monotonicity and buildability.
ae_kappa_solve <- function(truth, w, sigma_wor, kappa_max, target = PANEL_KAPPA, n_mono = 64L) {
  if (!length(sigma_wor) || any(!is.finite(sigma_wor)) || any(sigma_wor <= 0)) stop("sigma_wor must be positive", call. = FALSE)
  scale <- unname(stats::quantile(1 / sigma_wor, 0.9, type = 7))        # Q_0.9(T2/sigma) = T2 * scale
  kk <- seq(0, kappa_max, length.out = n_mono)
  t2 <- vapply(kk, ae_t2_map, numeric(1), truth = truth, w = w)
  mono <- all(diff(t2) > 0)
  at_max <- t2[n_mono] * scale
  if (!mono || at_max < target)
    return(list(buildable = FALSE, kappa = NA_real_, T2 = NA_real_, stat = NA_real_, at_max = at_max,
                monotone = mono, scale = scale, kappa_max = kappa_max, target = target,
                reason = if (!mono) "T2(kappa) not strictly increasing on [0, kappa_max]" else
                  sprintf("Q_0.9 at kappa_max is %.4f < target %.2f within the clip bound", at_max, target)))
  r <- stats::uniroot(function(k) ae_t2_map(k, truth, w) * scale - target, c(0, kappa_max), tol = 1e-13)
  T2 <- ae_t2_map(r$root, truth, w)
  list(buildable = TRUE, kappa = r$root, T2 = T2, stat = T2 * scale, at_max = at_max, monotone = mono,
       scale = scale, kappa_max = kappa_max, target = target, reason = "")
}

# ---- C-AE-3: the single frozen calibration and its provenance ------------------------------

# The JSON record (all numbers at full precision); content_sha256 covers every other field.
.ae_canonical <- function(x) jsonlite::toJSON(x[setdiff(names(x), "content_sha256")], auto_unbox = TRUE, digits = I(17), na = "null", null = "null")

ae_frozen_write <- function(path, rec) {
  if (file.exists(path)) stop("the frozen kappa_AE record already exists (", path, "): it is derived once (C-AE-3)", call. = FALSE)
  rec$content_sha256 <- digest::digest(as.character(.ae_canonical(rec)), algo = "sha256", serialize = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(jsonlite::toJSON(rec, auto_unbox = TRUE, digits = I(17), na = "null", null = "null", pretty = TRUE), path, useBytes = TRUE)
  invisible(rec)
}

# Loads and verifies the frozen record: its content hash, the rule, the root, the master, the complete
# population set, the cell, the specification identifier and (if given) the P0 grid hash the runner uses. Any mismatch
# is an error: a job never recalibrates and never runs the AE arm on an unverified value.
ae_frozen_load <- function(path, spec_sha = NULL, grid_hash = NULL) {
  if (!file.exists(path)) stop("no frozen kappa_AE record at ", path, " -- the round-3 preparation job (runs/run_prep.R --task panel) derives it once (C-AE-3)", call. = FALSE)
  rec <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  h <- digest::digest(as.character(.ae_canonical(rec)), algo = "sha256", serialize = FALSE)
  bad <- character(0)
  if (!identical(h, rec$content_sha256)) bad <- c(bad, "content hash")
  if (!identical(rec$rule, AE_RULE)) bad <- c(bad, "rule")
  if (!identical(rec$root, SIM_PROD_ROOT)) bad <- c(bad, "root")
  if (!identical(rec$master_id, "FS")) bad <- c(bad, "master")
  if (!identical(as.integer(rec$population_ids), as.integer(AE_PANEL_POPS))) bad <- c(bad, "population set")
  if (!identical(rng_encode_cell(lapply(rec$cell, function(v) if (is.numeric(v)) as.integer(v) else v)), rng_encode_cell(AE_PANEL_CELL))) bad <- c(bad, "cell")
  if (!is.null(spec_sha) && !identical(rec$spec_sha256, spec_sha)) bad <- c(bad, "spec hash")
  if (!is.null(grid_hash) && !identical(rec$p0_primary_hash, grid_hash)) bad <- c(bad, "P0 grid hash")
  if (length(bad)) stop("the frozen kappa_AE record fails verification (", paste(bad, collapse = ", "), "): ", path, call. = FALSE)
  rec
}

# ---- C1 shift control and the g-formula-keyed Z -----------------------------------

# Z keyed to the fixed artificial shift Psi_N + kappa sigma_m,wor (the graded predicted-to-fail
# form) and to Psi_N^gf (reported, no verdict), on the host arm's psi and se.
panel_shift_z <- function(psi, se, Psi_N, Psi_N_gf, sigma_wor, kappa = PANEL_KAPPA)
  c(Z_shift = (psi - (Psi_N + kappa * sigma_wor)) / se, Z_gf = (psi - Psi_N_gf) / se)

# ---- U10 (G8; its own oracle cell) ---------------------------------------------

u10_estimate <- function(sample) {
  d <- attr(sample, "design")
  if (is.null(d)) stop("u10_estimate needs draw_sample() output", call. = FALSE)
  w1 <- sample$M_pop_h / as.numeric(d$m_h)[sample$h]                  # stage-1 weight only (C-U10-1)
  uc <- hajek_uc(sample$phi0, w1, sample$h, sample$k)
  list(arm = ARM_ORACLE_U10, psi = uc$psi, se2 = uc$se2, se = sqrt(uc$se2), N_hat = uc$N_hat, w1 = w1)
}

# Census algebra for U10 (exact; C-U10-1): Psi_U10, sigma2_m,U10 = E[se2_U10], N_hat1 = n0 M_pop.
u10_census <- function(pop, m_h, n0 = 20) {
  psu <- pop$psu; units <- pop$units
  idx <- rep(seq_len(nrow(psu)), psu$M)
  strata <- sort(unique(psu$h)); sh <- match(psu$h, strata)
  M_pop_h <- as.vector(table(factor(psu$h, levels = strata)))
  if (length(m_h) == 1L) m_h <- rep(m_h, length(strata))
  if (any(psu$M < n0)) stop("u10_census needs n0 <= M_k (fixed-n0 design)", call. = FALSE)
  M_pop <- nrow(psu); N1 <- n0 * M_pop
  phibar <- as.vector(rowsum(units$phi0, idx, reorder = TRUE)) / psu$M
  S2_phi <- as.vector(rowsum((units$phi0 - phibar[idx])^2, idx, reorder = TRUE)) / (psu$M - 1)
  Y <- n0 * phibar
  S2_Y <- as.vector(tapply(Y, sh, stats::var))
  within <- as.vector(rowsum(n0 * (1 - n0 / psu$M) * S2_phi, sh, reorder = TRUE))
  sigma2 <- sum((M_pop_h^2 / m_h) * S2_Y + (M_pop_h / m_h) * within) / N1^2
  list(Psi_U10 = mean(phibar), sigma2_m_u10 = sigma2, N_hat1 = N1, Psi_N = mean(units$phi0),
       bias = mean(phibar) - mean(units$phi0))
}

# ---- predicted values for G6 limb (a), boundary reading ----------------------------

# The panel's controls on one population: predicted (census, first order for the fitted host)
# value of the statistic each control corrupts, with its compliant region (rows):
#   C1 shift          mean Z_shift      region [-0.10, 0.10]    predicted -kappa sigma_wor / sigma_m
#   wrong clustering  mean R_var        region [0.85, 1.15]     (gradeable fitted limb) census limits
#   aligned error     Q_0.9(|T2|/sw)    region [0, 0.10]        kappa_AE's predicted statistic
# U10 (own cell)      mean Z_U10        region [-0.10, 0.10]    (Psi_U10 - Psi_N) / sigma_m,U10
panel_predictions <- function(cen, wc, ae_stat = NA_real_, kappa = PANEL_KAPPA) {
  rv <- wc$R_var
  data.frame(
    control = c("C1 shift: mean Z keyed to Psi_N + 0.30 sigma_m,wor",
                "C2 wrong clustering: no m_h/(m_h-1), mean R_var",
                "C2 wrong clustering: unit-level clusters, mean R_var",
                "C2 wrong clustering: no stratification, mean R_var",
                "C4 aligned error: Q_0.9(|T2| / sigma_m,wor)"),
    statistic = c("mean Z", "mean R_var", "mean R_var", "mean R_var", "Q_0.9(|T2|/sigma_m,wor)"),
    predicted = c(-kappa * sqrt(cen$sigma2_m_wor / cen$sigma2_m), rv[["nomh"]], rv[["unit"]], rv[["nostrat"]], ae_stat),
    lo = c(-0.10, 0.85, 0.85, 0.85, 0), hi = c(0.10, 1.15, 1.15, 1.15, 0.10),
    stringsAsFactors = FALSE)
}
