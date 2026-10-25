# Hajek point estimate, ultimate-cluster (UC) variance and the U8 survey cross-check.
#
# ESTIMATOR (the paper's ultimate-cluster variance estimator):
#   psi   = sum_i w_i phi_i / N_hat, N_hat = sum_i w_i  -- the Hajek ratio with the
#           ESTIMATED N_hat in the denominator, never forced to N (U2)
#   z_i   = w_i (phi_i - psi) / N_hat                   -- linearised Hajek mean
#   t_hj  = sum_{i in PSU (h,j)} z_i
#   se2   = sum_h m_h/(m_h - 1) sum_j (t_hj - tbar_h)^2
# m_h is the number of SAMPLED PSUs in stratum h. NO finite-population correction:
# se2 targets the UC/no-FPC idealisation sigma_m^2 := E_des[se2] at eta0,
# which carries the m_h/(m_h-1) factor and no census-PSU-count factor (U8).
#
# U8, three limbs:
#   single-stage, fpc NOT declared  -- check_u8(): se2 must reproduce
#       survey::svymean(~phi, svydesign(ids=~psu, strata=~stratum, weights=~w)).
#   two-stage, correctly declared  -- se2_two_stage_wor(): ids=~psu+unit,
#       fpc=~M_pop_h+M_k; its design MEAN reproduces sigma2_m_wor. That is a
#       dry-run object over redraws; only the per-sample helper is here.
#   negative (blanket single-stage fpc) -- se2_single_stage_fpc(): its design
#       mean is sigma2_m_wor - sum_h f_h V_h, strictly below sigma2_m_wor (the census deficit,
#       census_fpc_summary()); checked over redraws in the calibration cell.
#
# TECHNICAL CONVENTIONS (labelled):
#   TC-V1 a PSU is the (stratum, psu) PAIR, so PSU ids may repeat across strata
#         (survey's nest = TRUE).
#   TC-V2 a stratum with fewer than 2 sampled PSUs is a hard error: m_h/(m_h-1) is
#         undefined, and every frozen cell has m_h >= 2. No lonely-PSU rule
#         is applied, here or in the survey comparison (survey.lonely.psu = "fail").
#   TC-V3 U8's tolerance: the design states U8 as "numeric" with no number.
#         check_u8() uses relative difference <= 1e-12; both
#         sides are the same algebra, so the observed gap is rounding (~1e-15).

VARIANCE_DEV_LABEL <- "UC variance (spec v9.8)"

.check_var_inputs <- function(phi, w, stratum, psu) {
  n <- length(phi)
  if (!n || length(w) != n || length(stratum) != n || length(psu) != n)
    stop("phi, w, stratum and psu must have one common, positive length", call. = FALSE)
  if (!is.numeric(phi) || any(!is.finite(phi))) stop("phi must be finite", call. = FALSE)
  if (!is.numeric(w) || any(!is.finite(w)) || any(w <= 0)) stop("w must be finite and > 0", call. = FALSE)
  if (anyNA(stratum) || anyNA(psu)) stop("stratum and psu must not contain NA", call. = FALSE)
  invisible(TRUE)
}

# Hajek mean and its UC variance (eq:ucvar). Returns list(psi, N_hat, se2, m_h).
hajek_uc <- function(phi, w, stratum, psu) {
  .check_var_inputs(phi, w, stratum, psu)
  N_hat <- sum(w)
  psi <- sum(w * phi) / N_hat
  z <- w * (phi - psi) / N_hat
  s_code <- match(stratum, unique(stratum))
  p_code <- match(psu, unique(psu))
  pair <- (s_code - 1) * max(p_code) + p_code         # TC-V1: PSU = (stratum, psu)
  t <- rowsum(z, pair, reorder = FALSE)[, 1]          # groups in order of first appearance
  s_of_t <- s_code[!duplicated(pair)]
  m_h <- tabulate(s_of_t)
  if (any(m_h < 2L))
    stop("every stratum needs >= 2 sampled PSUs (UC factor m_h/(m_h-1), TC-V2)", call. = FALSE)
  tbar <- rowsum(t, s_of_t)[, 1] / m_h                # rowsum sorts by code 1..S
  dev <- t - tbar[s_of_t]
  se2 <- sum((m_h / (m_h - 1))[s_of_t] * dev^2)
  list(psi = psi, N_hat = N_hat, se2 = se2,
       m_h = stats::setNames(m_h, as.character(unique(stratum))))
}

# Survey designs are built with lonely PSUs set to fail and ultimate-cluster
# shortcuts off, so a global option cannot silently change what is compared.
.with_survey_options <- function(expr) {
  op <- options(survey.lonely.psu = "fail", survey.ultimate.cluster = FALSE)
  on.exit(options(op), add = TRUE)
  force(expr)
}

# U8 single-stage limb: se2 against survey::svymean with fpc NOT declared.
check_u8 <- function(phi, w, stratum, psu, tol = 1e-12) {
  est <- hajek_uc(phi, w, stratum, psu)
  dat <- data.frame(phi = phi, w = w, stratum = stratum, psu = psu)
  sv <- .with_survey_options({
    des <- survey::svydesign(ids = ~psu, strata = ~stratum, weights = ~w, data = dat, nest = TRUE)
    survey::svymean(~phi, des)
  })
  se2_survey <- as.numeric(attr(sv, "var"))
  psi_survey <- as.numeric(sv)
  rel_diff <- abs(est$se2 - se2_survey) / abs(se2_survey)
  psi_rel_diff <- abs(est$psi - psi_survey) / abs(psi_survey)
  list(se2 = est$se2, se2_survey = se2_survey, rel_diff = rel_diff,
       psi = est$psi, psi_survey = psi_survey, psi_rel_diff = psi_rel_diff, tol = tol,
       pass = is.finite(rel_diff) && is.finite(psi_rel_diff) && rel_diff <= tol && psi_rel_diff <= tol,
       label = VARIANCE_DEV_LABEL)
}

# U8 two-stage limb: per-sample two-stage WOR variance estimate with both
# stages' fpc declared (M_pop_h PSUs in the stratum, M_k units in the PSU). Its
# mean over redraws is the dry-run comparator for sigma2_m_wor; no gate here.
se2_two_stage_wor <- function(phi, w, stratum, psu, unit, M_pop_h, M_k) {
  .check_var_inputs(phi, w, stratum, psu)
  dat <- data.frame(phi = phi, w = w, stratum = stratum, psu = psu, unit = unit,
                    M_pop_h = M_pop_h, M_k = M_k)
  sv <- .with_survey_options({
    des <- survey::svydesign(ids = ~psu + unit, strata = ~stratum, fpc = ~M_pop_h + M_k,
                             weights = ~w, data = dat, nest = TRUE)
    survey::svymean(~phi, des)
  })
  list(psi = as.numeric(sv), se2_wor = as.numeric(attr(sv, "var")),
       label = paste(VARIANCE_DEV_LABEL, "| U8 two-stage limb: dry-run object, not gated"))
}

# U8 negative limb: a BLANKET single-stage design with the stage-1 fpc declared
# (ids = ~psu, fpc = ~M_pop_h). Never used to "repair" U8; computed only to show that its
# design mean falls short of sigma2_m_wor by sum_h f_h V_h. Dry-run object, not gated per sample.
se2_single_stage_fpc <- function(phi, w, stratum, psu, M_pop_h) {
  .check_var_inputs(phi, w, stratum, psu)
  dat <- data.frame(phi = phi, w = w, stratum = stratum, psu = psu, M_pop_h = M_pop_h)
  sv <- .with_survey_options({
    des <- survey::svydesign(ids = ~psu, strata = ~stratum, fpc = ~M_pop_h, weights = ~w, data = dat, nest = TRUE)
    survey::svymean(~phi, des)
  })
  list(psi = as.numeric(sv), se2_fpc1 = as.numeric(attr(sv, "var")),
       label = paste(VARIANCE_DEV_LABEL, "| U8 negative limb: dry-run object, not gated"))
}

# DEV calibration gate: the census-centred LINEARISED score with KNOWN N.
# u_i = (phi_i - Psi_N) / N, and the HT total sum_i w_i u_i is LINEAR in the sample, so the
# variance estimators of the three declared designs have design expectations EXACTLY
#   single-stage, fpc not declared (U8's UC limb)   -> sigma2_m
#   two-stage, fpc = ~M_pop_h + M_k (two-stage limb) -> sigma2_m_wor
#   blanket single-stage fpc = ~M_pop_h (negative)   -> sigma2_m_wor - sum_h f_h V_h
# and E[sum_i w_i u_i] = 0, all by census algebra (proved by exact enumeration in
# tests/testthat/test-variance-linearised.R). The Hajek estimators (psi_hat and N_hat in the
# residuals: check_u8, se2_two_stage_wor, se2_single_stage_fpc) differ from these in finite
# samples by O(1/m), so only the linearised versions are gated over many redraws; the Hajek
# ones stay descriptive. No estimator is changed.
lin_var_designs <- function(phi, Psi_N, N, w, stratum, psu, unit, M_pop_h, M_k) {
  .check_var_inputs(phi, w, stratum, psu)
  dat <- data.frame(u = (phi - Psi_N) / N, w = w, stratum = stratum, psu = psu, unit = unit,
                    M_pop_h = M_pop_h, M_k = M_k)
  .with_survey_options({
    tot <- function(des) { t <- survey::svytotal(~u, des); c(est = as.numeric(t), var = as.numeric(attr(t, "var"))) }
    uc <- tot(survey::svydesign(ids = ~psu, strata = ~stratum, weights = ~w, data = dat, nest = TRUE))
    wor <- tot(survey::svydesign(ids = ~psu + unit, strata = ~stratum, fpc = ~M_pop_h + M_k, weights = ~w,
                                 data = dat, nest = TRUE))
    fpc1 <- tot(survey::svydesign(ids = ~psu, strata = ~stratum, fpc = ~M_pop_h, weights = ~w, data = dat, nest = TRUE))
    list(total = uc[["est"]], se2_uc = uc[["var"]], se2_wor = wor[["var"]], se2_fpc1 = fpc1[["var"]],
         label = paste(VARIANCE_DEV_LABEL, "| linearised score, known N (calibration gate)"))
  })
}

# ---- E2 WR and E3 control helpers -------------------------------------------

# E2 WR code-check: the census-centred LINEARISED score with KNOWN N under the
# with-replacement design. u_i = (phi_i - Psi_N)/N; every selection (h, j) is its own ultimate
# cluster (draw_sample_wr). The UC/WR variance of the total sum_i w_i u_i has design expectation
# EXACTLY sigma2_m,WR = sum_h [(1 - 1/M_pop,h) B_h + V_h] (census_wr_target), and E[total] = 0
# (exact enumeration, tests/testthat/test-sampling-wr.R). GATED over redraws. Two paths: survey's
# svytotal (ids = ~selection, strata, no fpc: a WR design) and the direct cluster-total formula.
lin_var_wr <- function(phi, Psi_N, N, w, stratum, sel) {
  .check_var_inputs(phi, w, stratum, sel)
  u <- (phi - Psi_N) / N
  dat <- data.frame(u = u, w = w, stratum = stratum, sel = sel)
  sv <- .with_survey_options({
    t <- survey::svytotal(~u, survey::svydesign(ids = ~sel, strata = ~stratum, weights = ~w, data = dat, nest = TRUE))
    c(est = as.numeric(t), var = as.numeric(attr(t, "var")))
  })
  direct <- .uc_var_of_totals(w * u, stratum, sel)
  list(total = sv[["est"]], se2 = sv[["var"]], se2_direct = direct, total_direct = sum(w * u),
       label = paste(VARIANCE_DEV_LABEL, "| WR code-check: linearised score, known N (gate)"))
}

# UC variance of a sum of per-unit scores z over clusters within strata:
# sum_h m_h/(m_h - 1) sum_c (t_c - tbar_h)^2, t_c the cluster total of z (TC-V1 pairing).
.uc_var_of_totals <- function(z, stratum, cluster, factor = TRUE) {
  s_code <- match(stratum, unique(stratum))
  c_code <- match(cluster, unique(cluster))
  pair <- (s_code - 1) * max(c_code) + c_code
  t <- rowsum(z, pair, reorder = FALSE)[, 1]
  s_of_t <- s_code[!duplicated(pair)]
  m_h <- tabulate(s_of_t)
  if (any(m_h < 2L)) stop("every stratum needs >= 2 clusters", call. = FALSE)
  tbar <- rowsum(t, s_of_t)[, 1] / m_h
  k <- if (factor) (m_h / (m_h - 1))[s_of_t] else 1
  sum(k * (t - tbar[s_of_t])^2)
}

# E3 power-control panel, wrong-clustering se2 control (host: R-a PSU-CF AIPW):
# the three "computable wrong limits" of the UC estimator on a score phi (the host's fitted IF),
# each on the Hajek linearised z_i = w_i (phi_i - psi_hat) / N_hat:
#   uc       the correct UC estimator (= hajek_uc()$se2)
#   nomh     no m_h/(m_h - 1) factor
#   unit     units as clusters within strata (unit-level clustering)
#   nostrat  PSUs as clusters, strata ignored (no stratification)
wrong_clustering_se2 <- function(phi, w, stratum, psu) {
  .check_var_inputs(phi, w, stratum, psu)
  N_hat <- sum(w); psi <- sum(w * phi) / N_hat
  wrong_clustering_from_z(w * (phi - psi) / N_hat, stratum, psu)
}

# The four estimators from per-unit scores z (Hajek: z = w (phi - psi_hat)/N_hat; linearised with
# known N: z = w (phi - Psi_N)/N, whose expectations wrong_clustering_census() gives exactly).
wrong_clustering_from_z <- function(z, stratum, psu) {
  c(uc = .uc_var_of_totals(z, stratum, psu),
    nomh = .uc_var_of_totals(z, stratum, psu, factor = FALSE),
    unit = .uc_var_of_totals(z, stratum, seq_along(z)),
    nostrat = .uc_var_of_totals(z, rep(1L, length(z)), paste(stratum, psu, sep = ":")))
}

# Independent path for the wrong-clustering variants: survey::svymean under the wrong designs
# (unit: ids = ~1 with strata; nostrat: ids = ~psu without strata). survey has no switch that
# drops m_h/(m_h - 1), so nomh is checked as uc * (m_h - 1)/m_h, valid when every stratum has the
# same m_h (true in the panel); otherwise NA.
wrong_clustering_survey <- function(phi, w, stratum, psu) {
  dat <- data.frame(phi = phi, w = w, stratum = stratum, psu = paste(stratum, psu, sep = ":"))
  .with_survey_options({
    v <- function(des) as.numeric(attr(survey::svymean(~phi, des), "var"))
    uc <- v(survey::svydesign(ids = ~psu, strata = ~stratum, weights = ~w, data = dat, nest = TRUE))
    m_h <- tapply(dat$psu, dat$stratum, function(x) length(unique(x)))
    c(uc = uc,
      nomh = if (length(unique(m_h)) == 1L) uc * (m_h[[1]] - 1) / m_h[[1]] else NA_real_,
      unit = v(survey::svydesign(ids = ~1, strata = ~stratum, weights = ~w, data = dat)),
      nostrat = v(survey::svydesign(ids = ~psu, weights = ~w, data = dat)))
  })
}

# Census predictions of the four estimators above, EXACT for the census-centred linearised score
# with known N at eta0 (u = (phi0 - Psi_N)/N) under two-stage SRS-WOR with fixed n0 (exact
# enumeration in tests/testthat/test-controls.R). The fitted host differs only to first order
# (the "computable wrong limits" of are oracle-score census quantities). R_var limit =
# value / sigma2_m. Per stratum h (M_pop,h PSUs, m_h sampled, n0 units each):
#   E sum_{i in s_h} z_i^2 = (M_pop,h/m_h) sum_k (M_k/n0) sum_{i in k} u_i^2
#   E (sum_{i in s_h} z_i)^2 = (1 - f_h) B_h + V_h + T_h^2,   T_h = sum_{i in h} u_i
#   unit:    sum_h n_h/(n_h - 1) [E sum z^2 - E (sum z)^2 / n_h],  n_h = m_h n0
#   nomh:    sum_h ((m_h - 1)/m_h)(B_h + V_h)
#   nostrat: m/(m - 1) [sum_h (M_pop,h/m_h) sum_k (Y_k^2 + M_k (M_k - n0) S2_u,k / n0) - (sigma2_m_wor + T^2)/m]
wrong_clustering_census <- function(pop, m_h, n0 = 20) {
  cen <- census_targets(pop, m_h = m_h, n0 = n0)
  psu <- pop$psu; units <- pop$units
  idx <- rep(seq_len(nrow(psu)), psu$M)
  u <- (units$phi0 - cen$Psi_N) / cen$N
  ps <- cen$per_stratum
  sh <- match(psu$h, ps$h)
  Y_k <- as.vector(rowsum(u, idx, reorder = TRUE))
  U2_k <- as.vector(rowsum(u^2, idx, reorder = TRUE))
  S2_k <- as.vector(rowsum((u - (Y_k / psu$M)[idx])^2, idx, reorder = TRUE)) / (psu$M - 1)
  mh <- ps$m_h[sh]; Mh <- ps$M_pop_h[sh]
  Ez2 <- as.vector(rowsum((Mh / mh) * (psu$M / n0) * U2_k, sh, reorder = TRUE))
  T_h <- as.vector(rowsum(Y_k, sh, reorder = TRUE))
  var_h <- (1 - ps$f_h) * ps$B_h + ps$V_h
  n_h <- ps$m_h * n0
  unit <- sum(n_h / (n_h - 1) * (Ez2 - (var_h + T_h^2) / n_h))
  nomh <- sum((ps$m_h - 1) / ps$m_h * (ps$B_h + ps$V_h))
  m <- sum(ps$m_h); T_all <- sum(u)
  Et2 <- sum((Mh / mh) * (Y_k^2 + psu$M * (psu$M - n0) * S2_k / n0))
  nostrat <- m / (m - 1) * (Et2 - (cen$sigma2_m_wor + T_all^2) / m)
  list(sigma2_m = cen$sigma2_m, sigma2_m_wor = cen$sigma2_m_wor,
       expect = c(uc = cen$sigma2_m, nomh = nomh, unit = unit, nostrat = nostrat),
       R_var = c(uc = 1, nomh = nomh, unit = unit, nostrat = nostrat) / c(1, cen$sigma2_m, cen$sigma2_m, cen$sigma2_m),
       cen = cen)
}
