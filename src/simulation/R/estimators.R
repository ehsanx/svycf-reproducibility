# Estimator arms: the oracle AIPW
# and the PSU-cross-fitted (PSU-CF) AIPW with the R-a learner, plus the G8 checks
# that bind them (U3, U4/U5 on the folds actually used, U7) and the T2/T1
# remainder decomposition.
#
# The fold stream is the level-3 key of the random-number contract (R/rng_streams.R),
# under a SYNTHETIC root (tests) or the production root. The
# runners label every output with their run mode (R/run_mode.R).
#
# Requires, sourced beforehand: R/quadrature.R, R/nuisance_truth.R, R/p0_grid.R,
# R/remainders.R (aipw_phi, check_omega), R/rng_streams.R, R/folds.R, R/variance.R
# (hajek_uc) and, for the R-a arm, R/learner_ra.R. Nothing here sources anything.
#
# DEFINITIONS (Supplement S1):
#   oracle   psi_oracle = sum_i w_i phi0_i / N_hat, the Hajek mean of the uncentred
#            EIF at the known marginal eta0; N_hat = sum_i w_i is an
#            ESTIMATOR, never forced to N.
#   PSU-CF   strict whole-PSU folds, V fixed; fold v's nuisance is fit
#            on the other folds' PSUs only and scores fold v's units;
#            phi_hat_i = phi(O_i; eta_hat_{-v(i)}) with the post-clip/post-truncation
#            nuisances.
#            psi_CF = sum_v omega_v (P_n^w)_v phi_hat, omega_v = N_hat_v / N_hat,
#            (P_n^w)_v f = sum_{i in v} w_i f_i / N_hat_v (U3, Supplement S1).
#   SE       UC/no-FPC hajek_uc() on the fitted IF phi_hat (U8).
#   U3       sum_v omega_v = 1 and sum_v N_hat_v = N_hat to <= 64 eps (N_hat), on a
#            separate code path.
#   U7       sum_v omega_v (P_n^w)_v f_v = psi_CF - psi_oracle to <= 64 eps max(|psi_oracle|, 1),
#            f_v = phi(.; eta_hat_{-v}) - phi(.; eta0), on a code path that does not
#            reuse psi_CF. AIPW arms only.
#   T2, T1   T2 = sum_v omega_v P0 f_v by t2_closed_form() on a frozen H-specific P0
#            grid; T1 = (psi_CF - psi_oracle) - T2 by subtraction.
#
# TECHNICAL CONVENTIONS:
#   TC-E1 fold draw (contract level 3): key (root, cell, population id,
#         rep id, fold-structure id); one uniform per sampled PSU from entity
#         stratum(h), in ascending frame index k; the strict builder ranks them.
#         fold_structure_id = "psu-s1" is the primary fold structure.
#   TC-E2 U7 is recomputed from scratch: every unit is re-scored with the fit of its
#         own fold label, fold masses are re-summed with which(), and psi_CF is not
#         read. So U7 tests the partition, held-out scoring, omega_v and the fold
#         averages together, and can fail (tested).
#   TC-E3 U3's separate path re-sums N_hat_v with rowsum() over the fold labels and
#         N_hat with sum(); the arm itself sums over fold_rows(). Agreement between
#         the two paths is checked to 64 eps N_hat as well.
#   TC-E4 R-a is a deterministic function of its training set (C-RA-7), so it draws
#         nothing from the level-4 learner stream; no learner key is used.
#   TC-E5 arm names: "oracle", "Ra-AIPW".

EST_DEV_LABEL <- "engine estimator arms (spec v9.8); the runner's run-mode label says what a run is"

ARM_ORACLE <- "oracle"          # TC-E5
ARM_RA_AIPW <- "Ra-AIPW"
FOLD_STRUCTURE_PRIMARY <- "psu-s1"   # TC-E1, draft "<psu|unit>-s<seed index>"

.assert_est_root <- function(root) assert_engine_root(root, "fold stream")

# ---- folds ------------------------------------------------------------------

# The level-3 fold draw for make_strict_psu_folds() (TC-E1). Returns a pure
# function(fold_key, psus): its only randomness is the keyed, scoped stream, so the
# builder's purity check passes and the caller's RNG is untouched.
dev_fold_draw <- function(root) {
  .assert_est_root(root)
  draw <- function(fold_key, psus) {
    h <- psus$stratum
    # The builder hands over its canonical (stratum, psu) table; the per-stratum
    # draw order below relies on it, so assert rather than assume.
    if (any(order(h, psus$psu, method = "radix") != seq_along(h)))
      stop("fold draw: psus must be in canonical (stratum, psu) order", call. = FALSE)
    ctx <- c(fold_key, list(root = root))
    hs <- unique(h)
    st <- rng_states(rng_entity_keys("fold", ctx, "stratum", list(h = hs)))
    score <- numeric(length(h))
    rng_scoped(for (j in seq_along(hs)) {
      i <- which(h == hs[j])            # ascending frame index k within stratum h
      rng_set_state(st[j, ])
      score[i] <- stats::runif(length(i))
    })
    score
  }
  attr(draw, "convention") <- paste0(
    RNG_DRAFT_VERSION, " level 3 (ADOPTED 2026-10-11, R2; docs/RNG-CONVENTIONS-DRAFT.md §7): one uniform per ",
    "sampled PSU from entity stratum(h), ascending k; root ", root)
  draw
}

# Strict V-fold PSU partition of one realised sample (U4 level 3).
# key is the sample key list(cell, population_id, rep_id, ...); only the level-3
# fields reach the fold stream (V never enters it).
dev_psu_folds <- function(sample, key, V, root, fold_structure_id = FOLD_STRUCTURE_PRIMARY) {
  fold_key <- list(cell = key$cell, population_id = key$population_id, rep_id = key$rep_id,
                   fold_structure_id = fold_structure_id)
  d <- attr(sample, "design")
  m_h_design <- if (is.null(d)) NULL else d$m_h
  make_strict_psu_folds(sample$h, sample$k, V, fold_key, dev_fold_draw(root),
                        m_h_design = m_h_design)
}

# ---- arms -------------------------------------------------------------------

# Oracle AIPW: eta0 is known, so the arm is the Hajek mean of phi0
# with its UC variance. The same number is psi_oracle inside every fitted
# replicate (level 2). psu: the ultimate-cluster id within stratum -- the PSU k, or the
# with-replacement selection j in the E2 WR code-check (draw_sample_wr).
oracle_aipw <- function(sample, psu = sample$k) {
  uc <- hajek_uc(sample$phi0, sample$w, sample$h, psu)
  list(arm = ARM_ORACLE, psi = uc$psi, se2 = uc$se2, se = sqrt(uc$se2), N_hat = uc$N_hat,
       phi = sample$phi0, label = EST_DEV_LABEL)
}

# PSU-CF AIPW. learner(train) returns list(Q = function(a, newdata),
# g = function(newdata)) already clipped/truncated as the arm uses it; for
# R-a pass function(train) fit_ra(train, a_cpsu_term = ...). fold is the output of
# make_strict_psu_folds() for this sample.
psu_cf_aipw <- function(sample, fold, learner, arm = ARM_RA_AIPW) {
  V <- .fold_V(fold)
  n <- nrow(sample)
  if (length(fold) != n) stop("fold must carry one label per sampled unit", call. = FALSE)
  # U4 on the vector this arm actually uses, independent recount.
  assert_strict_psu_folds(fold, sample$h, sample$k, V)
  rows <- fold_rows(fold)
  fits <- vector("list", V)
  fit_rows <- score_rows <- vector("list", V)
  Q1 <- Q0 <- g <- rep(NA_real_, n)
  clip_hits <- 0L
  for (v in seq_len(V)) {
    score <- rows[[v]]
    train <- sort(unlist(rows[-v], use.names = FALSE))
    # A learner that takes `fold` gets the fold index (R-b keys its level-4 stream by it);
    # R-a's learner(train) is called exactly as before.
    fits[[v]] <- if ("fold" %in% names(formals(learner))) learner(sample[train, , drop = FALSE], fold = v) else
      learner(sample[train, , drop = FALSE])
    nd <- sample[score, , drop = FALSE]
    Q1[score] <- fits[[v]]$Q(1, nd)
    Q0[score] <- fits[[v]]$Q(0, nd)
    g[score] <- fits[[v]]$g(nd)
    fit_rows[[v]] <- train
    score_rows[[v]] <- score
    # / G7: fitted clipping-hit count, from the pre-clip marginals if exposed.
    if (!is.null(fits[[v]]$raw)) {
      q0c <- fits[[v]]$conventions$q0; dc <- fits[[v]]$conventions$delta
      rq1 <- fits[[v]]$raw$Q(1, nd); rq0 <- fits[[v]]$raw$Q(0, nd); rg <- fits[[v]]$raw$g(nd)
      clip_hits <- clip_hits + sum(rq1 < q0c | rq1 > 1 - q0c | rq0 < q0c | rq0 > 1 - q0c |
                                     rg < dc | rg > 1 - dc)
    }
  }
  # U5: fitting and scored PSU sets disjoint, fold by fold. Hard error.
  check_u5_disjoint(sample$h, sample$k, fit_rows, score_rows, cf_scheme = "psu")
  phi <- aipw_phi(sample$Y, sample$A, Q1, Q0, g)

  # Fold masses and held-out Hajek fold averages (U3; Supplement S1).
  N_hat <- sum(sample$w)
  N_hat_v <- vapply(rows, function(r) sum(sample$w[r]), numeric(1))
  omega <- N_hat_v / N_hat
  fold_means <- vapply(seq_len(V), function(v) {
    r <- rows[[v]]
    sum(sample$w[r] * phi[r]) / N_hat_v[v]
  }, numeric(1))
  psi <- sum(omega * fold_means)
  # UC SE on the fitted IF (U8); its own psi is the pooled Hajek mean, which
  # equals psi to rounding (kept for the record, not used).
  uc <- hajek_uc(phi, sample$w, sample$h, sample$k)
  diag <- do.call(rbind, lapply(seq_len(V), function(v) {
    d <- fits[[v]]$diagnostics
    if (is.null(d)) NULL else cbind(fold = v, d)
  }))
  list(arm = arm, psi = psi, se2 = uc$se2, se = sqrt(uc$se2), psi_pooled = uc$psi,
       N_hat = N_hat, N_hat_v = N_hat_v, omega = omega, fold_means = fold_means,
       phi = phi, Q1 = Q1, Q0 = Q0, g = g, fits = fits, fit_rows = fit_rows,
       score_rows = score_rows, fold = fold, V = V, clip_hits = clip_hits,
       diagnostics = diag, label = EST_DEV_LABEL)
}

# ---- checks -----------------------------------------------------------------

.eps <- .Machine$double.eps

.check_rows <- function(check, value, reference, tolerance, pass) {
  data.frame(check = check, value = value, reference = reference,
             abs_diff = abs(value - reference), tolerance = tolerance,
             pass = !is.na(pass) & pass, stringsAsFactors = FALSE)
}

# U3 (bit-exact class, 64 eps): the arm's omega_v and N_hat_v against a
# separate path (TC-E3).
check_u3 <- function(w, fold, omega, N_hat_v) {
  V <- .fold_V(fold)
  lab <- factor(as.vector(fold), levels = seq_len(V))
  Nv_sep <- as.vector(rowsum(w, lab, reorder = TRUE))
  N_sep <- sum(w)
  tolN <- 64 * .eps * N_sep
  ok_fin <- all(is.finite(c(omega, N_hat_v, Nv_sep, N_sep)))
  .check_rows(
    check = c("U3 sum_v omega_v = 1 (64 eps)",
              "U3 sum_v N_hat_v = N_hat (64 eps N_hat)",
              "U3 N_hat_v reproduced on a separate path (64 eps N_hat)",
              "U3 omega_v = N_hat_v/N_hat on a separate path (64 eps)",
              "U3 omega_v >= 0 and N_hat_v > 0"),
    value = c(sum(omega), sum(N_hat_v), max(abs(N_hat_v - Nv_sep)), max(abs(omega - Nv_sep / N_sep)),
              as.numeric(all(omega >= 0) && all(N_hat_v > 0))),
    reference = c(1, N_sep, 0, 0, 1),
    tolerance = c(64 * .eps, tolN, tolN, 64 * .eps, 0),
    pass = ok_fin & c(abs(sum(omega) - 1) <= 64 * .eps,
                      abs(sum(N_hat_v) - N_sep) <= tolN,
                      max(abs(N_hat_v - Nv_sep)) <= tolN,
                      max(abs(omega - Nv_sep / N_sep)) <= 64 * .eps,
                      all(omega >= 0) && all(N_hat_v > 0)))
}

# U4 and U5 re-run on the arm's recorded folds and row sets; a
# violation is recorded as a failed check instead of stopping the run.
check_u4_u5 <- function(sample, cf) {
  msg <- character(2)
  u4 <- tryCatch({ assert_strict_psu_folds(cf$fold, sample$h, sample$k, cf$V); TRUE },
                 error = function(e) { msg[1] <<- paste(" --", conditionMessage(e)); FALSE })
  u5 <- tryCatch(identical(check_u5_disjoint(sample$h, sample$k, cf$fit_rows, cf$score_rows,
                                             cf_scheme = "psu"), "pass"),
                 error = function(e) { msg[2] <<- paste(" --", conditionMessage(e)); FALSE })
  # every unit scored exactly once (the scored sets partition the sample)
  u5b <- identical(sort(unlist(cf$score_rows, use.names = FALSE)), seq_len(nrow(sample)))
  .check_rows(
    check = c(paste0("U4 strict folds: V | m_h, m_vh = m_h/V per stratum, whole PSUs (independent recount)", msg[1]),
              paste0("U5 fitting and scored PSU sets disjoint in every fold", msg[2]),
              "U5 scored sets partition the sample (every unit scored once)"),
    value = as.numeric(c(u4, u5, u5b)), reference = 1, tolerance = 0,
    pass = c(u4, u5, u5b))
}

# U7: the P0-free fold identity, recomputed from scratch (TC-E2). fits[[v]]
# is the nuisance trained without fold v; psi_cf is used ONLY as the comparator.
check_u7 <- function(sample, fold, fits, psi_cf, psi_oracle) {
  V <- .fold_V(fold)
  lab <- as.vector(fold)
  f <- rep(NA_real_, nrow(sample))
  for (v in seq_len(V)) {
    i <- which(lab == v)
    nd <- sample[i, , drop = FALSE]
    f[i] <- aipw_phi(sample$Y[i], sample$A[i], fits[[v]]$Q(1, nd), fits[[v]]$Q(0, nd), fits[[v]]$g(nd)) -
      sample$phi0[i]
  }
  Nv <- vapply(seq_len(V), function(v) sum(sample$w[lab == v]), numeric(1))
  Pv <- vapply(seq_len(V), function(v) sum(sample$w[lab == v] * f[lab == v]) / Nv[v], numeric(1))
  lhs <- sum((Nv / sum(sample$w)) * Pv)
  rhs <- psi_cf - psi_oracle
  tol <- 64 * .eps * max(abs(psi_oracle), 1)          # floored at 1, so a near-zero psi cannot fail correct code
  out <- .check_rows("U7 P0-free fold identity: sum_v omega_v (P_n^w)_v f_v = psi_CF - psi_oracle (64 eps max(|psi_oracle|, 1), R17)",
                     lhs, rhs, tol, is.finite(lhs) && is.finite(rhs) && abs(lhs - rhs) <= tol)
  attr(out, "fold_averages") <- Pv
  out
}

# Every listed output finite (a NaN/Inf can never read as a pass).
check_finite <- function(what, values) {
  ok <- length(values) > 0L && all(is.finite(unlist(values, use.names = FALSE)))
  .check_rows(paste("finite outputs:", what), as.numeric(ok), 1, 0, ok)
}

# Exit status of a strict run: 0 only if every REQUIRED check passed. A missing
# or NA pass counts as a failure, as in u9_gate_status() (R/remainders.R).
required_status <- function(checks) {
  req <- checks[checks$required, , drop = FALSE]
  if (!nrow(req) || anyNA(req$pass) || !all(req$pass)) 1L else 0L
}

# ---- E3 arms: individual-CF and no-CF --------
#
# individual-CF (the adversarial cluster-leakage stress): V folds at the UNIT level,
# so a PSU's units fall in several folds and its other units train the nuisance that scores it.
# U4's unit analogue applies (every unit exactly one label); U5 is false BY DESIGN
# and is asserted the other way round (the leakage is real); U3 and U7 apply unchanged.
# no-CF: one fit on the whole sample scores that same sample (fold 0, RNG draft vocabulary);
# U3/U4/U5/U7 are vacuous or false by construction and are skipped, never asserted.
#
# TC-E6: unit folds -- within each stratum every sampled unit
# gets one uniform from the level-3 fold stream (key: cell, population id, rep id,
# fold-structure id "unit-s1"; entity stratum(h)), units taken in canonical (k, i) order; the
# units are ordered by their uniforms (ties by canonical position) and the ordered list is
# cut into V consecutive groups whose sizes differ by at most one. Scoped stream.

FOLD_STRUCTURE_UNIT <- "unit-s1"     # RNG draft: "<psu|unit>-s<seed index>"

dev_unit_folds <- function(sample, key, V, root, fold_structure_id = FOLD_STRUCTURE_UNIT) {
  .assert_est_root(root)
  V <- as.integer(V)
  if (length(V) != 1L || is.na(V) || V < 2L) stop("V must be an integer >= 2", call. = FALSE)
  ctx <- list(cell = key$cell, population_id = key$population_id, rep_id = key$rep_id,
              fold_structure_id = fold_structure_id, root = root)
  hs <- sort(unique(sample$h))
  st <- rng_states(rng_entity_keys("fold", ctx, "stratum", list(h = hs)))
  ord_all <- order(sample$h, sample$k, sample$i, method = "radix")
  fold <- integer(nrow(sample))
  rng_scoped(for (j in seq_along(hs)) {
    idx <- ord_all[sample$h[ord_all] == hs[j]]          # canonical (k, i) order within stratum
    rng_set_state(st[j, ])
    u <- stats::runif(length(idx))
    o <- idx[order(u, seq_along(idx))]
    sizes <- diff(round(seq(0, length(o), length.out = V + 1L)))
    fold[o] <- rep(seq_len(V), times = sizes)
  })
  attr(fold, "fold_scheme") <- "unit"
  attr(fold, "V") <- V
  attr(fold, "fold_key") <- ctx
  attr(fold, "draw_convention") <- paste0(RNG_DRAFT_VERSION, " level 3, fold structure ", fold_structure_id,
                                          ": one uniform per unit per stratum, V consecutive groups (TC-E6, R15)")
  fold
}

# Individual-CF AIPW: psu_cf_aipw's scoring and fold algebra on unit-level folds, without the
# whole-PSU assertions. learner(train, fold = v) as for PSU-CF.
indiv_cf_aipw <- function(sample, fold, learner, arm) {
  if (!identical(attr(fold, "fold_scheme"), "unit")) stop("indiv_cf_aipw needs unit-level folds (dev_unit_folds)", call. = FALSE)
  V <- attr(fold, "V")
  n <- nrow(sample)
  if (length(fold) != n) stop("fold must carry one label per sampled unit", call. = FALSE)
  rows <- fold_rows(fold)
  fits <- vector("list", V)
  fit_rows <- score_rows <- vector("list", V)
  Q1 <- Q0 <- g <- rep(NA_real_, n)
  clip_hits <- 0L
  for (v in seq_len(V)) {
    score <- rows[[v]]
    train <- sort(unlist(rows[-v], use.names = FALSE))
    fits[[v]] <- if ("fold" %in% names(formals(learner))) learner(sample[train, , drop = FALSE], fold = v) else
      learner(sample[train, , drop = FALSE])
    nd <- sample[score, , drop = FALSE]
    Q1[score] <- fits[[v]]$Q(1, nd)
    Q0[score] <- fits[[v]]$Q(0, nd)
    g[score] <- fits[[v]]$g(nd)
    fit_rows[[v]] <- train
    score_rows[[v]] <- score
    if (!is.null(fits[[v]]$raw)) {
      q0c <- fits[[v]]$conventions$q0; dc <- fits[[v]]$conventions$delta
      rq1 <- fits[[v]]$raw$Q(1, nd); rq0 <- fits[[v]]$raw$Q(0, nd); rg <- fits[[v]]$raw$g(nd)
      clip_hits <- clip_hits + sum(rq1 < q0c | rq1 > 1 - q0c | rq0 < q0c | rq0 > 1 - q0c | rg < dc | rg > 1 - dc)
    }
  }
  phi <- aipw_phi(sample$Y, sample$A, Q1, Q0, g)
  N_hat <- sum(sample$w)
  N_hat_v <- vapply(rows, function(r) sum(sample$w[r]), numeric(1))
  omega <- N_hat_v / N_hat
  fold_means <- vapply(seq_len(V), function(v) { r <- rows[[v]]; sum(sample$w[r] * phi[r]) / N_hat_v[v] }, numeric(1))
  psi <- sum(omega * fold_means)
  uc <- hajek_uc(phi, sample$w, sample$h, sample$k)
  diag <- do.call(rbind, lapply(seq_len(V), function(v) {
    d <- fits[[v]]$diagnostics
    if (is.null(d)) NULL else cbind(fold = v, d)
  }))
  list(arm = arm, psi = psi, se2 = uc$se2, se = sqrt(uc$se2), psi_pooled = uc$psi,
       N_hat = N_hat, N_hat_v = N_hat_v, omega = omega, fold_means = fold_means,
       phi = phi, Q1 = Q1, Q0 = Q0, g = g, fits = fits, fit_rows = fit_rows,
       score_rows = score_rows, fold = fold, V = V, clip_hits = clip_hits,
       diagnostics = diag, label = EST_DEV_LABEL)
}

# no-CF AIPW: one fit on the whole sample (learner(sample, fold = 0) for stochastic learners)
# scores that same sample; psi is the Hajek mean of phi_hat, SE the UC estimator (U8 applies).
nocf_aipw <- function(sample, learner, arm) {
  fit <- if ("fold" %in% names(formals(learner))) learner(sample, fold = 0L) else learner(sample)
  Q1 <- fit$Q(1, sample); Q0 <- fit$Q(0, sample); g <- fit$g(sample)
  clip_hits <- 0L
  if (!is.null(fit$raw)) {
    q0c <- fit$conventions$q0; dc <- fit$conventions$delta
    rq1 <- fit$raw$Q(1, sample); rq0 <- fit$raw$Q(0, sample); rg <- fit$raw$g(sample)
    clip_hits <- sum(rq1 < q0c | rq1 > 1 - q0c | rq0 < q0c | rq0 > 1 - q0c | rg < dc | rg > 1 - dc)
  }
  phi <- aipw_phi(sample$Y, sample$A, Q1, Q0, g)
  uc <- hajek_uc(phi, sample$w, sample$h, sample$k)
  list(arm = arm, psi = uc$psi, se2 = uc$se2, se = sqrt(uc$se2), N_hat = sum(sample$w),
       phi = phi, Q1 = Q1, Q0 = Q0, g = g, fits = list(fit), clip_hits = clip_hits,
       diagnostics = if (is.null(fit$diagnostics)) NULL else cbind(fold = 0L, fit$diagnostics),
       label = EST_DEV_LABEL)
}

# U4's unit analogue and U5-false-by-design for individual-CF.
check_indcf <- function(sample, cf) {
  f <- as.vector(cf$fold); V <- cf$V
  one_label <- length(f) == nrow(sample) && !anyNA(f) && all(f %in% seq_len(V))
  sizes <- table(sample$h, factor(f, levels = seq_len(V)))
  balanced <- all(apply(sizes, 1L, function(x) max(x) - min(x) <= 1L))
  key <- paste(sample$h, sample$k)
  shared <- vapply(seq_len(V), function(v)
    length(intersect(unique(key[cf$fit_rows[[v]]]), unique(key[cf$score_rows[[v]]]))), integer(1))
  scored_once <- identical(sort(unlist(cf$score_rows, use.names = FALSE)), seq_len(nrow(sample)))
  .check_rows(
    check = c("U4 (unit analogue): every unit carries exactly one fold label in 1..V",
              "U4 (unit analogue): per-stratum unit fold sizes differ by at most one (TC-E6)",
              "U5 (unit analogue): scored sets partition the sample (every unit scored once)",
              sprintf("U5 false by design (:343): fold %d training and scored PSU sets overlap (%d PSUs shared)", seq_len(V), shared)),
    value = c(as.numeric(one_label), as.numeric(balanced), as.numeric(scored_once), shared),
    reference = c(1, 1, 1, rep(NA_real_, V)), tolerance = c(0, 0, 0, rep(NA_real_, V)),
    pass = c(one_label, balanced, scored_once, shared > 0L))
}

ARM_RA_NOCF <- "Ra-AIPW-nocf"         # E3 R-a no-CF negative control
