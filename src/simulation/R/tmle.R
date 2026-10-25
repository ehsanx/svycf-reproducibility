# The PSU-cross-fitted TMLE arm, paired with the PSU-CF AIPW arm
# of R/estimators.R, and the conformance checks that bind it.
# The method is the paper's TMLE variant (Corollary on TMLE; its proof in Supplement S1,
# pooled score).
#
# It shares the AIPW arm's settings: the random-number contract, a_cpsu = 1.25 and R-a's
# configuration (R/learner_ra.R header); the runner's run mode labels every output.
#
# Requires, sourced beforehand: R/remainders.R (aipw_phi, t2_closed_form),
# R/variance.R (hajek_uc) and R/estimators.R (psu_cf_aipw, .check_rows, .fold_V,
# fold_rows). Nothing here sources anything.
#
# DEFINITION (the paper's TMLE variant):
#   initial   the AIPW arm's held-out fits, Q_{-v} clipped to [q0, 1-q0] and g_{-v}
#             truncated to [delta, 1-delta] BEFORE targeting (condition (v)).
#   tilt      ONE pooled logistic fluctuation common to all folds, holding g fixed:
#             Q*(a,c) = expit{logit Q_{-v}(a,c) + eps' H(a,c)}, H = (a/g, (1-a)/(1-g)),
#             eps solving the pooled survey-weighted score P_n^w[H (Y - Q*_A)] = 0,
#             i.e. weighted logistic regression of Y on H with offset logit Q_{-v}(A,C).
#             No post-fluctuation clip (T3).
#   estimate  psi_TMLE = sum_v omega_v (P_n^w)_v phi(.; Q*, g) = P_n^w phi(.; eta_cf)
#             (cor:tmle), the same fold-weighted Hajek form as the AIPW arm.
#   SE        UC/no-FPC hajek_uc() on phi(.; Q*, g) (U8; Supplement S1).
#   T2, T1    T2^TMLE = sum_v omega_v P0{phi(.; Q*_v, g_v) - phi(.; eta0)} by
#             t2_closed_form() on the targeted fold fits; T1^TMLE by subtraction.
#
# TECHNICAL CONVENTIONS (labelled; none changes a graded object unless flagged in
#   TC-T1 Clever covariate H = (A/g, (1-A)/(1-g)): the two-column choice that the
#         contract's "exactly 2 columns" requires. Because A(1-A) = 0,
#         the Gram matrix is diagonal and its lambda_min is min(P_n^w[A/g^2],
#         P_n^w[(1-A)/(1-g)^2]).
#   TC-T2 The targeting call is stats::glm.fit(quasibinomial, offset, no intercept)
#         with weights w / mean(w) and glm.control(epsilon = 1e-14, maxit = 100);
#         weight scaling does not change eps. The conformance refit is a separate
#         Newton-Raphson on the score equations, written without glm and with
#         Hajek-normalised weights; Q* must agree to 1e-8 at every unit and arm.
#   TC-T3 Pairing: the TMLE consumes the AIPW arm's own out-of-fold nuisances (same
#         sample, same folds, same fits). R-a is a deterministic function of its
#         training set (TC-E4), so sharing the fit is exact here. The forest arm R-b
#         shares its AIPW arm's fits in the same way (C-RB-6, R/learner_rb.R).
#   TC-T4 Arm label "Ra-TMLE" (TC-E5).

ARM_RA_TMLE <- "Ra-TMLE"            # TC-T4
TMLE_C_H <- 0.01                    # Gram floor c_H on the Hajek-weighted mean Gram
TMLE_CONFORMANCE_TOL <- 1e-8
TMLE_SCORE_TOL <- 1e-9              # targeted score equations, Hajek-normalised (convergence check)

# Clever covariate (TC-T1): exactly two columns.
tmle_clever <- function(A, g) cbind(H1 = A / g, H0 = (1 - A) / (1 - g))

# Fluctuated outcome regression for both arms (no post-fluctuation clip, T3).
tmle_qstar <- function(Q1, Q0, g, eps) {
  list(Q1 = stats::plogis(stats::qlogis(Q1) + eps[1] / g),
       Q0 = stats::plogis(stats::qlogis(Q0) + eps[2] / (1 - g)))
}

# The targeting call (TC-T2): one pooled weighted logistic fit of Y on H with offset
# logit Q_initial(A, C). Returns the 2-vector eps and the design actually used.
tmle_fluctuate <- function(Y, A, Q1, Q0, g, w) {
  H <- tmle_clever(A, g)
  off <- stats::qlogis(ifelse(A == 1, Q1, Q0))
  fit <- suppressWarnings(stats::glm.fit(
    x = H, y = Y, weights = w / mean(w), offset = off, family = stats::quasibinomial(),
    intercept = FALSE, control = stats::glm.control(epsilon = 1e-14, maxit = 100)))
  list(eps = unname(fit$coefficients), converged = isTRUE(fit$converged), iter = fit$iter,
       H = H, offset = off)
}

# independent refit (TC-T2): Newton-Raphson on the pooled weighted score, with
# its own H, offset and normalised weights; shares no code with tmle_fluctuate().
tmle_refit_newton <- function(Y, A, Q1, Q0, g, w, tol = 1e-13, maxit = 200L) {
  h1 <- A / g
  h0 <- (1 - A) / (1 - g)
  q <- ifelse(A == 1, Q1, Q0)
  off <- log(q / (1 - q))
  W <- w / sum(w)
  # Weighted logistic loss, computed stably: log(1 + exp(x)) = max(x, 0) + log1p(exp(-|x|)).
  loss <- function(e) {
    eta <- off + e[1] * h1 + e[2] * h0
    sum(W * (pmax(eta, 0) + log1p(exp(-abs(eta))) - Y * eta))
  }
  e <- c(0, 0)
  step <- c(Inf, Inf)
  it <- 0L
  while (it < maxit && max(abs(step)) >= tol) {
    it <- it + 1L
    p <- 1 / (1 + exp(-(off + e[1] * h1 + e[2] * h0)))
    s <- c(sum(W * h1 * (Y - p)), sum(W * h0 * (Y - p)))
    v <- W * p * (1 - p)
    info <- matrix(c(sum(v * h1 * h1), sum(v * h1 * h0), sum(v * h1 * h0), sum(v * h0 * h0)), 2L)
    step <- solve(info, s)
    # Damped Newton: halve the step until the convex loss does not increase.
    l0 <- loss(e)
    while (loss(e + step) > l0 + 1e-15 * abs(l0) && max(abs(step)) >= tol) step <- step / 2
    e <- e + step
  }
  Q1s <- 1 / (1 + exp(-(log(Q1 / (1 - Q1)) + e[1] / g)))
  Q0s <- 1 / (1 + exp(-(log(Q0 / (1 - Q0)) + e[2] / (1 - g))))
  list(eps = e, iter = it, converged = max(abs(step)) < tol, Q1 = Q1s, Q0 = Q0s)
}

# Hajek-weighted mean Gram of the clever covariates actually used:
# G = sum_i w_i H_i H_i' / N_hat, no per-column rescaling.
tmle_gram <- function(H, w) {
  G <- crossprod(H * sqrt(w)) / sum(w)
  list(G = G, lambda_min = min(eigen(G, symmetric = TRUE, only.values = TRUE)$values))
}

# Targeted fold fits for T2^TMLE: fold v's clipped initial fit composed with the
# pooled eps; g is unchanged.
tmle_targeted_fits <- function(fits, eps) {
  lapply(fits, function(f) {
    force(f)
    list(Q = function(a, newdata) {
      q <- f$Q(a, newdata); gg <- f$g(newdata)
      stats::plogis(stats::qlogis(q) + if (a == 1) eps[1] / gg else eps[2] / (1 - gg))
    }, g = f$g, eps = eps)
  })
}

# PSU-CF TMLE paired with the AIPW arm `cf` (the output of psu_cf_aipw() on this
# sample). q0 and delta are the frozen clip constants.
psu_cf_tmle <- function(sample, cf, q0, delta, arm = ARM_RA_TMLE) {
  n <- nrow(sample)
  if (length(cf$Q1) != n || length(cf$g) != n) stop("cf must be the AIPW arm on this sample", call. = FALSE)
  tg <- tmle_fluctuate(sample$Y, sample$A, cf$Q1, cf$Q0, cf$g, sample$w)
  qs <- tmle_qstar(cf$Q1, cf$Q0, cf$g, tg$eps)
  phi <- aipw_phi(sample$Y, sample$A, qs$Q1, qs$Q0, cf$g)
  rows <- fold_rows(cf$fold)
  V <- length(rows)
  N_hat <- sum(sample$w)
  N_hat_v <- vapply(rows, function(r) sum(sample$w[r]), numeric(1))
  omega <- N_hat_v / N_hat
  fold_means <- vapply(seq_len(V), function(v) {
    r <- rows[[v]]
    sum(sample$w[r] * phi[r]) / N_hat_v[v]
  }, numeric(1))
  psi <- sum(omega * fold_means)
  plugin <- sum(sample$w * (qs$Q1 - qs$Q0)) / N_hat
  uc <- hajek_uc(phi, sample$w, sample$h, sample$k)
  gram <- tmle_gram(tg$H, sample$w)
  list(arm = arm, psi = psi, se2 = uc$se2, se = sqrt(uc$se2), psi_pooled = uc$psi, psi_plugin = plugin,
       N_hat = N_hat, N_hat_v = N_hat_v, omega = omega, fold_means = fold_means,
       eps = tg$eps, converged = tg$converged, iter = tg$iter, H = tg$H, offset = tg$offset,
       Q1_star = qs$Q1, Q0_star = qs$Q0, g = cf$g, Q1_init = cf$Q1, Q0_init = cf$Q0,
       phi = phi, gram = gram$G, lambda_min = gram$lambda_min,
       targeted_fits = tmle_targeted_fits(cf$fits, tg$eps), q0 = q0, delta = delta,
       label = EST_DEV_LABEL)
}

# conformance and targeting checks for one replicate. `tm` is psu_cf_tmle()'s
# output, `cf` the AIPW arm it was paired with. Gram: a premise, recorded every time;
# a breach makes C6 UNGRADABLE (premise), not a build failure, so the
# caller records it as not required.
check_tmle <- function(sample, cf, tm) {
  n <- nrow(sample)
  lab <- as.vector(cf$fold)
  # (1) initial nuisances are the clipped out-of-fold values, recomputed from scratch.
  Q1r <- Q0r <- gr <- rep(NA_real_, n)
  for (v in seq_along(cf$fits)) {
    i <- which(lab == v)
    nd <- sample[i, , drop = FALSE]
    Q1r[i] <- cf$fits[[v]]$Q(1, nd); Q0r[i] <- cf$fits[[v]]$Q(0, nd); gr[i] <- cf$fits[[v]]$g(nd)
  }
  same_init <- identical(Q1r, tm$Q1_init) && identical(Q0r, tm$Q0_init) && identical(gr, tm$g)
  in_box <- all(tm$Q1_init >= tm$q0 & tm$Q1_init <= 1 - tm$q0 & tm$Q0_init >= tm$q0 &
                  tm$Q0_init <= 1 - tm$q0 & tm$g >= tm$delta & tm$g <= 1 - tm$delta)
  # (2) one pooled tilt: a 2-column design, one 2-vector eps for every fold.
  eps_shared <- all(vapply(tm$targeted_fits, function(f) identical(f$eps, tm$eps), logical(1)))
  pooled <- ncol(tm$H) == 2L && length(tm$eps) == 2L && eps_shared
  # (3) conformance: Q* against the independent Newton refit, both arms, every unit.
  ref <- tmle_refit_newton(sample$Y, sample$A, tm$Q1_init, tm$Q0_init, tm$g, sample$w)
  conf <- max(abs(c(tm$Q1_star - ref$Q1, tm$Q0_star - ref$Q0)))
  # (4) targeted score equations on the Hajek scale.
  QA <- ifelse(sample$A == 1, tm$Q1_star, tm$Q0_star)
  score <- colSums(sample$w * tm$H * (sample$Y - QA)) / sum(sample$w)
  # (5) one-step = plug-in once the score equations hold (cor:tmle).
  plug_tol <- 1e-8 * max(1, abs(tm$psi))
  .check_rows(
    check = c("TMLE §12.9 initial Q clipped to [q0,1-q0], g truncated to [delta,1-delta], equal to the out-of-fold fits (recomputed)",
              "TMLE §12.9 single pooled fluctuation: 2-column clever covariate, one 2-vector eps for all folds",
              "TMLE §12.9 targeting fit converged (glm.fit)",
              "TMLE §12.9 conformance: Q* = independent pooled weighted logistic refit (1e-8, both arms, every unit)",
              "TMLE score equations P_n^w[H (Y - Q*)] = 0 (Hajek scale, 1e-9)",
              "TMLE one-step psi = plug-in P_n^w[Q*1 - Q*0] (info, 1e-8)",
              "TMLE premise (T2, §12.9): Gram lambda_min >= c_H = 0.01 (breach -> C6 UNGRADABLE, not a build failure)"),
    value = c(as.numeric(same_init && in_box), as.numeric(pooled), as.numeric(tm$converged), conf,
              max(abs(score)), tm$psi, tm$lambda_min),
    reference = c(1, 1, 1, 0, 0, tm$psi_plugin, TMLE_C_H),
    tolerance = c(0, 0, 0, TMLE_CONFORMANCE_TOL, TMLE_SCORE_TOL, plug_tol, NA_real_),
    pass = c(same_init && in_box, pooled, tm$converged,
             is.finite(conf) && conf <= TMLE_CONFORMANCE_TOL,
             all(is.finite(score)) && max(abs(score)) <= TMLE_SCORE_TOL,
             is.finite(tm$psi - tm$psi_plugin) && abs(tm$psi - tm$psi_plugin) <= plug_tol,
             is.finite(tm$lambda_min) && tm$lambda_min >= TMLE_C_H))
}
