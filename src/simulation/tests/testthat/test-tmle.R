# Tests of the PSU-CF TMLE arm (R/tmle.R) and its
# conformance checks. Part A uses synthetic vectors only (no population, no learner);
# Part B pairs the arm with the AIPW arm on the same small fixture as
# test-estimators.R (SYNTHETIC root, a_cpsu = 1.25 passed explicitly).
# The method is the paper's TMLE variant (Corollary on TMLE).

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R",
            "folds.R", "dgp.R", "census.R", "sampling.R", "variance.R", "learner_ra.R", "estimators.R",
            "tmle.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

Q0C <- 0.005; DEL <- 0.01           # frozen clip constants

# ---- Part A: synthetic vectors -------------------------------------------------
syn <- withr::with_seed(20261006, {
  n <- 600
  h <- rep(1:30, each = 20); k <- rep(1:60, each = 10)
  w <- stats::runif(n, 5, 40)
  C <- stats::rnorm(n)
  g <- pmin(pmax(stats::plogis(-0.3 + 0.8 * C), DEL), 1 - DEL)
  A <- stats::rbinom(n, 1, g)
  Y <- stats::rbinom(n, 1, stats::plogis(-0.2 + 0.5 * A + 0.7 * C))
  # deliberately misspecified initial Q, clipped
  Q1 <- pmin(pmax(stats::plogis(0.1 + 0.4 * C), Q0C), 1 - Q0C)
  Q0 <- pmin(pmax(stats::plogis(-0.5 + 0.2 * C), Q0C), 1 - Q0C)
  list(n = n, h = h, k = k, w = w, A = A, Y = Y, Q1 = Q1, Q0 = Q0, g = g)
})

test_that("the targeting call solves the pooled weighted score equations", {
  tg <- with(syn, tmle_fluctuate(Y, A, Q1, Q0, g, w))
  expect_true(tg$converged)
  expect_identical(dim(tg$H), c(as.integer(syn$n), 2L))                     # exactly 2 columns
  expect_length(tg$eps, 2L)
  qs <- with(syn, tmle_qstar(Q1, Q0, g, tg$eps))
  QA <- ifelse(syn$A == 1, qs$Q1, qs$Q0)
  score <- colSums(syn$w * tg$H * (syn$Y - QA)) / sum(syn$w)
  expect_lt(max(abs(score)), 1e-10)
  expect_gt(max(abs(tg$eps)), 1e-3)                                # the misspecified Q is moved
})

test_that("§12.9 conformance: glm.fit and the independent Newton refit agree to 1e-8", {
  tg <- with(syn, tmle_fluctuate(Y, A, Q1, Q0, g, w))
  qs <- with(syn, tmle_qstar(Q1, Q0, g, tg$eps))
  ref <- with(syn, tmle_refit_newton(Y, A, Q1, Q0, g, w))
  expect_true(ref$converged)
  expect_lt(max(abs(c(qs$Q1 - ref$Q1, qs$Q0 - ref$Q0))), 1e-8)
  expect_lt(max(abs(tg$eps - ref$eps)), 1e-8)
  # weight scale does not change eps (Hajek invariance)
  tg2 <- with(syn, tmle_fluctuate(Y, A, Q1, Q0, g, 1000 * w))
  expect_lt(max(abs(tg2$eps - tg$eps)), 1e-10)
})

test_that("Gram is the Hajek-weighted mean of H H', diagonal, lambda_min as stated", {
  H <- tmle_clever(syn$A, syn$g)
  gr <- tmle_gram(H, syn$w)
  N <- sum(syn$w)
  d1 <- sum(syn$w * syn$A / syn$g^2) / N
  d0 <- sum(syn$w * (1 - syn$A) / (1 - syn$g)^2) / N
  expect_equal(unname(gr$G), matrix(c(d1, 0, 0, d0), 2), tolerance = 1e-13)
  expect_equal(gr$lambda_min, min(d1, d0), tolerance = 1e-13)
  # unnormalised cross-product would be N times larger: the floor is on the mean
  expect_equal(max(eigen(crossprod(H * sqrt(syn$w)), only.values = TRUE)$values) / N,
               max(d1, d0), tolerance = 1e-12)
})

test_that("no post-fluctuation clip: Q* is the untruncated expit and may leave [q0, 1 - q0]", {
  qs <- tmle_qstar(1 - Q0C, Q0C, 0.5, c(2, -2))
  expect_identical(qs$Q1, stats::plogis(stats::qlogis(1 - Q0C) + 2 / 0.5))
  expect_gt(qs$Q1, 1 - Q0C)
  expect_lt(qs$Q0, Q0C)
})

test_that("a hard-clipped Q* fails conformance (mutation)", {
  # One unit in five sits at the clip boundary, the rest far too low: the pooled tilt
  # must raise Q1 (eps1 > 0) and so carries the boundary units outside the box.
  Q1b <- ifelse(seq_len(syn$n) %% 5L == 0L, 1 - Q0C, 0.01)
  tg <- with(syn, tmle_fluctuate(Y, A, Q1b, Q0, g, w))
  expect_gt(tg$eps[1], 0)
  qs <- with(syn, tmle_qstar(Q1b, Q0, g, tg$eps))
  ref <- with(syn, tmle_refit_newton(Y, A, Q1b, Q0, g, w))
  clipped <- pmin(pmax(qs$Q1, Q0C), 1 - Q0C)
  expect_true(any(qs$Q1 > 1 - Q0C) || any(qs$Q1 < Q0C))
  expect_gt(max(abs(clipped - ref$Q1)), 1e-8)
})

test_that("a Q that already solves the score is a fixed point: eps = 0 and TMLE = AIPW", {
  H <- tmle_clever(syn$A, syn$g)
  c1 <- sum(syn$w * H[, 1] * syn$Y) / sum(syn$w * H[, 1])
  c0 <- sum(syn$w * H[, 2] * syn$Y) / sum(syn$w * H[, 2])
  Q1c <- rep(c1, syn$n); Q0c <- rep(c0, syn$n)
  tg <- with(syn, tmle_fluctuate(Y, A, Q1c, Q0c, g, w))
  expect_lt(max(abs(tg$eps)), 1e-10)
  qs <- with(syn, tmle_qstar(Q1c, Q0c, g, tg$eps))
  expect_equal(aipw_phi(syn$Y, syn$A, qs$Q1, qs$Q0, syn$g), aipw_phi(syn$Y, syn$A, Q1c, Q0c, syn$g),
               tolerance = 1e-9)
})

# ---- Part B: paired with the AIPW arm on the estimator-test fixture -------------
ROOT <- "SYNTHETIC-DEV-E2E"
CO <- dgp_coefficients(1.25)   # explicit, never defaulted
GH <- gauss_hermite(64L)
CELL <- list(regime = "many", H = 25, m_h = 2, J = 20, m = 50, f = "1/10", e2c = "none")
POP <- build_population(25, 20, CO, "MSTEST", 1, ROOT)
KEY <- list(cell = CELL, population_id = 1, rep_id = 1)
S <- draw_sample(POP, 2, 20, KEY, ROOT)
FOLD <- dev_psu_folds(S, KEY, 2, ROOT)
OR <- oracle_aipw(S)
CTX <- p0_context(25L, CO, primary = list(rule = "gauss-legendre", n_c = 3L, n_psu = 1L),
                  second = list(rule = "clenshaw-curtis", n_c = 3L, n_psu = 2L))
shift_learner <- function(train) {
  d <- mean(train$C1 * train$Y)
  clip_fit(list(Q = function(a, nd) Q0_marginal(a, nd, nd$KB, CO, GH) + 0.05 * d + 0.02 * a,
                g = function(nd) 0.9 * g0_marginal(nd, nd$KB, CO, GH) + 0.05 + 0.1 * d))
}
CF <- psu_cf_aipw(S, FOLD, shift_learner)
TM <- psu_cf_tmle(S, CF, Q0C, DEL)

test_that("the TMLE is paired with the AIPW arm: same sample, folds and initial fits", {
  expect_identical(TM$Q1_init, CF$Q1); expect_identical(TM$Q0_init, CF$Q0); expect_identical(TM$g, CF$g)
  expect_identical(TM$omega, CF$omega); expect_identical(TM$N_hat_v, CF$N_hat_v)
  expect_identical(TM$arm, ARM_RA_TMLE)
  expect_match(TM$label, "engine estimator arms")
})

test_that("psi_TMLE is the fold-weighted Hajek mean of phi(Q*, g), with its UC SE", {
  phi <- aipw_phi(S$Y, S$A, TM$Q1_star, TM$Q0_star, TM$g)
  expect_identical(TM$phi, phi)
  expect_equal(TM$psi, sum(S$w * phi) / sum(S$w), tolerance = 1e-14)
  expect_identical(TM$se2, hajek_uc(phi, S$w, S$h, S$k)$se2)
  u8 <- check_u8(TM$phi, S$w, S$h, S$k)
  expect_true(u8$pass)
})

test_that("check_tmle: every required row passes on a correct arm", {
  ch <- check_tmle(S, CF, TM)
  expect_identical(nrow(ch), 7L)
  expect_true(all(ch$pass[1:6]))
  expect_true(is.finite(ch$value[7]))                              # Gram recorded every time
  t2 <- t2_closed_form(CTX, TM$targeted_fits, TM$omega)
  expect_true(is.finite(t2$T2))
})

test_that("check_tmle catches tampering", {
  bad <- TM; bad$Q1_init[1] <- bad$Q1_init[1] * 0.9                # initial not the out-of-fold fit
  expect_false(check_tmle(S, CF, bad)$pass[1])
  bad <- TM; bad$targeted_fits[[2]]$eps <- bad$eps + c(1e-3, 0)     # a fold-indexed tilt
  expect_false(check_tmle(S, CF, bad)$pass[2])
  bad <- TM; bad$Q1_star <- pmin(bad$Q1_star, stats::quantile(bad$Q1_star, 0.9, names = FALSE))  # a post-fluctuation clip
  expect_false(check_tmle(S, CF, bad)$pass[4])
  bad <- TM; bad$Q0_star <- bad$Q0_star + 1e-4                      # score no longer solved
  expect_false(check_tmle(S, CF, bad)$pass[5])
})

test_that("targeted fold fits apply the pooled eps to each fold's clipped initial fit", {
  nd <- S[1:5, , drop = FALSE]
  for (v in seq_along(CF$fits)) {
    f0 <- CF$fits[[v]]; f1 <- TM$targeted_fits[[v]]
    expect_identical(f1$g(nd), f0$g(nd))
    expect_equal(f1$Q(1, nd), stats::plogis(stats::qlogis(f0$Q(1, nd)) + TM$eps[1] / f0$g(nd)), tolerance = 1e-15)
    expect_equal(f1$Q(0, nd), stats::plogis(stats::qlogis(f0$Q(0, nd)) + TM$eps[2] / (1 - f0$g(nd))), tolerance = 1e-15)
  }
})
