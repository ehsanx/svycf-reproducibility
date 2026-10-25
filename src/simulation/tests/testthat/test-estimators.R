# Tests of the estimator arms (R/estimators.R):
# oracle AIPW, PSU-CF AIPW, the level-3 fold draw, U3, U4/U5, U7 and T2/T1.
# Every population is a test population: the random-number conventions
# under a SYNTHETIC root; a_cpsu = 1.25 is the study's value,
# passed explicitly. U3,
# U4, U5, U6, U7, U8, U9.

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R",
            "folds.R", "dgp.R", "census.R", "sampling.R", "variance.R", "learner_ra.R", "estimators.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

ROOT <- "SYNTHETIC-DEV-E2E"
CO <- dgp_coefficients(1.25)   # explicit, never defaulted
GH <- gauss_hermite(64L)
# The slice's cell, on a test-only master id so no run's population is reused.
CELL <- list(regime = "many", H = 25, m_h = 2, J = 20, m = 50, f = "1/10", e2c = "none")
POP <- build_population(25, 20, CO, "MSTEST", 1, ROOT)
KEY <- list(cell = CELL, population_id = 1, rep_id = 1)
S <- draw_sample(POP, 2, 20, KEY, ROOT)
FOLD <- dev_psu_folds(S, KEY, 2, ROOT)
OR <- oracle_aipw(S)
# Small grids: these tests exercise the plumbing, the strict gate uses the defaults.
CTX <- p0_context(25L, CO, primary = list(rule = "gauss-legendre", n_c = 3L, n_psu = 1L),
                  second = list(rule = "clenshaw-curtis", n_c = 3L, n_psu = 2L))

# Cheap deterministic learners (no glmer). truth_learner returns eta0 itself;
# shift_learner depends continuously on its training set, so the two fold fits differ.
truth_learner <- function(train)
  list(Q = function(a, nd) Q0_marginal(a, nd, nd$KB, CO, GH), g = function(nd) g0_marginal(nd, nd$KB, CO, GH))
shift_learner <- function(train) {
  d <- mean(train$C1 * train$Y)
  clip_fit(list(Q = function(a, nd) Q0_marginal(a, nd, nd$KB, CO, GH) + 0.05 * d + 0.02 * a,
                g = function(nd) 0.9 * g0_marginal(nd, nd$KB, CO, GH) + 0.05 + 0.1 * d))
}
CF <- psu_cf_aipw(S, FOLD, shift_learner)

test_that("the fold draw is the keyed level-3 stream: pure, per-stratum, SYNTHETIC only", {
  expect_identical(attr(FOLD, "V"), 2L)
  expect_true(all(attr(FOLD, "m_vh") == 1L))                      # m_h/V = 1 in every stratum
  expect_match(attr(FOLD, "draw_convention"), "level 3 \\(ADOPTED 2026-10-11, R2")
  fk <- attr(FOLD, "fold_key")
  expect_identical(fk$fold_structure_id, FOLD_STRUCTURE_PRIMARY)
  # Scores equal rng_streams' own fold-level entity stream, stratum by stratum.
  tab <- attr(FOLD, "psu_table")
  sc <- dev_fold_draw(ROOT)(fk, tab[, c("stratum", "psu")])
  for (h in c(1, 2, 25))
    expect_identical(sc[tab$stratum == h], rng_uniforms("fold", c(fk, list(root = ROOT)), 2,
                                                        list(type = "stratum", h = h)))
  # Deterministic; fields outside the level-3 key (arm, fold, ...) are ignored.
  expect_identical(dev_psu_folds(S, c(KEY, list(arm = "Ra-AIPW", fold = 2)), 2, ROOT), FOLD)
  # The fold-structure id moves the partition (P(no stratum changes) = 2^-25).
  other <- dev_psu_folds(S, KEY, 2, ROOT, fold_structure_id = "psu-s2")
  expect_false(identical(as.vector(other), as.vector(FOLD)))
  # The caller's RNG is untouched.
  withr::with_seed(7, {
    s0 <- get(".Random.seed", envir = globalenv())
    dev_psu_folds(S, KEY, 2, ROOT)
    expect_identical(get(".Random.seed", envir = globalenv()), s0)
  })
  expect_error(dev_fold_draw("DEV-E2E"), "SYNTHETIC")
  expect_error(dev_psu_folds(S, KEY, 3, ROOT), "incompatible fold request")   # V is never reduced
})

test_that("the oracle is the Hajek mean of phi0 with its UC variance", {
  expect_identical(OR$psi, sum(S$w * S$phi0) / sum(S$w))
  expect_identical(OR$se2, hajek_uc(S$phi0, S$w, S$h, S$k)$se2)
  expect_identical(OR$arm, ARM_ORACLE)
  expect_match(OR$label, "engine estimator arms")
})

test_that("with eta0 as the learner, PSU-CF reproduces the oracle and T2 vanishes", {
  cf <- psu_cf_aipw(S, FOLD, truth_learner)
  expect_equal(cf$phi, S$phi0, tolerance = 1e-14)
  expect_lte(abs(cf$psi - OR$psi), 64 * .Machine$double.eps * abs(OR$psi))
  t2 <- t2_closed_form(CTX, cf$fits, cf$omega)
  expect_lte(abs(t2$T2), 1e-15)
  expect_true(check_u7(S, cf$fold, cf$fits, cf$psi, OR$psi)$pass)
})

test_that("psi_CF is the omega-weighted fold average; U3, U4/U5 and U7 pass for correct code", {
  expect_identical(CF$omega, CF$N_hat_v / sum(S$w))
  expect_identical(CF$psi, sum(CF$omega * CF$fold_means))
  # = the pooled Hajek mean of phi_hat in algebra (Supplement S1).
  expect_lte(abs(CF$psi - sum(S$w * CF$phi) / sum(S$w)), 64 * .Machine$double.eps * abs(CF$psi))
  expect_identical(CF$se2, hajek_uc(CF$phi, S$w, S$h, S$k)$se2)
  expect_true(all(check_u3(S$w, CF$fold, CF$omega, CF$N_hat_v)$pass))
  expect_true(all(check_u4_u5(S, CF)$pass))
  # held-out scoring: every unit is scored by the fit that did NOT train on its PSU
  for (v in 1:2) {
    i <- CF$score_rows[[v]]
    expect_identical(CF$g[i], CF$fits[[v]]$g(S[i, ]))
    expect_false(any(paste(S$h, S$k)[i] %in% paste(S$h, S$k)[CF$fit_rows[[v]]]))
  }
  u7 <- check_u7(S, CF$fold, CF$fits, CF$psi, OR$psi)
  expect_true(u7$pass)
  expect_lte(u7$abs_diff, 64 * .Machine$double.eps * abs(OR$psi))
  # T1 by subtraction restates the paired difference; U7 is the real test.
  t2 <- t2_closed_form(CTX, CF$fits, CF$omega)
  expect_true(is.finite(t2$T2) && abs(t2$T2) > 1e-6)
  t1 <- t1_by_subtraction(CF$psi, OR$psi, t2$T2)
  expect_equal(t1 + t2$T2, CF$psi - OR$psi, tolerance = 1e-15)
})

test_that("U7 and U3 can fail: wrong-fold scoring, perturbed psi_CF, perturbed masses", {
  # scoring every unit with the fit trained on its OWN fold (fits swapped)
  expect_false(check_u7(S, CF$fold, rev(CF$fits), CF$psi, OR$psi)$pass)
  # the same fits, but a relabelled partition: the identity sees the mismatch
  expect_false(check_u7(S, relabel_folds(CF$fold, c(2L, 1L)), CF$fits, CF$psi, OR$psi)$pass)
  # a psi_CF off by 1e-12 is outside 64 eps |psi_oracle|
  expect_false(check_u7(S, CF$fold, CF$fits, CF$psi + 1e-12, OR$psi)$pass)
  # masses that still sum to 1 but disagree with the separate path
  bad <- check_u3(S$w, CF$fold, CF$omega + c(1e-13, -1e-13), CF$N_hat_v)
  expect_true(bad$pass[1])
  expect_false(all(bad$pass))
  expect_false(all(check_u3(S$w, CF$fold, CF$omega * (1 + 1e-12), CF$N_hat_v)$pass))
})

test_that("U4/U5 violations are recorded as failed checks, not skipped", {
  leaky <- CF
  leaky$fit_rows[[1]] <- seq_len(nrow(S))                         # fold 1's fit saw its own PSUs
  r <- check_u4_u5(S, leaky)
  expect_false(r$pass[2])
  expect_match(r$check[2], "U5 violated")
  split_psu <- CF
  split_psu$fold[1] <- 3L - split_psu$fold[1]                    # one unit leaves its PSU's fold
  expect_false(check_u4_u5(S, split_psu)$pass[1])
  missing_unit <- CF
  missing_unit$score_rows[[1]] <- missing_unit$score_rows[[1]][-1]
  expect_false(check_u4_u5(S, missing_unit)$pass[3])
})

test_that("R-a end to end on one replicate: post-clip fits, U3/U7/U8/U9 pass, U6 invariance", {
  ra <- function(train) fit_ra(train, a_cpsu_term = TRUE)         # a_cpsu = 1.25 != 0
  cf <- psu_cf_aipw(S, FOLD, ra)
  expect_identical(cf$arm, ARM_RA_AIPW)
  expect_true(all(cf$Q1 >= RA_Q_CLIP & cf$Q1 <= 1 - RA_Q_CLIP))
  expect_true(all(cf$g >= RA_G_TRUNC & cf$g <= 1 - RA_G_TRUNC))
  expect_identical(nrow(cf$diagnostics), 8L)                      # Q and g, per class, per fold
  expect_true(all(check_u3(S$w, cf$fold, cf$omega, cf$N_hat_v)$pass))
  expect_true(all(check_u4_u5(S, cf)$pass))
  expect_true(check_u7(S, cf$fold, cf$fits, cf$psi, OR$psi)$pass)
  expect_true(check_u8(cf$phi, S$w, S$h, S$k)$pass)
  # U9 brute-force limb on the arm's own fold fits; grid orders as in test-learner-ra.R
  # (the 3-node test grid above is too coarse for the 1e-6 tolerance).
  ctx_u9 <- p0_context(25L, CO, primary = list(rule = "gauss-legendre", n_c = 6L, n_psu = 2L),
                       second = list(rule = "clenshaw-curtis", n_c = 7L, n_psu = 3L))
  xc <- u9_cross_check(ctx_u9, cf$fits)
  expect_true(all(xc$pass))
  # U6 invariance limb binds R-a: a global relabel leaves psi_CF
  # bit-identical, because R-a is a deterministic function of its training set.
  cf_rel <- psu_cf_aipw(S, relabel_folds(FOLD, c(2L, 1L)), ra)
  expect_identical(cf_rel$psi, cf$psi)
  expect_identical(cf_rel$se2, cf$se2)
})

test_that("required_status: a failed, NA or erroring required check fails the run; info never does", {
  ok <- data.frame(check = c("a", "b"), pass = c(TRUE, TRUE), required = c(TRUE, FALSE))
  expect_identical(required_status(ok), 0L)
  expect_identical(required_status(transform(ok, pass = c(TRUE, FALSE))), 0L)   # info row only
  expect_identical(required_status(transform(ok, pass = c(FALSE, TRUE))), 1L)
  expect_identical(required_status(transform(ok, pass = c(NA, TRUE))), 1L)
  expect_identical(required_status(ok[0, ]), 1L)
  expect_false(check_finite("x", list(1, c(2, NaN)))$pass)
  expect_false(check_finite("x", list(Inf))$pass)
  expect_true(check_finite("x", list(1, 2:3))$pass)
})
