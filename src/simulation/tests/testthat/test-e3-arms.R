# Tests of the E3 pieces: R-c (fit_rc, the full-depth foil), the
# individual-CF arm with its unit-level folds, and the no-CF arm. Small fixed-strata fixture,
# SYNTHETIC root, a_cpsu = 1.25 explicit. Implementation checks only.

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R",
            "folds.R", "dgp.R", "census.R", "sampling.R", "variance.R", "learner_ra.R", "estimators.R",
            "learner_rb.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

ROOT <- "SYNTHETIC-DEV-E2E"
CO <- dgp_coefficients(1.25)
CELL <- list(regime = "fixed", H = 5, m_h = 8, J = 80, m = 40, f = "1/10", e2c = "none")
POP <- build_population(5, 80, CO, "E3TEST", 1, ROOT)
KEY <- list(cell = CELL, population_id = 1, rep_id = 1)
S <- draw_sample(POP, 8, 20, KEY, ROOT)
LRA <- function(train) fit_ra(train, a_cpsu_term = TRUE)

test_that("unit-level folds: one label per unit, balanced within strata, PSUs split, keyed and scoped", {
  uf <- dev_unit_folds(S, KEY, 2, ROOT)
  expect_identical(attr(uf, "fold_scheme"), "unit")
  expect_true(all(as.vector(uf) %in% 1:2) && length(uf) == nrow(S))
  sz <- table(S$h, as.vector(uf))
  expect_true(all(abs(sz[, 1] - sz[, 2]) <= 1))
  split_psus <- tapply(as.vector(uf), paste(S$h, S$k), function(x) length(unique(x)) > 1)
  expect_gt(sum(split_psus), 0)                                   # the leakage stress is real
  expect_identical(dev_unit_folds(S, KEY, 2, ROOT), uf)
  expect_false(identical(as.vector(dev_unit_folds(S, list(cell = CELL, population_id = 1, rep_id = 2), 2, ROOT)), as.vector(uf)))
  expect_identical(fold_rows(uf), unname(split(seq_along(uf), factor(as.vector(uf), levels = 1:2))))
  withr::with_seed(4, {
    s0 <- get(".Random.seed", envir = globalenv())
    dev_unit_folds(S, KEY, 2, ROOT)
    expect_identical(get(".Random.seed", envir = globalenv()), s0)
  })
  expect_error(assert_strict_psu_folds(uf, S$h, S$k, 2), "split across folds")   # never a PSU fold
})

test_that("individual-CF AIPW: U3 and U7 hold; U4's unit analogue holds; U5 is false by design", {
  uf <- dev_unit_folds(S, KEY, 2, ROOT)
  cf <- indiv_cf_aipw(S, uf, LRA, arm = "Ra-AIPW-indcf")
  or <- oracle_aipw(S)
  expect_true(all(check_u3(S$w, cf$fold, cf$omega, cf$N_hat_v)$pass))
  expect_true(check_u7(S, cf$fold, cf$fits, cf$psi, or$psi)$pass)
  ch <- check_indcf(S, cf)
  expect_true(all(ch$pass))
  expect_true(all(ch$value[grepl("false by design", ch$check)] > 0))
  expect_error(indiv_cf_aipw(S, dev_psu_folds(S, KEY, 2, ROOT), LRA, "x"), "unit-level folds")
})

test_that("no-CF AIPW: one full-sample fit; psi is the Hajek mean of phi_hat; U8 reproduces se2", {
  nc <- nocf_aipw(S, LRA, ARM_RA_NOCF)
  phi <- aipw_phi(S$Y, S$A, nc$Q1, nc$Q0, nc$g)
  expect_identical(nc$phi, phi)
  expect_equal(nc$psi, sum(S$w * phi) / sum(S$w), tolerance = 1e-14)
  expect_true(check_u8(nc$phi, S$w, S$h, S$k)$pass)
  expect_length(nc$fits, 1L)
  expect_identical(nc$arm, "Ra-AIPW-nocf")
})

TR1 <- S[fold_rows(dev_psu_folds(S, KEY, 2, ROOT))[[2]], , drop = FALSE]

test_that("R-c: full-depth settings, own level-4 keys, interpolates its training data", {
  sd_c <- rb_seeds(ROOT, CELL, 1, 1, 1, arm = ARM_RC_PSUCF)
  expect_false(identical(sd_c, rb_seeds(ROOT, CELL, 1, 1, 1)))                     # R-b's key differs
  expect_false(identical(rb_seeds(ROOT, CELL, 1, 1, 0, arm = ARM_RC_NOCF), sd_c))  # fold 0 = no-CF
  fit <- fit_rc(TR1, sd_c)
  expect_identical(fit$conventions$learner, "R-c")
  expect_identical(fit$conventions$min.node.size, 1L)
  expect_false(fit$conventions$replace)
  expect_true(all(is.na(fit$diagnostics$oob_mse)))                                 # no OOB units
  i1 <- TR1$A == 1
  expect_identical(fit$raw$Q(1, TR1[i1, , drop = FALSE]), as.numeric(TR1$Y[i1]))     # interpolation
  expect_identical(fit$raw$g(TR1), as.numeric(TR1$A))
  chk <- rb_prediction_check(fit, S)                                               # explicit-seed identity
  expect_true(all(chk$identical))
  withr::with_seed(8, {
    s0 <- get(".Random.seed", envir = globalenv())
    f2 <- fit_rc(TR1, sd_c); invisible(f2$Q(1, S)); invisible(f2$g(S))
    expect_identical(get(".Random.seed", envir = globalenv()), s0)
  })
})

test_that("R-c in the three E3 arms: PSU-CF passes U3-U7, no-CF runs on fold 0, individual-CF leaks", {
  lc <- function(arm) function(train, fold) fit_rc(train, rb_seeds(ROOT, CELL, 1, 1, fold, arm = arm))
  fold <- dev_psu_folds(S, KEY, 2, ROOT)
  or <- oracle_aipw(S)
  cp <- psu_cf_aipw(S, fold, lc(ARM_RC_PSUCF), arm = ARM_RC_PSUCF)
  expect_true(all(check_u3(S$w, cp$fold, cp$omega, cp$N_hat_v)$pass))
  expect_true(all(check_u4_u5(S, cp)$pass))
  expect_true(check_u7(S, cp$fold, cp$fits, cp$psi, or$psi)$pass)
  cn <- nocf_aipw(S, lc(ARM_RC_NOCF), ARM_RC_NOCF)
  expect_identical(cn$diagnostics$fold, rep(0L, 3))
  expect_identical(cn$diagnostics$seed, rb_seeds(ROOT, CELL, 1, 1, 0, arm = ARM_RC_NOCF))
  ci <- indiv_cf_aipw(S, dev_unit_folds(S, KEY, 2, ROOT), lc(ARM_RC_INDCF), ARM_RC_INDCF)
  expect_true(all(check_indcf(S, ci)$pass))
  expect_true(check_u7(S, ci$fold, ci$fits, ci$psi, or$psi)$pass)
})

# ---- the leakage-neutral sensitivity's C_psu coarsening ----------
test_that("coarsen_cpsu: quartile bins of [-1, 1] at their midpoints, boundaries closed on the left", {
  d <- data.frame(Cpsu = c(-1, -0.75, -0.5, -1e-9, 0, 0.4999, 0.5, 0.99, 1))
  expect_identical(coarsen_cpsu(d)$Cpsu, c(-0.75, -0.75, -0.25, -0.25, 0.25, 0.25, 0.75, 0.75, 0.75))
})

test_that("fit_rc_coarse: trains and predicts on coarsened C_psu; shares Rc-indcf's seeds; prediction identity", {
  seeds <- rb_seeds(ROOT, CELL, 1, 1, 1, arm = ARM_RC_INDCF)
  fc <- fit_rc_coarse(TR1, seeds)
  expect_lte(fc$n_cpsu_values_train, 4L)
  nd <- S[1:30, ]
  nd2 <- nd; nd2$Cpsu <- coarsen_cpsu(nd)$Cpsu
  expect_identical(fc$Q(1, nd), fc$Q(1, nd2))                    # only the bin of C_psu reaches the forest
  expect_identical(fc$raw$g(nd), fc$inner$raw$g(coarsen_cpsu(nd)))
  pc <- rb_prediction_check(fc$inner, coarsen_cpsu(nd))
  expect_true(all(pc$identical))
  expect_identical(fc$diagnostics$seed, seeds)
})
