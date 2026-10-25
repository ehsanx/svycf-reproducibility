# DEV / PROVISIONAL tests of the R-b random-forest learner (R/learner_rb.R)
# and of R-b PSU-CF AIPW with its paired TMLE on the estimator-test fixture
# (SYNTHETIC root, a_cpsu = 1.25 passed explicitly). No forest T1/T2
# is computed. The prediction seed (C-RB-8) keeps the ambient RNG untouched; tested below.

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R",
            "folds.R", "dgp.R", "census.R", "sampling.R", "variance.R", "learner_ra.R", "estimators.R",
            "tmle.R", "learner_rb.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

ROOT <- "SYNTHETIC-DEV-E2E"
CO <- dgp_coefficients(1.25)   # explicit, never defaulted
CELL <- list(regime = "many", H = 25, m_h = 2, J = 20, m = 50, f = "1/10", e2c = "none")
POP <- build_population(25, 20, CO, "MSTEST", 1, ROOT)
KEY <- list(cell = CELL, population_id = 1, rep_id = 1)
S <- draw_sample(POP, 2, 20, KEY, ROOT)
FOLD <- dev_psu_folds(S, KEY, 2, ROOT)
TRAIN <- S[fold_rows(FOLD)[[1]], , drop = FALSE]

test_that("rb_seeds: the keyed level-4 stream, three seeds, fold- and rep-specific, SYNTHETIC only", {
  s1 <- rb_seeds(ROOT, CELL, 1, 1, 1)
  expect_type(s1, "integer"); expect_length(s1, 3L)
  expect_true(all(s1 >= 1L & s1 <= .Machine$integer.max))
  expect_identical(rb_seeds(ROOT, CELL, 1, 1, 1), s1)
  expect_false(identical(rb_seeds(ROOT, CELL, 1, 1, 2), s1))
  expect_false(identical(rb_seeds(ROOT, CELL, 1, 2, 1), s1))
  expect_error(rb_seeds("DEV", CELL, 1, 1, 1), "SYNTHETIC")
  withr::with_seed(3, {
    s0 <- get(".Random.seed", envir = globalenv())
    rb_seeds(ROOT, CELL, 1, 1, 1)
    expect_identical(get(".Random.seed", envir = globalenv()), s0)
  })
})

FIT <- fit_rb(TRAIN, rb_seeds(ROOT, CELL, 1, 1, 1))

test_that("fit_rb uses the frozen §12.11 settings on the six covariates, Q by arm", {
  expect_identical(FIT$conventions$num.trees, 500L)
  expect_identical(FIT$conventions$mtry, 3L)
  expect_identical(FIT$conventions$min.node.size, 25L)
  expect_true(FIT$conventions$replace)
  expect_identical(FIT$conventions$respect.unordered.factors, "order")
  expect_identical(FIT$diagnostics$component, c("Q1", "Q0", "g"))
  expect_identical(FIT$diagnostics$n_train, c(sum(TRAIN$A == 1), sum(TRAIN$A == 0), nrow(TRAIN)))
  expect_true(all(is.finite(FIT$diagnostics$oob_mse)))
})

test_that("predictions: raw in [0, 1], clipped into [q0, 1-q0] and [delta, 1-delta], reproducible", {
  nd <- S[1:200, , drop = FALSE]
  r1 <- FIT$raw$Q(1, nd); r0 <- FIT$raw$Q(0, nd); rg <- FIT$raw$g(nd)
  expect_true(all(r1 >= 0 & r1 <= 1 & r0 >= 0 & r0 <= 1 & rg >= 0 & rg <= 1))
  q1 <- FIT$Q(1, nd); g <- FIT$g(nd)
  expect_true(all(q1 >= 0.005 & q1 <= 0.995 & g >= 0.01 & g <= 0.99))
  expect_false(identical(r1, r0))
  again <- fit_rb(TRAIN, rb_seeds(ROOT, CELL, 1, 1, 1))
  expect_identical(again$raw$Q(1, nd), r1); expect_identical(again$raw$g(nd), rg)
  other <- fit_rb(TRAIN, rb_seeds(ROOT, CELL, 1, 1, 2))
  expect_false(identical(other$raw$g(nd), rg))
})

test_that("OI-RB-3 fixed: fitting and every prediction leave the ambient RNG untouched", {
  nd <- S[1:200, , drop = FALSE]
  genv <- globalenv()
  withr::with_seed(11, {
    s0 <- get(".Random.seed", envir = genv)
    f <- fit_rb(TRAIN, rb_seeds(ROOT, CELL, 1, 1, 1))
    invisible(list(f$Q(1, nd), f$Q(0, nd), f$g(nd), f$raw$Q(1, nd), f$raw$g(nd)))
    expect_identical(get(".Random.seed", envir = genv), s0)
  })
  # the sensitivity control: ranger's own predict() with seed = NULL does move it (the leak)
  withr::with_seed(11, {
    s0 <- get(".Random.seed", envir = genv)
    invisible(stats::predict(FIT$forests$g, data = .rb_x(nd), num.threads = 1L))
    expect_false(identical(get(".Random.seed", envir = genv), s0))
  })
  expect_error(.rb_predict(FIT$forests$g, .rb_x(nd), NULL), "explicit prediction seed")
})

test_that("OI-RB-3 fix leaves predictions unchanged: explicit training seed == seed = NULL, identical()", {
  chk <- rb_prediction_check(FIT, S)
  expect_identical(chk$component, c("Q1", "Q0", "g"))
  expect_identical(chk$n, rep(nrow(S), 3L))
  expect_true(all(chk$identical))
  expect_identical(chk$max_abs_diff, c(0, 0, 0))
  # an unrelated explicit seed gives the same regression-forest predictions as well
  expect_identical(stats::predict(FIT$forests$g, data = .rb_x(S), seed = 1L, num.threads = 1L)$predictions,
                   FIT$raw$g(S))
  expect_s3_class(FIT$forests$Q1, "ranger"); expect_equal(FIT$forests$Q1$num.trees, 500)
  # the check's own seed = NULL calls are scoped
  withr::with_seed(5, {
    s0 <- get(".Random.seed", envir = globalenv())
    rb_prediction_check(FIT, S[1:50, , drop = FALSE])
    expect_identical(get(".Random.seed", envir = globalenv()), s0)
  })
})

test_that("a training set without one treatment arm is a hard error (a learner failure)", {
  expect_error(fit_rb(TRAIN[TRAIN$A == 0, ], rb_seeds(ROOT, CELL, 1, 1, 1)), "no training units")
  expect_error(fit_rb(TRAIN, 1:2), "three seeds")
})

test_that("R-b PSU-CF AIPW passes U3/U4/U5/U7 and its paired TMLE passes §12.9", {
  learner_b <- function(train, fold) fit_rb(train, rb_seeds(ROOT, CELL, 1, 1, fold))
  cf <- psu_cf_aipw(S, FOLD, learner_b, arm = ARM_RB_AIPW)
  expect_identical(cf$arm, ARM_RB_AIPW)
  # each fold's forest used its own fold-keyed seeds
  expect_identical(cf$fits[[1]]$diagnostics$seed, rb_seeds(ROOT, CELL, 1, 1, 1))
  expect_identical(cf$fits[[2]]$diagnostics$seed, rb_seeds(ROOT, CELL, 1, 1, 2))
  or <- oracle_aipw(S)
  expect_true(all(check_u3(S$w, cf$fold, cf$omega, cf$N_hat_v)$pass))
  expect_true(all(check_u4_u5(S, cf)$pass))
  expect_true(check_u7(S, cf$fold, cf$fits, cf$psi, or$psi)$pass)
  tm <- psu_cf_tmle(S, cf, 0.005, 0.01, arm = ARM_RB_TMLE)
  expect_identical(tm$arm, ARM_RB_TMLE)
  ch <- check_tmle(S, cf, tm)
  expect_true(all(ch$pass[1:6]))
  expect_true(is.finite(ch$value[7]))
})
