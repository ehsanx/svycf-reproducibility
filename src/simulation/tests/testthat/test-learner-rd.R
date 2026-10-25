# DEV tests of the R-d learner (R/learner_rd.R): the frozen whole-PSU training
# throttle on top of R-a, on a fixed-strata fixture (H = 5, m = 200, m_h = 40, J = 400;
# SYNTHETIC root, a_cpsu = 1.25 passed explicitly). Implementation checks
# only: no nuisance rate is measured or claimed (G5 needs the full ladder).

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R",
            "folds.R", "dgp.R", "census.R", "sampling.R", "variance.R", "learner_ra.R", "estimators.R",
            "learner_rd.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

ROOT <- "SYNTHETIC-DEV-E2E"
CO <- dgp_coefficients(1.25)   # explicit, never defaulted
CELL <- list(regime = "fixed", H = 5, m_h = 40, J = 400, m = 200, f = "1/10", e2c = "none")
POP <- build_population(5, 400, CO, "FSTEST", 1, ROOT)
KEY <- list(cell = CELL, population_id = 1, rep_id = 1)
S <- draw_sample(POP, 40, 20, KEY, ROOT)
FOLD <- dev_psu_folds(S, KEY, 2, ROOT)
TRAIN <- S[fold_rows(FOLD)[[2]], , drop = FALSE]      # fold 1's training set = fold 2's units
psu_key <- function(d) paste(d$h, d$k)

test_that("q_m(alpha) reproduces the frozen §1.7a table exactly (layer (i)); arm names; frozen grid only", {
  for (m in RD_M_GRID) for (a in RD_ALPHAS)
    expect_identical(rd_q(m, a), as.integer(RD_Q_TABLE[as.character(m), sprintf("%.2f", a)]))
  share <- vapply(RD_M_GRID, function(m) max(vapply(RD_ALPHAS, function(a) rd_q(m, a), integer(1))) / (m / 2), numeric(1))
  expect_identical(share[1], 0.3)                         # the cap, attained only at m = 200
  expect_true(all(diff(share) < 0))                       # anchored form: the share falls with m
  expect_identical(vapply(RD_ALPHAS, rd_arm, ""), c("Rd-AIPW-a020", "Rd-AIPW-a025", "Rd-AIPW-a030", "Rd-AIPW-a035"))
  expect_error(rd_q(100, 0.25), "frozen grid")
  expect_error(rd_q(200, 0.4), "frozen")
  g3 <- rbind(rd_g3_rows(200, 40, 2), rd_g3_rows(800, 160, 2))
  expect_true(all(g3$pass))
  expect_identical(sum(g3$required), 24L)
})

test_that("the subsample: q_m whole PSUs, all from the training fold, keyed by alpha/fold/rep, RNG scoped", {
  q <- rd_q(200, 0.35)
  sub <- rd_subsample(TRAIN, q, ROOT, CELL, 1, 1, 0.35, 1)
  expect_identical(sub$q, 30L)
  expect_identical(sub$pool, 100L)                        # m/2 at V = 2
  expect_identical(nrow(unique(sub$psus)), 30L)
  expect_true(all(psu_key(sub$psus) %in% psu_key(TRAIN)))
  expect_identical(sub$rows, which(psu_key(TRAIN) %in% psu_key(sub$psus)))   # every unit of every drawn PSU
  expect_identical(rd_subsample(TRAIN, q, ROOT, CELL, 1, 1, 0.35, 1), sub)
  expect_false(identical(rd_subsample(TRAIN, q, ROOT, CELL, 1, 1, 0.30, 1)$psus, sub$psus))   # alpha in the arm key
  expect_false(identical(rd_subsample(TRAIN, q, ROOT, CELL, 1, 1, 0.35, 2)$psus, sub$psus))
  expect_false(identical(rd_subsample(TRAIN, q, ROOT, CELL, 1, 2, 0.35, 1)$psus, sub$psus))
  expect_error(rd_subsample(TRAIN, q, "DEV", CELL, 1, 1, 0.35, 1), "SYNTHETIC")
  withr::with_seed(9, {
    s0 <- get(".Random.seed", envir = globalenv())
    rd_subsample(TRAIN, q, ROOT, CELL, 1, 1, 0.35, 1)
    expect_identical(get(".Random.seed", envir = globalenv()), s0)
  })
})

test_that("no clipping: q_m above the training pool is a hard error, never a fallback to R-a", {
  expect_error(rd_subsample(TRAIN, 101, ROOT, CELL, 1, 1, 0.35, 1), "exceeds the training pool of 100 PSUs -- hard invariant")
  small <- TRAIN[TRAIN$h <= 1, , drop = FALSE]            # 20 training PSUs < q_m = 30
  expect_error(fit_rd(small, 0.25, 200, ROOT, CELL, 1, 1, 1, a_cpsu_term = TRUE), "hard invariant, no clipping")
  expect_identical(rd_subsample(TRAIN, 100, ROOT, CELL, 1, 1, 0.35, 1)$pool, 100L)   # q = pool is allowed
})

FIT <- fit_rd(TRAIN, 0.35, 200, ROOT, CELL, 1, 1, 1, a_cpsu_term = TRUE)

test_that("fit_rd is exactly R-a fitted on the drawn whole PSUs", {
  sub <- rd_subsample(TRAIN, 30, ROOT, CELL, 1, 1, 0.35, 1)
  ra <- fit_ra(TRAIN[sub$rows, , drop = FALSE], a_cpsu_term = TRUE)
  nd <- S[1:300, , drop = FALSE]
  expect_identical(FIT$Q(1, nd), ra$Q(1, nd)); expect_identical(FIT$Q(0, nd), ra$Q(0, nd))
  expect_identical(FIT$g(nd), ra$g(nd))
  expect_identical(FIT$throttle$q, 30L); expect_identical(FIT$throttle$pool, 100L)
  expect_identical(FIT$throttle$n_units_fit, length(sub$rows))
  expect_lt(FIT$throttle$n_units_fit, nrow(TRAIN))
  expect_identical(FIT$throttle$arm, "Rd-AIPW-a035")
  expect_identical(sum(FIT$diagnostics$n_psu), 60L)       # 30 PSUs, each in one class, Q and g fits
})

test_that("R-d PSU-CF AIPW passes U3/U4/U5/U7 and the throttle checks; mutants fail them", {
  learner_d <- function(train, fold) fit_rd(train, 0.35, 200, ROOT, CELL, 1, 1, fold, a_cpsu_term = TRUE)
  cf <- psu_cf_aipw(S, FOLD, learner_d, arm = rd_arm(0.35))
  expect_identical(cf$arm, "Rd-AIPW-a035")
  or <- oracle_aipw(S)
  expect_true(all(check_u3(S$w, cf$fold, cf$omega, cf$N_hat_v)$pass))
  expect_true(all(check_u4_u5(S, cf)$pass))
  expect_true(check_u7(S, cf$fold, cf$fits, cf$psi, or$psi)$pass)
  ch <- check_rd_throttle(S, cf, 0.35, 200)
  expect_identical(nrow(ch), 12L)
  expect_true(all(ch$pass))
  # mutant 1: a drawn PSU replaced by a PSU scored in that fold -> membership and disjointness fail
  bad <- cf
  scored <- unique(S[cf$score_rows[[1]], c("h", "k")])
  bad$fits[[1]]$throttle$psus[1, ] <- scored[1, ]
  chb <- check_rd_throttle(S, bad, 0.35, 200)
  expect_false(chb$pass[grepl("fold 1 : every drawn PSU is in", chb$check)])
  expect_false(chb$pass[grepl("fold 1 : no drawn PSU is scored", chb$check)])
  # mutant 2: the full training set fitted (a fallback to R-a) -> no-fallback and whole-PSU checks fail
  bad2 <- cf
  bad2$fits[[2]]$throttle$n_units_fit <- length(cf$fit_rows[[2]])
  chb2 <- check_rd_throttle(S, bad2, 0.35, 200)
  expect_false(chb2$pass[grepl("fold 2 : fewer units", chb2$check)])
  expect_false(chb2$pass[grepl("fold 2 : whole PSUs", chb2$check)])
  # mutant 3: a wrong q -> the count check fails
  expect_false(all(check_rd_throttle(S, cf, 0.20, 800)$pass))
})
