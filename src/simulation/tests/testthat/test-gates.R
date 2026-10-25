# Tests of the pre-run gate machinery (R/gates.R): G3's enumeration
# of every frozen cell, G1's resize rule, G5's CR2 rate regression and r_m. Synthetic
# input only; nothing here grades a claim.

for (f in c("rng_streams.R", "learner_rd.R", "gates.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

test_that("G3: every frozen cell passes (a-i), (b) and (c); a bad cell fails", {
  cells <- g3_cells()
  expect_true(all(c("E1 many-strata", "E1 fixed-strata", "E2 FPC sweep", "E2 genuinely-strict",
                    "E2 vanishing-fraction", "E3 R-c foil + R-a controls", "E3 V-sensitivity",
                    "E4 deterministic array", "E3 power-control panel", "G8 U10 cell", "E2 WR code-check",
                    "E2 unequal-allocation pair") %in% cells$cell))
  expect_identical(nrow(cells), 39L)                               # 37 + the U10 cell + the bridge rung (A2)
  term <- cells[cells$cell == "E2 vanishing-fraction (bridge-only terminal rung)", ]
  expect_identical(c(term$H, term$J), c(200L, 800L)); expect_identical(term$master, "FPC")
  expect_identical(unname(GATE_MASTERS$FPC), c(200L, 800L))        # master FPC is 200 x 800
  expect_false(any(nzchar(cells$pending)))                            # every coordinate and count is fixed
  g3 <- g3_check(cells)
  expect_true(all(g3$pass[g3$required]))
  expect_true(all(g3$pass))                                        # the (a-ii) margin holds too
  ai <- g3[grepl("^G3\\(a-i\\)", g3$check), ]
  expect_identical(nrow(ai), 20L)                                  # R-d rungs m = 200..3200 x 4 alphas
  expect_true(any(grepl("q_m = 208 <= pool 1600", ai$check)))      # alpha .35 at m = 3200
  bad <- cells[1, ]; bad$J <- 30L
  expect_false(all(g3_check(bad)$pass))                            # J beyond the MS master
  bad2 <- cells[cells$cell == "E3 V-sensitivity", ][1, ]; bad2$m_h <- "42"; bad2$V <- 5L
  expect_false(all(g3_check(bad2)$pass))                           # V does not divide m_h
})

test_that("G1 resize (§12.8): double B to a cap of 4x, then UNGRADABLE", {
  mc <- function(B) 1 / sqrt(B)
  r <- g1_resize(1, mc, 4)
  expect_identical(r$status, "GRADEABLE"); expect_identical(r$B, 16); expect_identical(r$multiplier, 4L)
  expect_identical(g1_resize(2, mc, 4)$B, 4)                         # clears at the planned count
  expect_identical(g1_resize(0.5, mc, 4)$status, "UNGRADABLE")
})

test_that("CR2 equals HC2 when every cluster is a single observation", {
  set.seed(1)
  x <- seq(1, 4, length.out = 12); X <- cbind(1, x); y <- 2 - 0.3 * x + c(0.1, -0.2, 0.05, 0.3, -0.1, 0, 0.2, -0.3, 0.15, -0.05, 0.1, -0.1)
  fit <- stats::lm.fit(X, y); e <- fit$residuals
  h <- rowSums((X %*% solve(crossprod(X))) * X)
  hc2 <- solve(crossprod(X)) %*% t(X) %*% diag(e^2 / (1 - h)) %*% X %*% solve(crossprod(X))
  expect_equal(.cr2_vcov(X, e, seq_along(y)), hc2, tolerance = 1e-12)
})

test_that("G5: recovers a known slope, PASSes the true alpha and fails a wrong one", {
  m <- rep(c(200, 400, 800, 1600, 3200), times = 10); pop <- rep(1:10, each = 5)
  # population-specific slopes -0.30 + d_p (sum d_p = 0) and intercept shifts: the pooled OLS
  # slope is exactly -0.30 on this balanced design, and CR2 sees the between-population spread.
  d <- c(0.02, -0.02, 0.01, -0.01, 0.03, -0.03, 0.0, 0.0, 0.015, -0.015)
  a <- c(0.02, -0.03, 0.01, 0.04, -0.02, 0.0, 0.03, -0.01, -0.04, 0.02)
  r <- exp(-1 + a[pop] + (-0.30 + d[pop]) * log(m))
  g <- g5_fit(r, m, pop, alpha = 0.30)
  expect_equal(g$slope, -0.30, tolerance = 1e-12)
  expect_identical(g$df, 9L)
  expect_true(g$pass)
  expect_false(g5_fit(r, m, pop, alpha = 0.10)$pass)
  expect_error(g5_fit(c(r[-1], -1), m, pop, 0.3), "finite and > 0")
})

test_that("r_m (ruling R16, proposed): zero at the truth, max over arms of the census L2 error", {
  u <- data.frame(Q1 = c(0.2, 0.4, 0.6), Q0 = c(0.1, 0.3, 0.5), g = c(0.5, 0.5, 0.5))
  truth <- list(Q = function(a, nd) if (a == 1) nd$Q1 else nd$Q0, g = function(nd) nd$g)
  expect_identical(unname(nuisance_l2(truth, u)["r"]), 0)
  off <- list(Q = function(a, nd) (if (a == 1) nd$Q1 + 0.1 else nd$Q0), g = function(nd) nd$g - 0.05)
  l2 <- nuisance_l2(off, u)
  expect_equal(unname(l2[c("Q1", "Q0", "Q", "g", "r")]), c(0.1, 0, 0.1, 0.05, 0.1), tolerance = 1e-12)
})
