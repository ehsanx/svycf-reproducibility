# Exact-enumeration regression test for the DEV calibration gate: the census-centred
# LINEARISED score with KNOWN N (lin_var_designs, R/variance.R). On a tiny two-stage population
# every SRS-WOR sample is enumerated with its exact probability; the design expectations of the
# three declared designs' variance estimators must equal the census targets EXACTLY (to rounding):
# sigma2_m, sigma2_m_wor and sigma2_m_wor - sum_h f_h V_h, and E[sum w u] = 0. The Hajek UC
# estimator (estimated psi and N_hat in its residuals) does NOT reproduce sigma2_m exactly -- the
# reason the gate moved off it (review of d570258).

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "census.R", "variance.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

# Fixture: H = 2 strata x J = 3 PSUs, sizes 2..3, m_h = 2 (f = 2/3), n0 = 2.
PSU <- data.frame(h = rep(1:2, each = 3), k = rep(1:3, 2), M = c(2L, 3L, 3L, 3L, 2L, 3L))
PHI <- c(0.3, -0.7, 1.9, 0.2, -1.1, 2.4, 0.8, 1.5, -0.4, 0.05, 3.1, -2.2, 0.6, 1.7, -0.9, 0.35)
idx <- rep(seq_len(nrow(PSU)), PSU$M)
UNITS <- data.frame(h = PSU$h[idx], k = PSU$k[idx], i = sequence(PSU$M), phi0 = PHI, Q1 = PHI, Q0 = 0)
POP <- list(psu = PSU, units = UNITS, meta = list(label = "enumeration fixture"))
M_H <- 2L; N0 <- 2L; J <- 3L

enumerate <- function() {
  pairs <- utils::combn(J, M_H, simplify = FALSE)
  acc <- c(p = 0, total = 0, uc = 0, wor = 0, fpc1 = 0, hajek_uc = 0)
  for (a in pairs) for (b in pairs) {
    sel <- rbind(PSU[PSU$h == 1 & PSU$k %in% a, ], PSU[PSU$h == 2 & PSU$k %in% b, ])
    subs <- lapply(sel$M, function(Mk) utils::combn(Mk, N0, simplify = FALSE))
    grid <- expand.grid(lapply(subs, seq_along))
    p1 <- (1 / choose(J, M_H))^2
    for (g in seq_len(nrow(grid))) {
      parts <- lapply(seq_len(nrow(sel)), function(j) {
        u <- UNITS[UNITS$h == sel$h[j] & UNITS$k == sel$k[j], ][subs[[j]][[grid[g, j]]], ]
        u$M_k <- sel$M[j]; u$M_pop_h <- J; u$w <- (J / M_H) * (sel$M[j] / N0); u
      })
      s <- do.call(rbind, parts)
      p <- p1 * prod(1 / choose(sel$M, N0))
      lv <- lin_var_designs(s$phi0, mean(UNITS$phi0), nrow(UNITS), s$w, s$h, s$k, s$i, s$M_pop_h, s$M_k)
      acc <- acc + p * c(1, lv$total, lv$se2_uc, lv$se2_wor, lv$se2_fpc1, hajek_uc(s$phi0, s$w, s$h, s$k)$se2)
    }
  }
  acc
}
E <- enumerate()
CEN <- census_targets(POP, m_h = M_H, n0 = N0)
FS <- census_fpc_summary(CEN)

test_that("the enumeration is complete: probabilities sum to 1", {
  expect_equal(unname(E["p"]), 1, tolerance = 1e-14)
})

test_that("linearised known-N estimators reproduce the three census targets EXACTLY in expectation", {
  expect_lt(abs(E[["total"]]), 1e-15)                                         # E[sum w u] = 0
  expect_equal(E[["uc"]], CEN$sigma2_m, tolerance = 1e-12)                     # U8 single-stage limb
  expect_equal(E[["wor"]], CEN$sigma2_m_wor, tolerance = 1e-12)                # U8 two-stage limb
  expect_equal(E[["fpc1"]], CEN$sigma2_m_wor - FS$single_stage_fpc_deficit, tolerance = 1e-12)  # negative limb
  expect_gt(FS$single_stage_fpc_deficit, 0)
})

test_that("the Hajek UC estimator does not: its finite-sample expectation differs from sigma2_m", {
  expect_gt(abs(E[["hajek_uc"]] / CEN$sigma2_m - 1), 1e-6)
})
