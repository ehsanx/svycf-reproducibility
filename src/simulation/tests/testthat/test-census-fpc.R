# Tests of census_fpc_summary() (R/census.R): the E2 FPC-sweep census algebra.
# Pure algebra on census_targets(); SYNTHETIC root, a_cpsu = 1.25 explicit.
# Checks: sigma2_m, sigma2_m_wor, Delta_m, U8 (negative-limb deficit),
# (C3b overlay), (H = 50, m_h = 2, J in {100, 20, 8, 4}).

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R",
            "dgp.R", "census.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

ROOT <- "SYNTHETIC-DEV-E2E"
CO <- dgp_coefficients(1.25)
MASTER <- build_population(50, 100, CO, "FPCTEST", 1, ROOT)
eps64 <- 64 * .Machine$double.eps

test_that("per-stratum components add up to the census totals, and Delta_m = sigma2_m - sigma2_m_wor", {
  for (J in c(100L, 20L, 8L, 4L)) {
    cen <- census_targets(population_prefix(MASTER, 50, J), m_h = 2, n0 = 20)
    fs <- census_fpc_summary(cen)
    ps <- fs$per_stratum
    expect_true(all(ps$f_h == 2 / J))
    expect_equal(fs$fpc_bar, 1 - 2 / J, tolerance = 1e-15)
    expect_lte(abs(sum(ps$sigma2_m_h) - cen$sigma2_m), eps64 * cen$sigma2_m)
    expect_lte(abs(sum(ps$sigma2_m_wor_h) - cen$sigma2_m_wor), eps64 * cen$sigma2_m)
    expect_lte(abs(sum(ps$Delta_h) - cen$Delta_m), eps64 * cen$sigma2_m)
    expect_lte(fs$identity_gap, eps64 * cen$sigma2_m)
    expect_gt(fs$sigma_ratio, 1)                                 # f > 0 and B_h > 0 here
    expect_equal(fs$Delta_over_sigma2_m, 1 - cen$sigma2_m_wor / cen$sigma2_m, tolerance = 1e-12)
    # U8 negative limb: a blanket single-stage fpc targets sum_h (1 - f_h) sigma2_m,h,
    # which is sigma2_m_wor minus the census deficit sum_h f_h V_h.
    expect_equal(sum((1 - ps$f_h) * ps$sigma2_m_h), cen$sigma2_m_wor - fs$single_stage_fpc_deficit, tolerance = 1e-12)
    expect_true(fs$beta_pooled > 0 && fs$beta_pooled < 1)
  }
})

test_that("the C3b overlay is the nominal level when sigma_m = sigma_m,wor and exceeds it otherwise", {
  ps <- data.frame(h = 1:3, M_pop_h = 10, m_h = 2, f_h = 0, B_h = c(1, 2, 3) * 1e-4, V_h = c(3, 2, 1) * 1e-4)
  cen0 <- list(per_stratum = ps, sigma2_m = sum(ps$B_h + ps$V_h), sigma2_m_wor = sum(ps$B_h + ps$V_h), Delta_m = 0)
  fs0 <- census_fpc_summary(cen0)
  expect_equal(unname(fs0$overlay), c(0.80, 0.90, 0.95), tolerance = 1e-12)
  expect_identical(names(fs0$overlay), c("0.80", "0.90", "0.95"))
  ps$f_h <- 0.5
  cen1 <- list(per_stratum = ps, sigma2_m = sum(ps$B_h + ps$V_h), sigma2_m_wor = sum(0.5 * ps$B_h + ps$V_h),
               Delta_m = sum(0.5 * ps$B_h))
  fs1 <- census_fpc_summary(cen1)
  expect_true(all(fs1$overlay > c(0.80, 0.90, 0.95)))
  expect_equal(unname(fs1$overlay[3]), 2 * stats::pnorm(stats::qnorm(0.975) * sqrt(cen1$sigma2_m / cen1$sigma2_m_wor)) - 1)
})
