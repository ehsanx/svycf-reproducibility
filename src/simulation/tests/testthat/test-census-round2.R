# Tests of the census algebra for the E2 unequal-allocation pair, the
# Psi-bridge statistic of the vanishing-fraction cell and the WR target's gap.

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R", "dgp.R", "census.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

ROOT <- "SYNTHETIC-DEV-E2E"
CO <- dgp_coefficients(1.25)
POP <- build_population(10, 20, CO, "C2TEST", 1, ROOT)

test_that("unequal-allocation pair: same population, equal total m, Delta_m invariant to floating point", {
  up <- census_unequal_pair(POP, 6L, 2L, 20)
  expect_true(up$m_equal)
  expect_identical(sum(up$m_h_hh), 5L * 6L + 5L * 2L)
  expect_identical(up$m_h_hh[up$KB == 2L], rep(6L, 5)); expect_identical(up$m_h_hl[up$KB == 2L], rep(2L, 5))
  expect_lte(up$delta_gap, 64 * .Machine$double.eps * up$high_on_high$Delta_m)
  expect_false(isTRUE(all.equal(up$high_on_high$sigma2_m, up$high_on_low$sigma2_m)))   # the contrast lives in sigma2_m
  expect_identical(up$high_on_high$Psi_N, up$high_on_low$Psi_N)
  expect_true(all(is.finite(up$contrast)))
})

test_that("bridge statistic and WR gap", {
  cen <- census_targets(population_prefix(POP, 10, 10), m_h = 2, n0 = 20)
  expect_equal(census_bridge(cen, cen$Psi_N + 0.001), 0.001 / sqrt(cen$sigma2_m_wor), tolerance = 1e-12)
  wr <- census_wr_target(cen)
  expect_equal(wr$wr_gap, sum(cen$per_stratum$B_h) / 10, tolerance = 1e-14)        # J = M_pop,h = 10
  expect_lt(wr$sigma2_m_wr, cen$sigma2_m)
})

test_that("Amendment A2: at the full J the sub-block of a frame equals its direct build (entity addressing)", {
  wide <- build_population(3, 800, CO, "FPC-A2TEST", 1, ROOT)
  expect_identical(population_prefix(wide, 2, 800), build_population(2, 800, CO, "FPC-A2TEST", 1, ROOT))
  expect_identical(population_prefix(wide, 3, 100), build_population(3, 100, CO, "FPC-A2TEST", 1, ROOT))  # the old 100-wide block
})
