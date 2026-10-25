# Tests of E4, the model-layer-free deterministic array
# with its unit layer (one unit per PSU; U2 VACUOUS). Pure census algebra plus one
# sampled draw per rung through the engine's own sampler.

for (f in c("rng_streams.R", "census.R", "sampling.R", "variance.R", "e4.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

ROOT <- "SYNTHETIC-DEV-E2E"

test_that("§12.7 recipe: s(h), a(h), amplitude, and the uniform bound C_t = 1.25", {
  expect_equal(e4_s(1:5), c(0.5, 2 / 3, 5 / 6, 1, 0.5), tolerance = 1e-15)
  expect_equal(e4_a(1:2), c(-0.25, 0.25))
  arr <- e4_array(25)
  expect_identical(nrow(arr), 25L * 20L)
  expect_equal(arr$t[arr$h == 4 & arr$k == 20], 1 * 1 * cos(2 * pi) + 0.25, tolerance = 1e-15)
  for (H in E4_H) expect_lte(max(abs(e4_array(H)$t)), E4_CT)
  dom <- e4_array(100, dominant = TRUE)
  expect_equal(max(abs(dom$t[dom$h == 1])), sqrt(100) + 0.25, tolerance = 1e-12)
  expect_identical(dom[dom$h > 1, ], e4_array(100)[e4_array(100)$h > 1, ])
})

test_that("census: N = HJ, Psi_N = mean a(h), cond:nd holds, Lindeberg falls like 1/H; the control is dominant at every rung", {
  for (H in E4_H) {
    reg <- e4_census(e4_population(e4_array(H)))
    dom <- e4_census(e4_population(e4_array(H, dominant = TRUE)))
    expect_identical(reg$N, H * 20L)
    expect_equal(reg$Psi_N, mean(e4_a(seq_len(H))), tolerance = 1e-14)
    expect_true(reg$m_sigma2_wor >= 0.01 && reg$m_sigma2_wor <= 100)
    expect_equal(reg$m_sigma2_wor, 0.9 * (10 / 19) * mean(e4_s(seq_len(H))^2), tolerance = 1e-12)  # 0.282895 when 4 | H (S7)
    expect_lt(reg$lindeberg, E4_DELTA_L)
    expect_gte(dom$lindeberg, E4_DELTA_L)
    expect_gt(dom$max_abs_t, E4_CT)
    g <- e4_gate_rows(H, reg, dom)
    expect_true(all(g$pass))
    expect_equal(reg$sigma2_m - reg$sigma2_m_wor, reg$Delta_m, tolerance = 1e-14)
  }
  L <- vapply(E4_H, function(H) e4_census(e4_population(e4_array(H)))$lindeberg, numeric(1))
  expect_true(all(diff(L) < 0))
  expect_lt(diff(range(L * E4_H)) / mean(L * E4_H), 0.10)                # H x ratio ~ constant
})

test_that("R11 unit layer: one sampled draw has w = M_pop,h/m_h = 10, N_hat = N exactly; U1 and U8 hold", {
  pop <- e4_population(e4_array(25))
  cell <- list(regime = "many", H = 25, m_h = 2, J = 20, m = 50, f = "1/10", e2c = "e4")
  s <- draw_sample(pop, 2, 1, list(cell = cell, population_id = 0, rep_id = 1), ROOT, psu_cols = character(0))
  expect_true(all(s$w == 10))
  expect_identical(sum(s$w), as.numeric(pop$meta$N))                      # U2's target is identically zero
  expect_true(check_u1(s, pop)$pass)
  expect_true(check_u8(s$phi0, s$w, s$h, s$k)$pass)
  expect_identical(nrow(s), 50L)
})

test_that("array hash is a stable content hash; ks_normal is a Kolmogorov distance", {
  h1 <- e4_population(e4_array(50))$meta$hash
  expect_match(h1, "^[0-9a-f]{64}$")
  expect_identical(e4_population(e4_array(50))$meta$hash, h1)
  expect_false(identical(e4_population(e4_array(50, dominant = TRUE))$meta$hash, h1))
  z <- stats::qnorm(stats::ppoints(2000))
  expect_lt(ks_normal(z), 0.001)
  expect_gt(ks_normal(z + 0.5), 0.15)
})
