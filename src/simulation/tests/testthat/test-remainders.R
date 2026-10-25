# T1/T2 closed form and the U9 gate (R/remainders.R). Synthetic fixtures only:
# fitted nuisances are explicit functions; nothing is drawn from a production
# population. a_cpsu = 0 is an earlier coefficient table; 1.25 is the study's value.

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

# Small-but-accurate grids keep the suite fast; the check script uses the defaults.
ctx_for <- local({
  cache <- list()
  function(a, H = 5L) {
    key <- paste(a, H)
    if (is.null(cache[[key]]))
      cache[[key]] <<- p0_context(H, dgp_coefficients(a),
                                  primary = list(rule = "gauss-legendre", n_c = 6L, n_psu = 2L),
                                  second = list(rule = "clenshaw-curtis", n_c = 7L, n_psu = 3L))
    cache[[key]]
  }
})

truth_fit <- function(co, gh_n = 64L) {
  gh <- gauss_hermite(gh_n)
  list(Q = function(a, nd) Q0_marginal(a, nd, nd$KB, co, gh), g = function(nd) g0_marginal(nd, nd$KB, co, gh))
}
smooth_fit <- function(co) {
  list(Q = function(a, nd) plogis(lp_Q(a, nd, nd$KB, co) + 0.3 * nd$C1 - 0.2 * nd$Cpsu + 0.1),
       g = function(nd) plogis(0.8 * lp_g(nd, nd$KB, co) - 0.2 * nd$C2 * nd$C4))
}

test_that("quadrature rules are exact to their polynomial degree", {
  gl <- gauss_legendre(5); expect_equal(sum(gl$w * gl$x^8), 2 / 9, tolerance = 1e-14)
  expect_equal(sum(gl$w), 2, tolerance = 1e-14)
  cc <- clenshaw_curtis(9); expect_equal(sum(cc$w * cc$x^8), 2 / 9, tolerance = 1e-14)
  sub <- gauss_legendre(4, 0.25, 0.75); expect_equal(sum(sub$w * sub$x^3), (0.75^4 - 0.25^4) / 4, tolerance = 1e-14)
  gh <- gauss_hermite(20); z <- sqrt(2) * gh$x
  expect_equal(sum(gh$w * z^4) / sqrt(pi), 3, tolerance = 1e-12)
  expect_equal(sum(gh$w * z^6) / sqrt(pi), 15, tolerance = 1e-12)
})

test_that("size-weighted law and class shares match §12.3", {
  expect_equal(size_law_mean(), 34, tolerance = 1e-14)
  cl <- size_law_cells()
  expect_equal(range(cl$M), c(20, 48)); expect_equal(sum(cl$hi - cl$lo), 2, tolerance = 1e-14)
  expect_equal(p0_class_shares(5), c(3, 2) / 5); expect_equal(p0_class_shares(50), c(1, 1) / 2)
  g <- p0_grid(5, "gauss-legendre", 3L, 2L)
  expect_equal(sum(g$w), 1, tolerance = 1e-13)
  expect_equal(sum(g$w * g$C1), 0, tolerance = 1e-14)
  expect_equal(sum(g$w * g$C2^2), 1 / 3, tolerance = 1e-13)
  # size-biased C_psu: E_P0[M] = E[M^2]/E[M] = (65.5 + 34^2)/34 (exact rounded law)
  Mg <- round(34 + 14 * g$Cpsu)
  expect_equal(sum(g$w * Mg), (65.5 + 34^2) / 34, tolerance = 1e-12)
  expect_equal(sum(g$w[g$KB == 1]), 3 / 5, tolerance = 1e-13)
  # reweighting to another H changes only the class shares
  g50 <- p0_grid_reweight(g, 50L)
  expect_equal(sum(g50$w[g50$KB == 1]), 1 / 2, tolerance = 1e-13)
  expect_identical(g50[, 1:6], g[, 1:6])
})

test_that("Gauss-Hermite marginals agree with adaptive integration", {
  co <- dgp_coefficients(1.25); gh <- gauss_hermite(64)
  pts <- data.frame(C1 = c(-0.9, 0.3, 0.8), C2 = c(0.5, -0.7, 0.1), C3 = c(0.2, 0.9, -0.6),
                    C4 = c(-0.3, 0.4, 0.95), Cpsu = c(-0.95, 0.1, 0.99), KB = c(1L, 2L, 2L))
  for (i in seq_len(nrow(pts))) {
    p <- pts[i, ]; kb <- p$KB
    for (a in 0:1) {
      lp <- lp_Q(a, p, kb, co); s <- co$sigma_Y[kb]
      ref <- integrate(function(u) plogis(lp + u) * dnorm(u, 0, s), -Inf, Inf, rel.tol = 1e-13)$value
      expect_equal(Q0_marginal(a, p, kb, co, gh), ref, tolerance = 1e-10)
    }
    lp <- lp_g(p, kb, co); s <- co$sigma_A[kb]
    ref <- integrate(function(u) plogis(lp + u) * dnorm(u, 0, s), -Inf, Inf, rel.tol = 1e-13)$value
    expect_equal(g0_marginal(p, kb, co, gh), ref, tolerance = 1e-10)
  }
})

test_that("a_cpsu (S11) must be stated explicitly", {
  expect_error(dgp_coefficients(), "a_cpsu is required")
  expect_identical(dgp_coefficients(1.25)$a_cpsu, 1.25)
})

test_that("T2 vanishes when either nuisance is exact (double robustness)", {
  for (a in c(0, 1.25)) {
    ctx <- ctx_for(a); co <- ctx$co; tf <- truth_fit(co); sf <- smooth_fit(co)
    expect_equal(t2_closed_form(ctx, list(tf), 1)$T2, 0, tolerance = 1e-15)
    expect_equal(t2_closed_form(ctx, list(list(Q = tf$Q, g = sf$g)), 1)$T2, 0, tolerance = 1e-15)
    expect_equal(t2_closed_form(ctx, list(list(Q = sf$Q, g = tf$g)), 1)$T2, 0, tolerance = 1e-15)
    expect_lt(abs(t2_closed_form(ctx, list(sf), 1)$T2), 1e-2)
    expect_gt(abs(t2_closed_form(ctx, list(sf), 1)$T2), 1e-5)
  }
})

test_that("sign: Qhat1 = Q0_1 + d, ghat = 0.8 g0 gives T2 = -d/4 exactly", {
  # P0 f = E[(Qhat1 - Q0_1)(ghat - g0)/ghat] = d * (0.8 - 1)/0.8 = -d/4, for any g0.
  ctx <- ctx_for(1.25); co <- ctx$co; tf <- truth_fit(co); d <- 0.01
  fit <- list(Q = function(a, nd) tf$Q(a, nd) + if (a == 1) d else 0, g = function(nd) 0.8 * tf$g(nd))
  expect_equal(t2_closed_form(ctx, list(fit), 1)$T2, -d / 4, tolerance = 1e-14)
  expect_equal(p0_f_bruteforce(fit, ctx), -d / 4, tolerance = 1e-9)
  # The propensity factor written as (g_{a,0} - ghat_a), with + for a = 1 and - for a = 0, has the opposite sign.
  g <- ctx$primary; tr <- ctx$primary_truth
  Q1h <- fit$Q(1, g); Q0h <- fit$Q(0, g); gh <- fit$g(g)
  literal <- sum(g$w * ((tr$g - gh) / gh * (Q1h - tr$Q1) - ((1 - tr$g) - (1 - gh)) / (1 - gh) * (Q0h - tr$Q0)))
  expect_equal(literal, +d / 4, tolerance = 1e-14)
})

test_that("closed form matches brute-force P0 f_v within U9's 1e-6, per fold", {
  for (a in c(0, 1.25)) {
    ctx <- ctx_for(a); co <- ctx$co
    fits <- list(smooth_fit(co),
                 list(Q = function(a, nd) plogis(0.7 * lp_Q(a, nd, nd$KB, co) + 0.5 * nd$C3^2 - 0.1),
                      g = function(nd) plogis(lp_g(nd, nd$KB, co) + 0.4 * sin(2 * nd$C1) + 0.3 * nd$Cpsu^2)))
    xc <- u9_cross_check(ctx, fits)
    expect_true(all(xc$pass), info = paste(format(xc$abs_diff), collapse = " "))
  }
})

test_that("T2 is the omega-weighted sum of per-fold remainders; omega is checked", {
  ctx <- ctx_for(0); co <- ctx$co; tf <- truth_fit(co); sf <- smooth_fit(co)
  fits <- list(sf, list(Q = sf$Q, g = tf$g), list(Q = function(a, nd) pmin(sf$Q(a, nd) + 0.02, 0.99), g = sf$g))
  om <- c(0.30, 0.45, 0.25)
  r <- t2_closed_form(ctx, fits, om)
  expect_equal(r$T2, sum(om * r$per_fold), tolerance = 1e-16)
  expect_equal(r$per_fold[2], 0, tolerance = 1e-15)
  expect_error(t2_closed_form(ctx, fits, c(0.3, 0.45, 0.2)), "sum to 1")
  expect_error(t2_closed_form(ctx, fits, c(0.5, 0.5)))
  expect_error(t2_closed_form(ctx, fits, c(-0.1, 0.85, 0.25)))
})

test_that("clip_fit applies the frozen constants and T2 sees post-clip values", {
  raw <- list(Q = function(a, nd) ifelse(nd$C1 > 0, 1, 0),
              g = function(nd) ifelse(nd$C1 > 0, 0.9999, 0.0001))
  cf <- clip_fit(raw)
  nd <- data.frame(C1 = c(-0.5, 0.5), C2 = 0, C3 = 0, C4 = 0, Cpsu = 0, KB = 1L)
  expect_true(all(cf$g(nd) >= 0.01 & cf$g(nd) <= 0.99))
  expect_true(all(cf$Q(0, nd) >= 0.005 & cf$Q(1, nd) <= 0.995))
  ctx <- ctx_for(0)
  expect_error(t2_closed_form(ctx, list(list(Q = raw$Q, g = function(nd) rep(1, nrow(nd)))), 1), "clip it first")
  expect_true(is.finite(t2_closed_form(ctx, list(cf), 1)$T2))
})

test_that("U9 limb 1: Psi on the grid, GH64 vs GH128, second grid", {
  for (a in c(0, 1.25)) {
    u <- u9_psi(ctx_for(a))
    expect_true(all(u$pass), info = paste(format(u$abs_diff), collapse = " "))
  }
  # Psi depends on H only through the K_B shares; even H share one value.
  p50 <- u9_psi(p0_context_reweight(ctx_for(0), 50L))$value[1]
  p800 <- u9_psi(p0_context_reweight(ctx_for(0), 800L))$value[1]
  expect_identical(p50, p800)
})

test_that("T1 by subtraction and the T1 + T2 restatement", {
  expect_equal(t1_by_subtraction(0.20, 0.15, 0.01), 0.04)
  t2 <- -0.0031; t1 <- t1_by_subtraction(0.173, 0.161, t2)
  expect_equal(t1 + t2, 0.173 - 0.161, tolerance = 1e-15)
})

test_that("u9_select picks ~20 deterministic (rep, fold) pairs covering every fold", {
  s <- u9_select(400, 2L)
  expect_equal(nrow(s), 20L); expect_identical(s, u9_select(400, 2L))
  expect_setequal(s$fold, 1:2); expect_equal(range(s$rep), c(1, 400))
  expect_equal(nrow(u9_select(5, 2L)), 5L)
})

test_that("strict gate: only all-finite, all-within-tolerance results pass", {
  expect_identical(check_pass(c(1, 1, NA, NaN, Inf, 1), c(1, 1.1, 1, 1, Inf, 1), c(0, 0.05, 1, 1, 1, NA)),
                   c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
  ok <- data.frame(pass = c(TRUE, TRUE))
  expect_identical(u9_gate_status(ok), 0L)
  expect_identical(u9_gate_status(data.frame(pass = c(TRUE, FALSE))), 1L)
  expect_identical(u9_gate_status(data.frame(pass = c(TRUE, NA))), 1L)
  expect_identical(u9_gate_status(data.frame(pass = logical())), 1L)   # nothing checked is not a pass
})

test_that("grid-sizing target is 1% of Q_0.5(|T2|)", {
  expect_true(t2_grid_error_ok(-0.0031, -0.00309, 0.002)$ok)
  expect_false(t2_grid_error_ok(-0.0031, -0.0030, 0.002)$ok)
})
