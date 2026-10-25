# R-a learner (R/learner_ra.R): class-specific-variance RE model, GH-marginalised,
# clipped, prevalidated against eta0. Every sample is SYNTHETIC, drawn from
# the law by ra_synthetic_sample() under a SYNTHETIC root; a_cpsu = 1.25 is
# the study's value and is passed explicitly.

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R", "learner_ra.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

CO <- dgp_coefficients(a_cpsu = 1.25)

# Fits are cached: glmer is the slow part, and determinism is tested separately.
fit_cache <- local({
  cache <- list()
  function(n, seed) {
    key <- paste(n, seed)
    if (is.null(cache[[key]])) cache[[key]] <<- fit_ra(ra_synthetic_sample(CO, n, seed), a_cpsu_term = TRUE)
    cache[[key]]
  }
})
ctx_primary <- local({ ctx <- NULL; function() { if (is.null(ctx)) ctx <<- ra_p0_primary(2L, CO); ctx } })

test_that("the synthetic generator follows §12.3 and refuses non-SYNTHETIC roots", {
  d <- ra_synthetic_sample(CO, 30, seed = 3)
  expect_equal(nrow(d), 2 * 30 * RA_N0)
  expect_true(all(d$KB[d$h == 1] == 1) && all(d$KB[d$h == 2] == 2))   # parity,
  expect_true(all(abs(d[c("C1", "C2", "C3", "C4", "Cpsu")]) <= 1))
  expect_equal(max(tapply(d$Cpsu, paste(d$h, d$k), function(x) length(unique(x)))), 1)
  # entity streams: a smaller sample is the k <= n prefix of a larger one, same seed
  big <- ra_synthetic_sample(CO, 60, seed = 3)
  expect_identical(lapply(big[big$k <= 30, names(d)], identity), lapply(d, identity))
  expect_error(ra_synthetic_sample(CO, 5, seed = 1, root = "DEV-E2E"), "SYNTHETIC")
})

test_that("the S11 A:Cpsu term follows an explicit, required argument", {
  d <- ra_synthetic_sample(CO, 40, seed = 11)
  expect_error(fit_ra(d), "a_cpsu_term")
  expect_error(fit_ra(d, a_cpsu_term = 1), "TRUE or FALSE")
  with_t <- fit_cache(40, 11)
  without <- fit_ra(d, a_cpsu_term = FALSE)
  for (cls in c("1", "2")) {
    expect_true("A:Cpsu" %in% names(with_t$par[[cls]]$Q$beta))
    expect_false("A:Cpsu" %in% names(without$par[[cls]]$Q$beta))
    expect_false("A:Cpsu" %in% names(with_t$par[[cls]]$g$beta))
  }
})

test_that("predictions are marginal (GH64 at sigma_hat) and inside the clip bounds", {
  fit <- fit_cache(40, 11)
  gp <- ctx_primary()$primary
  q1 <- fit$Q(1, gp); q0 <- fit$Q(0, gp); g <- fit$g(gp)
  expect_true(all(is.finite(c(q1, q0, g))))
  expect_true(all(c(q1, q0) >= RA_Q_CLIP & c(q1, q0) <= 1 - RA_Q_CLIP))   # q0 = 0.005,
  expect_true(all(g >= RA_G_TRUNC & g <= 1 - RA_G_TRUNC))                 # delta = 0.01,
  # Off-support covariates force the clip to bind on both sides, and the clipped
  # values sit exactly on the bounds.
  # Rows "high" (every term, incl. C1*C3 = +81, pushes lp up) and "low" (C1*C3 = -81,
  # the rest pushes lp down), in each class.
  ext <- data.frame(C1 = c(9, -9), C2 = c(-9, 9), C3 = c(9, 9), C4 = c(-9, 9), Cpsu = c(9, -9))
  ext <- rbind(cbind(ext, KB = 1L), cbind(ext, KB = 2L))
  qe <- c(fit$Q(1, ext), fit$Q(0, ext)); ge <- fit$g(ext)
  expect_true(all(qe >= RA_Q_CLIP & qe <= 1 - RA_Q_CLIP) && all(ge >= RA_G_TRUNC & ge <= 1 - RA_G_TRUNC))
  expect_true(any(qe == RA_Q_CLIP) && any(qe == 1 - RA_Q_CLIP))
  expect_true(any(ge == RA_G_TRUNC) && any(ge == 1 - RA_G_TRUNC))
  # Marginal, not conditional-at-u=0: rebuild one point by hand with the class sigma.
  i <- which(gp$KB == 2L & gp$C1 == max(gp$C1) & gp$C3 == max(gp$C3) & gp$Cpsu == max(gp$Cpsu))[1]
  p <- fit$par[["2"]]$g; x <- gp[i, ]
  lp <- sum(p$beta * c(1, x$C1, x$C2, x$C3, x$C4, x$Cpsu, x$C1 * x$C3))
  gh <- gauss_hermite(64)
  by_hand <- sum(gh$w * plogis(lp + sqrt(2) * p$sigma * gh$x)) / sqrt(pi)
  expect_equal(fit$raw$g(x), by_hand, tolerance = 1e-12)
  expect_gt(abs(by_hand - plogis(lp)), 1e-4)
  # a newdata class the fit never saw is a hard error
  one <- ra_synthetic_sample(CO, 10, seed = 2); one <- one[one$KB == 1, ]
  expect_error(fit_ra(one, a_cpsu_term = TRUE)$g(gp[gp$KB == 2, ][1:3, ]), "not trained")
})

test_that("per-class random-intercept variances differ when the DGP's sigmas differ", {
  # The DGP's class sigmas: sigma_A = (0.50, 1.00), sigma_Y = (0.40, 0.90).
  fit <- fit_cache(300, 21)
  s <- function(cls, w) fit$par[[cls]][[w]]$sigma
  expect_gt(s("2", "g") - s("1", "g"), 0.25)
  expect_gt(s("2", "Q") - s("1", "Q"), 0.25)
  expect_equal(c(s("1", "g"), s("2", "g")), CO$sigma_A, tolerance = 0.15, scale = 1)
  expect_equal(c(s("1", "Q"), s("2", "Q")), CO$sigma_Y, tolerance = 0.15, scale = 1)
  # Control: with equal DGP sigmas the class estimates are close to each other.
  co_eq <- CO; co_eq$sigma_A <- c(0.75, 0.75); co_eq$sigma_Y <- c(0.65, 0.65)
  feq <- fit_ra(ra_synthetic_sample(co_eq, 300, seed = 21), a_cpsu_term = TRUE)
  expect_lt(abs(feq$par[["2"]]$g$sigma - feq$par[["1"]]$g$sigma), 0.15)
  expect_lt(abs(feq$par[["2"]]$Q$sigma - feq$par[["1"]]$Q$sigma), 0.15)
})

test_that("prevalidation error against eta0 shrinks from a small to a larger sample", {
  ctx <- ctx_primary()
  small <- prevalidate_ra(CO, 50, seed = 5, ctx = ctx)
  large <- prevalidate_ra(CO, 400, seed = 5, ctx = ctx)
  a <- small$errors[small$errors$KB == "all", ]; b <- large$errors[large$errors$KB == "all", ]
  expect_identical(a$nuisance, b$nuisance)
  expect_true(all(b$L2_P0 < a$L2_P0))
  expect_true(all(b$mean_abs_P0 < a$mean_abs_P0))
  expect_true(all(b$max_abs < a$max_abs))
  # 8x the PSUs: root-n would cut L2 by sqrt(8) = 2.8; require at least 1.5x.
  expect_true(all(a$L2_P0 / b$L2_P0 > 1.5))
  expect_true(all(b$L2_P0 < 0.03))
  expect_true(all(large$diagnostics$converged))
  expect_match(large$label, "implementation check, reported, not a gate")
})

test_that("fit_ra is deterministic", {
  d <- ra_synthetic_sample(CO, 40, seed = 11)
  f1 <- fit_ra(d, a_cpsu_term = TRUE); f2 <- fit_ra(d, a_cpsu_term = TRUE)
  gp <- ctx_primary()$primary
  expect_identical(f1$par, f2$par)
  expect_identical(f1$Q(1, gp), f2$Q(1, gp))
  expect_identical(f1$Q(0, gp), f2$Q(0, gp))
  expect_identical(f1$g(gp), f2$g(gp))
  expect_identical(f1$diagnostics, f2$diagnostics)
  # A function of the training SET (C-RA-7; U6 invariance limb binds R-a):
  # permuting the rows changes nothing, bit for bit.
  for (perm in list(rev(seq_len(nrow(d))), order(d$C2, -d$C4))) {
    fp <- fit_ra(d[perm, ], a_cpsu_term = TRUE)
    expect_identical(fp$par, f1$par)
    expect_identical(fp$Q(1, gp), f1$Q(1, gp))
    expect_identical(fp$g(gp), f1$g(gp))
  }
})

test_that("warnings and singular fits are recorded and the fit is kept (never dropped)", {
  # Singular: sigma = 0 in the DGP and few PSUs put sigma_hat on the boundary.
  co0 <- CO; co0$sigma_A <- c(0, 0); co0$sigma_Y <- c(0, 0)
  fs <- fit_ra(ra_synthetic_sample(co0, 15, seed = 1), a_cpsu_term = TRUE)
  sing <- fs$diagnostics[fs$diagnostics$singular, ]
  expect_gt(nrow(sing), 0L)
  expect_true(all(grepl("singular", sing$lme4_messages)) && all(!sing$clean) && all(sing$converged))
  expect_true(all(is.finite(fs$g(ctx_primary()$primary[1:50, ]))))
  # Real lme4 convergence warnings (a deliberately badly scaled predictor) are
  # recorded and the fit is kept; on_warning = "error" makes them fatal instead.
  d <- ra_synthetic_sample(CO, 30, seed = 4); d$C1 <- d$C1 * 1e3
  fw <- fit_ra(d, a_cpsu_term = TRUE)
  expect_true(all(fw$diagnostics$n_warnings > 0L))
  expect_true(all(grepl("failed to converge", fw$diagnostics$warnings)))
  expect_false(any(fw$diagnostics$converged))
  expect_true(all(is.finite(fw$Q(1, d[1:20, ]))))
  expect_error(fit_ra(d, a_cpsu_term = TRUE, on_warning = "error"), "very different scales")
  # A hard fit failure: one Nelder-Mead refit, then the classed ra_fit_failed condition naming the
  # model -- the caller records the arm FAILED, never a silent drop and never another learner.
  d$C1 <- d$C1 * 1e3
  expect_error(fit_ra(d, a_cpsu_term = TRUE), "R-a Q[K_B=1] FAILED after one Nelder-Mead refit", fixed = TRUE)
  expect_error(fit_ra(d, a_cpsu_term = TRUE), class = "ra_fit_failed")
})

test_that("an R-a fit plugs into the T2 closed form and passes U9's brute-force limb", {
  # U9's brute-force P0 f_v limb covers R-a. Small grids as in
  # test-remainders.R; the strict gate uses the defaults.
  ctx <- p0_context(2L, CO, primary = list(rule = "gauss-legendre", n_c = 6L, n_psu = 2L),
                    second = list(rule = "clenshaw-curtis", n_c = 7L, n_psu = 3L))
  fits <- list(fit_cache(40, 11), fit_cache(300, 21))
  cc <- u9_cross_check(ctx, fits)
  expect_true(all(cc$pass))
  t2 <- t2_closed_form(ctx, fits, omega = c(0.5, 0.5))
  expect_true(is.finite(t2$T2))
})
