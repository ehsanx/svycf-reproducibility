# R/t2_refine.R: the T2 refinement of fold fits whose clipped or truncated nuisance puts a kink
# in the integrand. Synthetic fits and small grids only (quadrature; no population, sample or fit), with exact
# answers where the integral has one. Also documents the two limitations: the trigger sees a fit only at the primary
# nodes, and a step between rungs is an error estimate, not a bound.

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "quad_rules_io.R", "t2_refine.R"))
  source(test_path("..", "..", "R", f), local = TRUE)
engine_real <- normalizePath(test_path("..", ".."), winslash = "/")
co <- dgp_coefficients(1.25)
small_ctx <- function(H = 5L) p0_context(H, co, primary = list(rule = "gauss-legendre", n_c = 4L, n_psu = 2L),
                                         second = list(rule = "clenshaw-curtis", n_c = 5L, n_psu = 3L))
const_truth <- function(nd) list(Q1 = rep(0.3, nrow(nd)), Q0 = rep(0.3, nrow(nd)), g = rep(0.5, nrow(nd)))
const_ctx <- function() { ctx <- small_ctx(); ctx$truth_fun <- const_truth; ctx$primary_truth <- const_truth(ctx$primary); ctx }
q0 <- 0.005; dl <- 0.01
# A fit that depends on C1 only: arm 1 is f1(C1), arm 0 equals the constant truth (so its term is 0), g is 0.6 against
# a truth of 0.5. Then P0 f = E[f1(X) - 0.3] (0.1 / 0.6), X ~ U(-1, 1): every other axis integrates to one.
c1_fit <- function(f1, g = 0.6) list(Q = function(a, nd) if (a == 1) f1(nd$C1) else rep(0.3, nrow(nd)),
                                     g = function(nd) rep(g, nrow(nd)), conventions = list(q0 = q0, delta = dl))
lse <- function(z) ifelse(z > 30, z, log1p(exp(z)))
# A smooth GH-based fit (not clipped anywhere): eta0 shifted on the logit scale.
gh_fit <- function(sQ, sg) list(
  Q = function(a, nd) pmin(pmax(stats::plogis(stats::qlogis(Q0_marginal(a, nd, nd$KB, co, gauss_hermite(16L))) + sQ), q0), 1 - q0),
  g = function(nd) pmin(pmax(stats::plogis(stats::qlogis(g0_marginal(nd, nd$KB, co, gauss_hermite(16L))) + sg), dl), 1 - dl),
  conventions = list(q0 = q0, delta = dl))

test_that("rung 0 is the frozen primary grid bit for bit, and an untriggered arm-replicate's T2 equals t2_closed_form()'s", {
  ctx <- small_ctx()
  fits <- list(gh_fit(0.2, -0.15), gh_fit(-0.1, 0.25))
  for (f in fits) {
    v <- t2r_closed_from_values(f$Q(1, ctx$primary), f$Q(0, ctx$primary), f$g(ctx$primary), ctx$primary, ctx$primary_truth)
    expect_identical(v, p0_f_closed(f, ctx$primary, ctx$primary_truth))
  }
  r <- t2_refined(ctx, fits, c(0.47, 0.53), T2_REFINE_POLICY_DEV)
  expect_identical(r$status, "primary")
  expect_identical(r$n_triggered, 0L)
  expect_identical(r$T2, t2_closed_form(ctx, fits, c(0.47, 0.53))$T2)
  expect_identical(r$per_fold, t2_closed_form(ctx, fits, c(0.47, 0.53))$per_fold)
})

test_that("the primary grid reshapes to its tensor (i1 fastest), so adjacent-node differences are along the axes", {
  ctx <- small_ctx(); g <- ctx$primary
  a1 <- t2r_as_tensor(g$C1, g); a4 <- t2r_as_tensor(g$C4, g); ap <- t2r_as_tensor(g$Cpsu, g); ak <- t2r_as_tensor(g$KB, g)
  expect_equal(dim(a1), c(4L, 4L, 4L, 4L, 58L, 2L))
  expect_identical(t2r_adj_max(a1, 2), 0); expect_gt(t2r_adj_max(a1, 1), 0)          # C1 varies along axis 1 only
  expect_identical(t2r_adj_max(a4, 1), 0); expect_gt(t2r_adj_max(a4, 4), 0)
  expect_identical(t2r_adj_max(ap, 1), 0); expect_gt(t2r_adj_max(ap, 5), 0)
  expect_equal(t2r_adj_max(ak, 5), 0); expect_equal(t2r_adj_max(ak, 6), 1)              # K_B is an integer column
  expect_null(t2r_as_tensor(g$C1[-1], g))                                               # wrong length: no tensor
})

test_that("the trigger: a clipped node triggers (distance 0); a smooth fit far from the bounds does not; near the bound within the margin triggers", {
  ctx <- small_ctx(); g <- ctx$primary
  far <- gh_fit(0, 0)
  t0 <- t2r_trigger(far$Q(1, g), far$Q(0, g), far$g(g), g, q0, dl)
  expect_false(t0$triggered); expect_gt(t0$dist_min, 1)
  clipped <- c1_fit(function(x) pmin(pmax(stats::plogis(-3 + 3 * x), q0), 1 - q0))
  t1 <- t2r_trigger(clipped$Q(1, g), clipped$Q(0, g), clipped$g(g), g, q0, dl)
  expect_true(t1$triggered); expect_identical(t1$dist_min, 0); expect_gt(t1$at_bound, 0)
  near <- c1_fit(function(x) stats::plogis(-4.8 + 0.3 * x))                            # never clipped, within 1 logit of q0
  t2 <- t2r_trigger(near$Q(1, g), near$Q(0, g), near$g(g), g, q0, dl)
  expect_true(t2$triggered); expect_identical(t2$at_bound, 0L); expect_lt(t2$dist_min, 1)
})

test_that("a triggered smooth fit is accepted at the first refined rung and matches the exact integral", {
  ctx <- const_ctx()
  a <- -4.8; b <- 0.3
  fit <- c1_fit(function(x) stats::plogis(a + b * x))
  exact <- ((lse(a + b) - lse(a - b)) / (2 * b) - 0.3) * (0.1 / 0.6)
  r <- t2_fold_refined(fit, ctx, T2_REFINE_POLICY_DEV)
  expect_identical(r$status, "refined"); expect_identical(r$rung, 1L)
  expect_lt(abs(r$value - exact), 1e-12)
  expect_lte(r$steps[1], T2_REFINE_POLICY_DEV$tol_conv)
})

test_that("a kink the ladder does not resolve is UNRESOLVED: T2 and T1 are NA, the rungs are kept, nothing is substituted", {
  ctx <- const_ctx()
  a <- -3; b <- 3
  fk <- c1_fit(function(x) pmin(pmax(stats::plogis(a + b * x), q0), 1 - q0))
  r <- t2_refined(ctx, list(fk, fk), c(0.5, 0.5), T2_REFINE_POLICY_DEV)
  expect_identical(r$status, "unresolved"); expect_true(is.na(r$T2)); expect_identical(r$n_unresolved, 2L)
  expect_true(is.finite(r$T2_last))
  expect_length(r$folds[[1]]$values, nrow(T2_REFINE_POLICY_DEV$ladder))                     # every rung computed and kept
  cols <- t2_refined_cols(r, psi_cf = 0.21, psi_oracle = 0.2, sig_wor = 0.01)
  expect_true(is.na(cols$T2)); expect_true(is.na(cols$T1)); expect_identical(cols$T2_status, "unresolved")
  expect_identical(cols$T2_unresolved_last, r$T2_last)
  dg <- t2_refined_diag(r, 1L, 1000001L, "Rd-AIPW-a025")
  expect_identical(nrow(dg), 2L); expect_true(all(dg$status == "unresolved")); expect_true(all(is.finite(as.matrix(dg[, c("rung0", "rung1", "rung2")]))))
  # resolved arm-replicate: T1 = (psi_cf - psi_oracle) - T2
  ok <- t2_refined(ctx, list(c1_fit(function(x) stats::plogis(-4.8 + 0.3 * x))), 1, T2_REFINE_POLICY_DEV)
  oc <- t2_refined_cols(ok, 0.21, 0.2, 0.01)
  expect_identical(oc$T2_status, "refined"); expect_identical(oc$T1, (0.21 - 0.2) - ok$T2)
})

test_that("LIMITATION 1: a step is an estimate -- on a one-dimensional kink Gauss-Legendre errors oscillate with the node count", {
  ctx <- const_ctx()
  a <- -3; b <- 3; xs <- (stats::qlogis(q0) - a) / b
  f1 <- function(x) pmin(pmax(stats::plogis(a + b * x), q0), 1 - q0)
  exact <- (0.5 * (q0 * (xs + 1) + (lse(a + b) - lse(a + b * xs)) / b) - 0.3) * (0.1 / 0.6)
  err <- vapply(c(4L, 5L, 6L), function(n) {
    tp <- t2r_tensor("gauss-legendre", n, 2L, 5L); tp$truth <- lapply(tp$chunks, function(ch) const_truth(t2r_chunk_data(tp, ch)$nd))
    t2r_rung_value(c1_fit(f1), tp) - exact
  }, 0)
  expect_equal(err, c(1.083e-6, 4.764e-6, 6.960e-7), tolerance = 1e-3)                    # not monotone: 4.8e-6 sits between 1.1e-6 and 7e-7
  # A looser tolerance accepts the first small step even though the accepted value is 4.8e-6 off: the reason the
  # acceptance rule halves U9's 1e-6, requires decreasing steps from the second rung on, and keeps U9's independent check.
  loose <- within(T2_REFINE_POLICY_DEV, tol_conv <- 5e-6)
  r <- t2_fold_refined(c1_fit(f1), ctx, loose)
  expect_identical(r$status, "refined"); expect_identical(r$rung, 1L)
  expect_gt(abs(r$value - exact), 4e-6)
  # the decreasing-step condition: a second step larger than the first is not convergence
  expect_false(t2r_accept(c(3.7e-6, 4.1e-6), 5e-6)); expect_true(t2r_accept(c(1.7e-5, 1.5e-6, 4.3e-7), 5e-7))
  expect_false(t2r_accept(c(1.7e-5, 1.5e-6), 5e-7)); expect_true(t2r_accept(4e-7, 5e-7))
})

test_that("LIMITATION 2: a clipped region between the primary nodes, invisible at them, does not trigger", {
  ctx <- const_ctx(); g <- ctx$primary
  dip <- function(x) pmin(pmax(0.3 - 0.3 * exp(-(x / 0.04)^2), q0), 1 - q0)             # clipped near C1 = 0; GL(4) nodes at +-0.34, +-0.86
  fit <- c1_fit(dip)
  tg <- t2r_trigger(fit$Q(1, g), fit$Q(0, g), fit$g(g), g, q0, dl)
  expect_false(tg$triggered)                                                              # every node reads 0.3: nothing to see
  ig <- function(lo, hi) stats::integrate(function(x) 0.5 * (dip(x) - 0.3) * (0.1 / 0.6), lo, hi, subdivisions = 2000L, rel.tol = 1e-13, abs.tol = 1e-16)$value
  ref <- ig(-1, -0.2) + ig(-0.2, 0) + ig(0, 0.2) + ig(0.2, 1)
  r <- t2_fold_refined(fit, ctx, T2_REFINE_POLICY_DEV)
  expect_identical(r$status, "primary")
  expect_gt(abs(r$value - ref), 1e-6)                                                     # the primary grid's value is wrong and is kept
})

test_that("a trigger passed in (the paired TMLE uses its AIPW fold's) is used instead of the fit's own", {
  ctx <- const_ctx()
  smooth <- c1_fit(function(x) stats::plogis(-1 + 0.3 * x))
  own <- t2_refined(ctx, list(smooth), 1, T2_REFINE_POLICY_DEV)
  expect_identical(own$status, "primary")
  forced <- t2_refined(ctx, list(smooth), 1, T2_REFINE_POLICY_DEV, triggers = list(list(triggered = TRUE, dist_min = 0, reach_max = 0, at_bound = 1L)))
  expect_identical(forced$status, "refined")
  expect_lt(abs(forced$T2 - own$T2), 1e-10)
})

test_that("the tensor code reproduces p0_grid() on a finer rule (summation order only)", {
  ctx <- small_ctx(); f <- gh_fit(0.3, -0.2)
  tp <- t2r_rung(new.env(parent = emptyenv()), ctx, "gauss-legendre", 5L, 2L, ctx$gh, 2e5)
  g5 <- p0_grid(5L, "gauss-legendre", 5L, 2L)
  expect_lt(abs(t2r_rung_value(f, tp, "closed") - p0_f_closed(f, g5, truth_at(g5, co, ctx$gh))), 1e-13)
  # the H-specific class shares enter the chunks exactly as p0_grid's weights do
  tp25 <- t2r_tensor("gauss-legendre", 5L, 2L, 25L)
  expect_equal(tp25$shares, p0_class_shares(25L))
})

test_that("U9 rows: an untriggered fold is round 3's row exactly; a refined fold uses the finer brute force; an unresolved fold is UNRESOLVED", {
  ctx <- small_ctx()
  smooth <- gh_fit(0.2, -0.1)
  r <- t2_refined(ctx, list(smooth, smooth), c(0.5, 0.5), T2_REFINE_POLICY_DEV)
  u <- u9_cross_check_refined(ctx, list(smooth, smooth), r, 1:2, T2_REFINE_POLICY_DEV)
  old <- u9_cross_check(ctx, list(smooth, smooth))
  one <- u[!grepl("identical nodes", u$check), ]
  expect_identical(one$check, paste0("U9 brute-force P0 f_v vs closed form (fold ", 1:2, ", 1e-6)"))
  expect_identical(one$value, old$closed_form); expect_identical(one$reference, old$brute_force); expect_identical(one$pass, old$pass)
  idn <- u[grepl("identical nodes", u$check), ]
  expect_true(all(idn$pass)); expect_true(all(idn$abs_diff <= 1e-12))
  # refined fold (forced trigger on a smooth fit): the brute force runs on the DEV policy's finer CC rule and agrees
  rf <- t2_refined(ctx, list(smooth), 1, T2_REFINE_POLICY_DEV, triggers = list(list(triggered = TRUE, dist_min = 0, reach_max = 0, at_bound = 1L)))
  ur <- u9_cross_check_refined(ctx, list(smooth), rf, 1L, T2_REFINE_POLICY_DEV)
  expect_true(any(grepl("refined fold, rung 1; brute force on CC\\(6,4\\)", ur$check)))
  expect_true(all(ur$pass))
  # unresolved fold (constant-truth kink): its 1e-6 row is UNRESOLVED and does not pass
  cctx <- const_ctx()
  fk <- c1_fit(function(x) pmin(pmax(stats::plogis(-3 + 3 * x), q0), 1 - q0))
  rk <- t2_refined(cctx, list(fk), 1, T2_REFINE_POLICY_DEV)
  uk <- u9_cross_check_refined(cctx, list(fk), rk, 1L, T2_REFINE_POLICY_DEV)
  row <- uk[grepl("UNRESOLVED", uk$check), ]
  expect_identical(nrow(row), 1L); expect_false(row$pass)
  expect_true(uk$pass[grepl("identical nodes", uk$check)])                                  # the code check still runs and passes
})

test_that("the refinement rules: a separately hashed artifact, merged after the frozen one, never overriding it", {
  QR <- quad_rules_dev(engine_real)
  RR <- t2_refine_rules(engine_real)
  expect_identical(RR$sha256, RR$manifest$sha256)
  quad_rules_activate(QR, required = TRUE)
  on.exit(quad_rules_deactivate())
  expect_error(quadrature_rule("gauss-legendre", 8L), "not in the active saved-rule artifact")
  added <- t2_refine_rules_extend(RR)
  expect_setequal(added, c("gauss-legendre:8", "gauss-legendre:3", "gauss-legendre:10", "gauss-legendre:4", "gauss-legendre:12",
                           "clenshaw-curtis:13", "clenshaw-curtis:7"))
  expect_identical(.QUAD_STORE$sha256, QR$sha256)                                           # the frozen artifact stays the active identity
  for (n in c(3L, 4L, 8L, 10L, 12L)) expect_equal(sum(gauss_legendre(n)$w), 2, tolerance = 1e-14)
  expect_equal(sum(clenshaw_curtis(13L)$w), 2, tolerance = 1e-14)
  # the frozen grids are untouched: GL(6,2) at H = 5 still equals the MANIFEST record
  gh <- QR$manifest$expected_grid_hashes
  expect_identical(p0_grid_hash(p0_grid(5L, "gauss-legendre", 6L, 2L)), gh$primary_sha256[gh$H == 5])
  bad <- RR; bad$rules[["gauss-legendre:6"]] <- list(x = c(0, 1), w = c(1, 1))
  expect_error(t2_refine_rules_extend(bad), "disagrees with the frozen artifact")
  # the production policy's cap, as documented
  expect_identical(t2_policy_max_points(T2_REFINE_POLICY), 150336 + 712704 + 2320000 + 7216128)
})
