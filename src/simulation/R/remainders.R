# Remainders T1, T2 and the U9 numerical gate.
#
# Definitions: f_v = phi(.; eta_hat_{-v}) - phi(.; eta0),
# T2 = sum_v omega_v P0 f_v, and T1 = (psi_CF - psi_oracle) - T2 by subtraction.
# T2 is computed in closed form on a frozen H-specific P0 grid from the post-clip
# nuisances; U9 checks the grid
# (P0 phi(.; eta0) = Psi, a second independent grid, GH64 vs GH128 <= 1e-8) and
# cross-checks P0 f_v against a brute-force evaluation (<= 1e-6).
#
# SIGN. The closed form here is the one in the frozen proofs (Supplement S1):
# arm a contributes s_a * integral (Qhat_a - Q0_a)(ghat_a - g0_a)/ghat_a dP0, with
# s_1 = +1, s_0 = -1, ghat_1 = ghat, ghat_0 = 1 - ghat. Written with the propensity
# factor as (g_{a,0} - ghat_a), the same integral is -P0 f_v.
# The definition of f_v is unambiguous, so this module follows it and the proofs; the
# brute-force check below evaluates P0 f_v from phi itself and would catch either
# sign.
#
# A fold's fitted nuisance is passed as list(Q = function(a, newdata), g = function(newdata)),
# already clipped/truncated as the arm uses it; clip_fit() applies the frozen
# constants for fixtures. newdata has columns C1..C4, Cpsu, KB.

# Efficient influence function of the ATE, without centring: P0 phi(.; eta0) = Psi.
aipw_phi <- function(y, a, Q1, Q0, g) {
  Q1 - Q0 + a / g * (y - Q1) - (1 - a) / (1 - g) * (y - Q0)
}

# Frozen constants: fitted outcome clip q0 = 0.005 and the standing-(v) propensity
# truncation delta = 0.01, both two-sided (Supplement S1).
clip_fit <- function(fit, q0 = 0.005, delta = 0.01) {
  list(Q = function(a, newdata) pmin(pmax(fit$Q(a, newdata), q0), 1 - q0),
       g = function(newdata) pmin(pmax(fit$g(newdata), delta), 1 - delta))
}

# Grids and eta0 on them, computed once per (H, DGP coefficients) and reused by
# every fold, replicate and arm. "primary" is the frozen T2 grid; "second" is the
# independent grid of another rule family; the brute-force path uses "second" with
# Gauss-Hermite at gh_bf nodes.
#
# Default orders (technical convention): GL(6,2) is 150,336
# points; on a deliberately rough smooth fixture its T2 error is ~1e-7, a 10x margin
# under U9's 1e-6, and U9 itself fails loudly if a fit needs more.
p0_context <- function(H, co,
                       primary = list(rule = "gauss-legendre", n_c = 6L, n_psu = 2L),
                       second = list(rule = "clenshaw-curtis", n_c = 9L, n_psu = 5L),
                       gh = 64L, gh_bf = 128L) {
  gp <- do.call(p0_grid, c(list(H = H), primary))
  gs <- do.call(p0_grid, c(list(H = H), second))
  list(H = H, co = co, gh = gh, gh_bf = gh_bf,
       primary = gp, primary_truth = truth_at(gp, co, gh), primary_hash = p0_grid_hash(gp),
       second = gs, second_truth = truth_at(gs, co, gh_bf), second_hash = p0_grid_hash(gs))
}

check_omega <- function(omega, n_folds) {
  stopifnot(length(omega) == n_folds, all(is.finite(omega)), all(omega >= 0))
  # U3 tolerance: fold masses sum to 1 to 64 eps.
  if (abs(sum(omega) - 1) > 64 * .Machine$double.eps) stop("fold masses must sum to 1 (U3), got ", format(sum(omega), digits = 17))
  invisible(TRUE)
}

# Per-fold P0 f_v by the closed form on a grid with its truth.
p0_f_closed <- function(fit, grid, truth) {
  Q1h <- fit$Q(1, grid); Q0h <- fit$Q(0, grid); gh <- fit$g(grid)
  if (any(gh <= 0 | gh >= 1)) stop("fitted propensity must lie strictly inside (0,1); clip it first")
  arm1 <- (Q1h - truth$Q1) * (gh - truth$g) / gh
  arm0 <- (Q0h - truth$Q0) * ((1 - gh) - (1 - truth$g)) / (1 - gh)
  sum(grid$w * (arm1 - arm0))
}

# T2 = sum_v omega_v P0 f_v on the frozen primary grid. grid = "second"
# evaluates the same closed form on the independent grid, for t2_grid_error_ok().
t2_closed_form <- function(ctx, fits, omega, grid = c("primary", "second")) {
  grid <- match.arg(grid)
  check_omega(omega, length(fits))
  per_fold <- vapply(fits, p0_f_closed, numeric(1), grid = ctx[[grid]], truth = ctx[[paste0(grid, "_truth")]])
  list(T2 = sum(omega * per_fold), per_fold = per_fold, grid_hash = ctx[[paste0(grid, "_hash")]])
}

# T1 by subtraction. A definition: T1 + T2 = psi_CF - psi_oracle tests
# nothing about T2; U7 and U9 are the real checks.
t1_by_subtraction <- function(psi_cf, psi_oracle, t2) (psi_cf - psi_oracle) - t2

# Brute-force P0 f_v: phi itself, enumerated over (A, Y) in {0,1}^2 with the DGP's
# conditional laws, on the independent grid at gh_bf nodes. It shares no algebra
# with p0_f_closed: no double-robust rearrangement, no arm sign convention. The
# DGP draws u^A and u^Y independently, so P(A=a, Y=y | c) factorises
# into the two marginal integrals.
p0_f_bruteforce <- function(fit, ctx) {
  grid <- ctx$second; tr <- ctx$second_truth
  Q1h <- fit$Q(1, grid); Q0h <- fit$Q(0, grid); gh <- fit$g(grid)
  if (any(gh <= 0 | gh >= 1)) stop("fitted propensity must lie strictly inside (0,1); clip it first")
  tot <- numeric(nrow(grid))
  for (a in 0:1) {
    pa <- if (a == 1) tr$g else 1 - tr$g
    Qa <- if (a == 1) tr$Q1 else tr$Q0
    for (y in 0:1) {
      py <- if (y == 1) Qa else 1 - Qa
      f <- aipw_phi(y, a, Q1h, Q0h, gh) - aipw_phi(y, a, tr$Q1, tr$Q0, tr$g)
      tot <- tot + pa * py * f
    }
  }
  sum(grid$w * tot)
}

# ---- U9 --------------------------------------------------------------------------

# Limb 1: P0 phi(.; eta0) = Psi, the GH64-vs-GH128 tolerance and a
# second independent grid. Psi here is the generating-law Psi of P_{0,H}.
u9_psi <- function(ctx, tol_gh = 1e-8, tol_grid = 1e-6) {
  gp <- ctx$primary; tp <- ctx$primary_truth
  psi <- sum(gp$w * (tp$Q1 - tp$Q0))
  tp128 <- truth_at(gp, ctx$co, ctx$gh_bf)
  psi128 <- sum(gp$w * (tp128$Q1 - tp128$Q0))
  psi2 <- sum(ctx$second$w * (ctx$second_truth$Q1 - ctx$second_truth$Q0))
  # P0 phi(.; eta0) by enumeration over (A, Y): equals Psi exactly in algebra.
  phi0 <- 0
  for (a in 0:1) for (y in 0:1) {
    pa <- if (a == 1) tp$g else 1 - tp$g
    Qa <- if (a == 1) tp$Q1 else tp$Q0
    py <- if (y == 1) Qa else 1 - Qa
    phi0 <- phi0 + sum(gp$w * pa * py * aipw_phi(y, a, tp$Q1, tp$Q0, tp$g))
  }
  out <- data.frame(
    check = c("GH64 vs GH128 (primary grid)", "primary vs second grid", "P0 phi(eta0) vs Psi"),
    value = c(psi, psi, phi0),
    reference = c(psi128, psi2, psi),
    abs_diff = c(abs(psi - psi128), abs(psi - psi2), abs(phi0 - psi)),
    tolerance = c(tol_gh, tol_grid, 1e-12),
    stringsAsFactors = FALSE
  )
  out$pass <- check_pass(out$value, out$reference, out$tolerance)
  out
}

# A check passes only if both sides are finite and within tolerance; NA, NaN and
# Inf fail, so an invalid computation can never read as a pass.
check_pass <- function(value, reference, tolerance) {
  ok <- is.finite(value) & is.finite(reference) & is.finite(tolerance)
  ok & abs(value - reference) <= tolerance
}

# Exit status for a strict gate: 0 only if every required check passed.
u9_gate_status <- function(results) {
  if (!nrow(results) || any(is.na(results$pass)) || !all(results$pass)) 1L else 0L
}

# Limb 2: closed-form P0 f_v vs brute force, per fold.
u9_cross_check <- function(ctx, fits, tol = 1e-6) {
  cf <- vapply(fits, p0_f_closed, numeric(1), grid = ctx$primary, truth = ctx$primary_truth)
  bf <- vapply(fits, p0_f_bruteforce, numeric(1), ctx = ctx)
  data.frame(fold = seq_along(fits), closed_form = cf, brute_force = bf,
             abs_diff = abs(cf - bf), tolerance = tol, pass = check_pass(cf, bf, tol))
}

# Which (replicate, fold) pairs get the "~20 brute-force P0 f_v cross-checks per m"
# (technical convention): k replicates evenly spaced over
# 1..n_reps, folds cycling 1..V, so every fold index is exercised.
u9_select <- function(n_reps, V, k = 20L) {
  k <- min(k, n_reps * V)
  reps <- unique(round(seq(1, n_reps, length.out = k)))
  data.frame(rep = reps, fold = ((seq_along(reps) - 1L) %% V) + 1L)
}

# Grid-sizing target: quadrature error <= 1% of the plotted Q_0.5(|T2|) at the
# largest m of an R-d trajectory. The error is estimated by primary vs second grid.
t2_grid_error_ok <- function(t2_primary, t2_second, q50_abs_t2) {
  err <- abs(t2_primary - t2_second)
  list(error = err, bound = 0.01 * q50_abs_t2, ok = check_pass(t2_primary, t2_second, 0.01 * q50_abs_t2))
}
