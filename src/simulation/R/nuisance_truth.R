# True marginal nuisances eta0 = (Q0, g0) of the frozen DGP, for the remainder
# machinery (R/remainders.R). Not a population builder: nothing here draws.
#
# eta0 is marginalised
# over the PSU effects u^A, u^Y by Gauss-Hermite at 64 nodes;
# 128 nodes is the quadrature-tolerance comparator.
#
# a_cpsu is the coefficient of A*C_psu in Q0's linear predictor; the study uses 1.25.
# It stays a required
# argument (never defaulted); 0 reproduces an earlier coefficient table (tests only).

dgp_coefficients <- function(a_cpsu) {
  if (missing(a_cpsu)) stop("a_cpsu is required: 1.25 = the value used in the study; 0 = an earlier coefficient table (tests only)")
  stopifnot(is.numeric(a_cpsu), length(a_cpsu) == 1L, is.finite(a_cpsu))
  list(
    # g0 table,
    g = c(icpt = 0, c1 = 0.50, c2 = -0.35, c3 = 0.45, c4 = -0.25, cpsu = 0.60, kb2 = 0.35, c1c3 = 0.40),
    sigma_A = c(0.50, 1.00),
    # Q0 table,
    Q = c(icpt = -0.15, a = 0.70, c1 = 0.45, c2 = -0.35, c3 = 0.40, c4 = -0.25, cpsu = 0.50, kb2 = 0.35,
          c1c3 = 0.40, kb2_csum = 0.25),
    sigma_Y = c(0.40, 0.90),
    a_cpsu = a_cpsu
  )
}

# C: data.frame (or list) with C1..C4, Cpsu; KB in {1, 2}.
lp_g <- function(C, KB, co) {
  b <- co$g
  kb2 <- as.numeric(KB == 2L)
  b[["icpt"]] + b[["c1"]] * C$C1 + b[["c2"]] * C$C2 + b[["c3"]] * C$C3 + b[["c4"]] * C$C4 +
    b[["cpsu"]] * C$Cpsu + b[["kb2"]] * kb2 + b[["c1c3"]] * C$C1 * C$C3
}

lp_Q <- function(a, C, KB, co) {
  b <- co$Q
  kb2 <- as.numeric(KB == 2L)
  csum <- C$C1 + C$C2 + C$C3 + C$C4
  b[["icpt"]] + b[["a"]] * a + b[["c1"]] * C$C1 + b[["c2"]] * C$C2 + b[["c3"]] * C$C3 +
    b[["c4"]] * C$C4 + b[["cpsu"]] * C$Cpsu + b[["kb2"]] * kb2 + b[["c1c3"]] * C$C1 * C$C3 +
    b[["kb2_csum"]] * kb2 * csum + co$a_cpsu * a * C$Cpsu
}

# E[expit(lp + U)], U ~ N(0, sigma^2), sigma per point; chunked so the
# n x n_gh matrix stays small.
gh_marginal <- function(lp, sigma, gh, chunk = 50000L) {
  out <- numeric(length(lp))
  z <- sqrt(2) * gh$x
  for (s in seq(1L, length(lp), by = chunk)) {
    i <- s:min(length(lp), s + chunk - 1L)
    out[i] <- drop(stats::plogis(lp[i] + outer(sigma[i], z)) %*% gh$w) / sqrt(pi)
  }
  out
}

g0_marginal <- function(C, KB, co, gh) gh_marginal(lp_g(C, KB, co), co$sigma_A[KB], gh)
Q0_marginal <- function(a, C, KB, co, gh) gh_marginal(lp_Q(a, C, KB, co), co$sigma_Y[KB], gh)

# eta0 at the points of a grid: list(Q1, Q0, g) (Q0 here is Q0(0, c), the a = 0 arm).
truth_at <- function(grid, co, gh_n = 64L) {
  gh <- gauss_hermite(gh_n)
  list(Q1 = Q0_marginal(1, grid, grid$KB, co, gh),
       Q0 = Q0_marginal(0, grid, grid$KB, co, gh),
       g  = g0_marginal(grid, grid$KB, co, gh))
}
