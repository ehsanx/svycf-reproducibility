# Quadrature rules for the P0 remainder machinery (R/remainders.R).
#
# Base R only, so the rules are part of the locked engine rather than a package
# version. Gauss rules come from the Golub-Welsch eigenproblem; Clenshaw-Curtis is
# a different rule family, used for the independent grid and the U9 brute-force
# evaluation so that a defect in one
# family cannot hide in the other.

# ---- saved-rule store --------------------------------------------
# The base rules' last bits depend on the platform's LAPACK (eigen) and libm (cos),
# so a grid hash computed from freshly built rules is not a cross-platform identity
# A saved rule artifact (R/quad_rules_io.R, grids/)
# stores every base rule a run uses, as exact hex floats with a sha256. Once a
# store is active, gauss_hermite / gauss_legendre / clenshaw_curtis take their BASE
# nodes and weights from it and scale them to [lo, hi] with the same code as a
# computed rule, so on the platform that saved the artifact the result is
# bit-identical. With required = TRUE a rule missing from the store is an error,
# never a silent recompute.
.QUAD_STORE <- new.env(parent = emptyenv())
.quad_store_get <- function(rule, n) {
  if (is.null(.QUAD_STORE$rules)) return(NULL)
  r <- .QUAD_STORE$rules[[paste0(rule, ":", n)]]
  if (is.null(r) && isTRUE(.QUAD_STORE$required))
    stop("quadrature rule ", rule, ":", n, " is not in the active saved-rule artifact ",
         .QUAD_STORE$sha256, call. = FALSE)
  r
}

# Base rules (unscaled). Gauss-Legendre and Clenshaw-Curtis bases live on [-1, 1];
# the scaling below is x = lo + half * (x_base + 1), w = half * w_base.
.gh_base <- function(n) {
  if (n == 1L) return(list(x = 0, w = sqrt(pi)))
  J <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1L) / 2)
  J[cbind(1:(n - 1L), 2:n)] <- off
  J[cbind(2:n, 1:(n - 1L))] <- off
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = sqrt(pi) * e$vectors[1L, o]^2)
}
.gl_base <- function(n) {
  k <- seq_len(n - 1L)
  off <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(1:(n - 1L), 2:n)] <- off
  J[cbind(2:n, 1:(n - 1L))] <- off
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  # 2 * v^2 then * half equals the former half * 2 * v^2 bit for bit (doubling is exact).
  list(x = e$values[o], w = 2 * e$vectors[1L, o]^2)
}
# Clenshaw-Curtis with n >= 2 points (Chebyshev extrema); explicit-weight formula
# (Trefethen, SIAM Rev. 2008, eq. 2.4 form).
.cc_base <- function(n) {
  N <- n - 1L
  j <- 0:N
  x <- cos(pi * j / N)
  w <- numeric(n)
  for (jj in j) {
    s <- 0
    for (k in seq_len(N %/% 2L)) {
      b <- if (2L * k == N) 1 else 2
      s <- s + b / (4 * k^2 - 1) * cos(2 * k * jj * pi / N)
    }
    c_j <- if (jj == 0L || jj == N) 1 else 2
    w[jj + 1L] <- c_j / N * (1 - s)
  }
  o <- order(x)
  list(x = x[o], w = w[o])
}
.scale_rule <- function(b, lo, hi) {
  half <- (hi - lo) / 2
  list(x = lo + half * (b$x + 1), w = half * b$w)
}

# Physicists' Gauss-Hermite: sum(w * f(x)) ~ integral of exp(-x^2) f(x) dx.
# E[f(Z)], Z ~ N(0, s^2), is then sum(w * f(sqrt(2) * s * x)) / sqrt(pi).
gauss_hermite <- function(n) {
  stopifnot(n >= 1L)
  st <- .quad_store_get("gauss-hermite", n)
  if (!is.null(st)) return(list(x = st$x, w = st$w))
  .gh_base(n)
}

# Gauss-Legendre on [lo, hi]; weights sum to hi - lo.
gauss_legendre <- function(n, lo = -1, hi = 1) {
  stopifnot(n >= 1L, hi > lo)
  b <- .quad_store_get("gauss-legendre", n)
  if (is.null(b)) {
    if (n == 1L) return(list(x = (lo + hi) / 2, w = hi - lo))
    b <- .gl_base(n)
  }
  .scale_rule(b, lo, hi)
}

# Clenshaw-Curtis with n >= 2 points (Chebyshev extrema) on [lo, hi]; weights
# sum to hi - lo.
clenshaw_curtis <- function(n, lo = -1, hi = 1) {
  stopifnot(n >= 2L, hi > lo)
  b <- .quad_store_get("clenshaw-curtis", n)
  if (is.null(b)) b <- .cc_base(n)
  .scale_rule(b, lo, hi)
}

quadrature_rule <- function(rule, n, lo = -1, hi = 1) {
  switch(rule,
         "gauss-legendre" = gauss_legendre(n, lo, hi),
         "clenshaw-curtis" = clenshaw_curtis(n, lo, hi),
         stop("unknown rule: ", rule))
}
