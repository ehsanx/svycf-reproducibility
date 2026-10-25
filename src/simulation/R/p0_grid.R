# Deterministic quadrature grids over the generating law P_{0,H} of the unit-level
# covariates, for T2's closed form and U9.
#
# The law: C1..C4 ~ U(-1,1) i.i.d.; C_psu ~ U(-1,1) per PSU;
# PSU size M_k = round(34 + 14*C_psu) on {20..48}; K_B(h) = 1 for odd h, 2 for even h.
# T2 integrates over the SIZE-WEIGHTED C-law: a unit drawn from P0 has
# C_psu density proportional to M(C_psu), and K_B = 1 with probability ceil(H/2)/H.
# That class share is why the grid is H-specific. u^C and b_h are fixed at 0,
# so nothing else enters.
#
# The size weight is a step function of C_psu, so the C_psu axis is integrated
# cell by cell between rounding breakpoints; within a cell everything is smooth.
# Ties at the breakpoints have measure zero and do not enter an integral.

size_law_cells <- function() {
  k <- 20:48
  lo <- pmax(-1, (k - 0.5 - 34) / 14)
  hi <- pmin(1, (k + 0.5 - 34) / 14)
  data.frame(M = k, lo = lo, hi = hi)
}

# E[M_k] under C_psu ~ U(-1,1); the design fixes it at 34.000.
size_law_mean <- function() {
  cl <- size_law_cells()
  sum(cl$M * (cl$hi - cl$lo)) / 2
}

p0_class_shares <- function(H) {
  stopifnot(H >= 1L, H == round(H))
  c(ceiling(H / 2) / H, floor(H / 2) / H)
}

# rule: "gauss-legendre" (primary) or "clenshaw-curtis" (independent);
# n_c nodes on each of C1..C4, n_psu nodes inside each of the 29 C_psu cells.
p0_grid <- function(H, rule = "gauss-legendre", n_c = 8L, n_psu = 4L) {
  pc <- quadrature_rule(rule, n_c)
  cl <- size_law_cells()
  Z <- sum(cl$M * (cl$hi - cl$lo))
  psu <- do.call(rbind, lapply(seq_len(nrow(cl)), function(i) {
    r <- quadrature_rule(rule, n_psu, cl$lo[i], cl$hi[i])
    data.frame(Cpsu = r$x, w_psu = r$w * cl$M[i] / Z)
  }))
  shares <- p0_class_shares(H)
  idx <- expand.grid(i1 = seq_len(n_c), i2 = seq_len(n_c), i3 = seq_len(n_c), i4 = seq_len(n_c),
                     ip = seq_len(nrow(psu)), KB = 1:2, KEEP.OUT.ATTRS = FALSE)
  idx <- idx[shares[idx$KB] > 0, , drop = FALSE]
  wc <- pc$w / 2
  g <- data.frame(C1 = pc$x[idx$i1], C2 = pc$x[idx$i2], C3 = pc$x[idx$i3], C4 = pc$x[idx$i4],
                  Cpsu = psu$Cpsu[idx$ip], KB = idx$KB)
  g$w <- wc[idx$i1] * wc[idx$i2] * wc[idx$i3] * wc[idx$i4] * psu$w_psu[idx$ip] * shares[idx$KB]
  attr(g, "definition") <- list(H = H, rule = rule, n_c = n_c, n_psu = n_psu,
                                law = "spec :1221-1243; size-weighted C_psu; K_B by parity")
  g
}

# The same grid for another H: only the K_B class shares change, so the nodes
# are kept and the weights rescaled (H with equal shares give identical grids).
p0_grid_reweight <- function(grid, H) {
  def <- attr(grid, "definition")
  old <- p0_class_shares(def$H); new <- p0_class_shares(H)
  if (any(old == 0 & new > 0)) stop("grid lacks a class that H = ", H, " needs")
  out <- grid
  out$w <- grid$w * (new / old)[grid$KB]
  def$H <- H
  attr(out, "definition") <- def
  out
}

p0_context_reweight <- function(ctx, H) {
  ctx$H <- H
  ctx$primary <- p0_grid_reweight(ctx$primary, H); ctx$primary_hash <- p0_grid_hash(ctx$primary)
  ctx$second <- p0_grid_reweight(ctx$second, H); ctx$second_hash <- p0_grid_hash(ctx$second)
  ctx
}

# Content hash of a grid (its definition plus nodes and weights), for recording
# the frozen grid with a run. Same R version assumed; never compared across OSes.
p0_grid_hash <- function(grid) {
  digest::digest(list(attr(grid, "definition"), as.matrix(grid)), algo = "sha256")
}
