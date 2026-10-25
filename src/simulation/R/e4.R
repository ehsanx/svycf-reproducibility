# E4 -- the model-layer-free design-layer stress (the C7 isolation witness): a deterministic array
# with one unit per PSU (U2 vacuous), 1 x 2000 draws per rung, G6 limb (b) at delta_L = 0.50.
#
# No model is posited or invoked: no DGP, no eta0, no A, no Y, no learner. The census PSU score
# totals t_hk are written down directly; the oracle score of the single unit of PSU (h, k) IS
# t_hk. MANY-STRATA ONLY: H in {25, 50, 100, 200, 400}, m_h = 2, J = 20 (f = 0.1), m = 2H.
#
# Unit layer: M_hk = 1 for every PSU, so w = M_pop,h / m_h = 10 and N_hat = N = H J exactly in
# every draw; U2's target is identically zero and U2 is recorded as VACUOUS, never as a pass.
# There is no stage-2 randomisation (n0 = 1 = M_hk): E4 is NOT a two-stage variance check and is
# never graded on C3a. Implementation check only (C7): no PASS or RED-FLAG.
#
# The design-regularity gate, census algebra at every rung:
#   uniform bound   max |t_hk| <= C_t = 1.25 (C_t = A (max s + max |a|), a consequence of A = 1)
#   Lindeberg       L = max_h Var_des(W_h) / sum_h Var_des(W_h), with Var_des(W_h) = (1 - f_h) B_h
#                   the design variance of stratum h's linearised contribution; a limit (-> 0), so
#                   it is reported per rung (with H L), not graded at one rung
#   cond:nd         c = 0.01 <= m sigma2_m,wor <= C = 100 (gate-validated constants)
# Dominant-stratum control (G6 limb (b)): s(1) = sqrt(H); L >= delta_L = 0.50 at EVERY rung
# (no monotonicity). It also breaches the uniform bound: stated, not hidden.

E4_H <- c(25L, 50L, 100L, 200L, 400L)
E4_J <- 20L; E4_MH <- 2L
E4_A <- 1; E4_CT <- 1.25; E4_DELTA_L <- 0.50; E4_ND_C <- c(lo = 0.01, hi = 100)

e4_s <- function(h) 1 / 2 + (1 / 2) * ((h - 1) %% 4) / 3            # in [0.5, 1.0]
e4_a <- function(h) (1 / 4) * (-1)^h

# The array: t_hk = A s(h) cos(2 pi k / J) + A a(h); the control sets s(1) = sqrt(H).
e4_array <- function(H, J = E4_J, A = E4_A, dominant = FALSE) {
  h <- rep(seq_len(H), each = J); k <- rep(seq_len(J), times = H)
  s <- e4_s(h)
  if (dominant) s[h == 1L] <- sqrt(H)
  data.frame(h = as.integer(h), k = as.integer(k), t = A * s * cos(2 * pi * k / J) + A * e4_a(h))
}

# Population in the engine's shape (one unit per PSU): psu (h, k, M = 1), units (h, k, i = 1,
# phi0 = t). Usable by sampling_frame()/draw_sample() with n0 = 1.
e4_population <- function(arr, label = "E4 deterministic array (§12.7), R11 unit layer") {
  psu <- data.frame(h = arr$h, k = arr$k, M = 1L)
  units <- data.frame(h = arr$h, k = arr$k, i = 1L, phi0 = arr$t)
  list(psu = psu, units = units,
       meta = list(H = max(arr$h), J = max(arr$k), N = nrow(units), label = label,
                   hash = rng_store_hash(units, c("h", "k", "i", "phi0"), "E4-array-v1")))
}

# First-stage census algebra (no stage 2): with z = t - Psi_N per unit (= PSU total),
# B_h = M_pop,h^2 S2_t,h / (m_h N^2), V_h = 0; sigma2_m = sum B_h; sigma2_m,wor = sum (1 - f_h) B_h.
e4_census <- function(pop, m_h = E4_MH) {
  u <- pop$units
  if (any(pop$psu$M != 1L)) stop("e4_census: the unit layer has one unit per PSU (M_hk = 1)", call. = FALSE)
  N <- nrow(u); Psi_N <- mean(u$phi0)
  strata <- sort(unique(u$h)); sh <- match(u$h, strata)
  M_pop_h <- as.vector(table(factor(u$h, levels = strata)))
  z <- u$phi0 - Psi_N
  zbar <- as.vector(rowsum(z, sh, reorder = TRUE)) / M_pop_h
  S2_t <- as.vector(rowsum((z - zbar[sh])^2, sh, reorder = TRUE)) / (M_pop_h - 1)
  B_h <- M_pop_h^2 * S2_t / (m_h * N^2)
  f_h <- m_h / M_pop_h
  var_w <- (1 - f_h) * B_h
  m <- m_h * length(strata)
  list(N = N, Psi_N = Psi_N, m = m, per_stratum = data.frame(h = strata, M_pop_h = M_pop_h, m_h = m_h, f_h = f_h,
                                                             S2_t = S2_t, B_h = B_h, var_W = var_w),
       sigma2_m = sum(B_h), sigma2_m_wor = sum(var_w), Delta_m = sum(f_h * B_h),
       lindeberg = max(var_w) / sum(var_w), m_sigma2_wor = m * sum(var_w), max_abs_t = max(abs(u$phi0)))
}

# The regularity gate and G6 limb (b) on one rung: one row per assertion. required = TRUE rows
# are the frozen numeric criteria; the Lindeberg ratio of the regular array is reported.
e4_gate_rows <- function(H, reg, dom) {
  data.frame(
    H = H,
    check = c(sprintf("E4 uniform bound: max |t_hk| <= C_t = %.2f (regular array)", E4_CT),
              sprintf("E4 cond:nd: %.2f <= m sigma2_m,wor <= %g (regular array)", E4_ND_C[["lo"]], E4_ND_C[["hi"]]),
              "E4 Lindeberg ratio max_h Var(W_h) / sum_h Var(W_h) (regular array; a limit -> 0, reported with H x ratio)",
              sprintf("G6 limb (b): dominant-stratum control's ratio >= delta_L = %.2f", E4_DELTA_L),
              sprintf("dominant-stratum control breaches the uniform bound C_t = %.2f (stated in §12.7)", E4_CT)),
    value = c(reg$max_abs_t, reg$m_sigma2_wor, reg$lindeberg, dom$lindeberg, dom$max_abs_t),
    reference = c(E4_CT, NA, NA, E4_DELTA_L, E4_CT),
    pass = c(reg$max_abs_t <= E4_CT, reg$m_sigma2_wor >= E4_ND_C[["lo"]] && reg$m_sigma2_wor <= E4_ND_C[["hi"]],
             TRUE, dom$lindeberg >= E4_DELTA_L, dom$max_abs_t > E4_CT),
    required = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE)
}

# Kolmogorov distance of a sample from N(0, 1) (C1's distributional QQ scalar; reported).
ks_normal <- function(z) {
  z <- sort(z); n <- length(z); p <- stats::pnorm(z)
  max(c(seq_len(n) / n - p, p - (seq_len(n) - 1) / n))
}
