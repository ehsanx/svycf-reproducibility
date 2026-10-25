# Census targets of a built population for one sampling cell:
# stage 1 SRS-WOR of m_h PSUs in each stratum, stage 2 SRS-WOR of n0 units in each
# sampled PSU (fixed n0 = 20). Census algebra on PSU
# totals only, no sampling.
#
# Requires: a population from build_population() / population_prefix() (R/dgp.R).
#
# TARGETS:
#   Psi_N    = P_N phi(.; eta0), the unit mean over the census, the primary target.
#   Psi_N_gf = P_N{Q0(1,C) - Q0(0,C)}, the finite-population g-formula parameter with
#              the marginal Q0 (the paper; the E3 power control keys to it).
#   The normings are those of the LINEARISED Hajek term L = N^-1 sum_i w_i z_i with
#   z = phi0 - Psi_N, w_i = (M_pop,h/m_h)(M_k/n0) (U1). With t_k the
#   census PSU total of z, S2_t,h its within-stratum variance (divisor M_pop,h - 1),
#   S2_z,k the within-PSU variance of z (divisor M_k - 1):
#     B_h = M_pop,h^2 S2_t,h / (m_h N^2)                                  between-PSU
#     V_h = M_pop,h / (m_h N^2) * sum_{k in h} M_k (M_k - n0) S2_z,k / n0   within-PSU
#   (Supplement S1: B_h = m_h S2 of the census PSU means mu_k, V_h = m_h E_K Var_2(e_K)).
#     sigma2_m     = sum_h (B_h + V_h)        = E_des[se^2], UC idealisation, no FPC
#     sigma2_m_wor = sum_h ((1-f_h) B_h + V_h) = Var_des(L) (U8; Supplement S1)
#     Delta_m      = sum_h f_h B_h            (= N^-2 sum_h M_pop,h S2_t,h)
#   f_h = m_h / M_pop,h is the first-stage fraction (Supplement S1).
#
# TECHNICAL CONVENTIONS (labelled):
#   TC-C1 m_h is a scalar or one value per stratum in h order (the unequal-allocation
#         pair allocates by K_B class; pass m_h[KB] expanded per stratum).
#   TC-C2 the census is the whole population passed in; to evaluate a cell on a larger
#         master, pass population_prefix(master, H, J) (prefix rule).
#   TC-C3 variances are two-pass (centre, then sum squares), for accuracy.

census_targets <- function(pop, m_h, n0 = 20) {
  psu <- pop$psu; units <- pop$units
  # Canonical order is what lets PSU totals be formed by grouping; assert it.
  if (any(order(psu$h, psu$k) != seq_len(nrow(psu))))
    stop("psu must be in canonical (h, k) order", call. = FALSE)
  idx <- rep(seq_len(nrow(psu)), psu$M)
  if (nrow(units) != length(idx) || any(units$h != psu$h[idx]) || any(units$k != psu$k[idx]))
    stop("units must be in canonical (h, k, i) order and match psu$M", call. = FALSE)

  strata <- sort(unique(psu$h))
  H <- length(strata)
  M_pop_h <- as.vector(table(factor(psu$h, levels = strata)))
  if (length(m_h) == 1L) m_h <- rep(m_h, H)
  if (length(m_h) != H || any(!is.finite(m_h)) || any(m_h != floor(m_h)))
    stop("m_h must be a whole-number scalar or one value per stratum (TC-C1)", call. = FALSE)
  # m_h >= 2: sigma2_m is E_des[se^2] with the m_h/(m_h-1) factor, undefined at 1.
  if (any(m_h < 2) || any(m_h > M_pop_h)) stop("need 2 <= m_h <= M_pop,h in every stratum", call. = FALSE)
  # pi_{i|k} = n0/M_k must lie in (0, 1].
  if (!is.numeric(n0) || length(n0) != 1L || n0 != floor(n0) || n0 < 1 || n0 > min(psu$M)) stop("need 1 <= n0 <= min M_k (pi_{i|k} <= 1)", call. = FALSE)

  N <- nrow(units)
  Psi_N <- mean(units$phi0)
  Psi_N_gf <- mean(units$Q1 - units$Q0)                   # the paper

  z <- units$phi0 - Psi_N
  t_k <- as.vector(rowsum(z, idx, reorder = TRUE))
  zbar <- t_k / psu$M
  S2_z <- as.vector(rowsum((z - zbar[idx])^2, idx, reorder = TRUE)) / (psu$M - 1)
  sh <- match(psu$h, strata)
  t_bar <- as.vector(rowsum(t_k, sh, reorder = TRUE)) / M_pop_h
  S2_t <- as.vector(rowsum((t_k - t_bar[sh])^2, sh, reorder = TRUE)) / (M_pop_h - 1)
  v2 <- as.vector(rowsum(psu$M * (psu$M - n0) * S2_z / n0, sh, reorder = TRUE))

  B_h <- M_pop_h^2 * S2_t / (m_h * N^2)
  V_h <- M_pop_h * v2 / (m_h * N^2)
  f_h <- m_h / M_pop_h

  list(Psi_N = Psi_N, Psi_N_gf = Psi_N_gf,
       per_stratum = data.frame(h = strata, M_pop_h = M_pop_h, m_h = m_h, f_h = f_h,
                                B_h = B_h, V_h = V_h, beta_h = B_h / (B_h + V_h)),
       sigma2_m = sum(B_h + V_h),
       sigma2_m_wor = sum((1 - f_h) * B_h + V_h),
       Delta_m = sum(f_h * B_h),
       N = N, M_pop = nrow(psu), m = sum(m_h), n0 = n0,
       label = paste("census targets;", pop$meta$label))
}

# E2 (FPC sweep) census summaries: pure algebra on census_targets() output,
# nothing estimated. Per stratum: the stage-1 FPC 1 - f_h, sigma2_m,h = B_h + V_h,
# sigma2_m,wor,h = (1 - f_h) B_h + V_h, Delta_h = f_h B_h and f_h V_h (U8).
# Totals: the G1 quantities (Delta_m / sigma2_m, sigma_m / sigma_m,wor), the pooled
# beta = sum B_h / sum (B_h + V_h), U8's negative-limb census deficit
# sum_h f_h V_h (a blanket single-stage fpc targets sigma2_m,wor - sum_h f_h V_h), the
# identity gap |(sigma2_m - sigma2_m,wor) - Delta_m|, and the C3b overlay value
# 2 Phi(z_{1-gamma/2} sigma_m / sigma_m,wor) - 1 at each level. All DESCRIPTIVE:
# no grade is attached here; C3a/C3b are graded on >= 100 populations.
census_fpc_summary <- function(cen, levels = c(0.80, 0.90, 0.95)) {
  ps <- cen$per_stratum
  ps$fpc_h <- 1 - ps$f_h
  ps$sigma2_m_h <- ps$B_h + ps$V_h
  ps$sigma2_m_wor_h <- (1 - ps$f_h) * ps$B_h + ps$V_h
  ps$Delta_h <- ps$f_h * ps$B_h
  ps$fV_h <- ps$f_h * ps$V_h
  ratio <- sqrt(cen$sigma2_m / cen$sigma2_m_wor)
  z <- stats::qnorm(1 - (1 - levels) / 2)
  list(per_stratum = ps, f_bar = mean(ps$f_h), fpc_bar = 1 - mean(ps$f_h),
       sigma_ratio = ratio, Delta_over_sigma2_m = cen$Delta_m / cen$sigma2_m,
       beta_pooled = sum(ps$B_h) / sum(ps$B_h + ps$V_h),
       identity_gap = abs((cen$sigma2_m - cen$sigma2_m_wor) - cen$Delta_m),
       single_stage_fpc_deficit = sum(ps$fV_h),
       overlay = stats::setNames(2 * stats::pnorm(z * ratio) - 1, sprintf("%.2f", levels)))
}

# ---- census algebra for the remaining E2 constructions ---------------

# E2 WR code-check target: the literal with-replacement design variance of the
# linearised total, sigma2_m,WR = sum_h [(1 - 1/M_pop,h) B_h + V_h] (census algebra). It differs
# from the WOR-defined sigma2_m = sum_h (B_h + V_h) by sum_h B_h / M_pop,h, which is why the
# check as first written would fail for correct code.
census_wr_target <- function(cen) {
  ps <- cen$per_stratum
  list(sigma2_m_wr = sum((1 - 1 / ps$M_pop_h) * ps$B_h + ps$V_h),
       wr_gap = sum(ps$B_h / ps$M_pop_h),                      # sigma2_m - sigma2_m,WR
       sigma2_m = cen$sigma2_m)
}

# E2 unequal-allocation invariance pair (master).
# The SAME population run twice: "high-on-high" gives the designated high-variance class K_B = 2
# m_h = m_high and K_B = 1 m_h = m_low; "high-on-low" swaps them. The class label is the frozen
# designation, never a ranking of realised B_h. Delta_m = N^-2 sum_h M_pop,h S2_t,h does not
# depend on the allocation, so the two Delta_m agree to floating point; the contrast lives in
# sigma2_m. Also returns G1's check that the intended variance contrast exists: the mean census
# S2_t,h (between-PSU variance of PSU totals) of the K_B = 2 strata exceeds that of K_B = 1.
census_unequal_pair <- function(pop, m_high = 6L, m_low = 2L, n0 = 20) {
  strata <- sort(unique(pop$psu$h))
  kb <- pop$psu$KB[match(strata, pop$psu$h)]
  if (anyNA(kb) || !all(kb %in% 1:2)) stop("every stratum needs a K_B class in {1, 2}", call. = FALSE)
  m_hh <- ifelse(kb == 2L, m_high, m_low)          # high-variance class gets the large allocation
  m_hl <- ifelse(kb == 2L, m_low, m_high)
  hh <- census_targets(pop, m_h = m_hh, n0 = n0)
  hl <- census_targets(pop, m_h = m_hl, n0 = n0)
  # S2_t,h recovered from B_h: B_h = M_pop,h^2 S2_t,h / (m_h N^2).
  s2t <- hh$per_stratum$B_h * hh$per_stratum$m_h * hh$N^2 / hh$per_stratum$M_pop_h^2
  list(high_on_high = hh, high_on_low = hl, m_h_hh = m_hh, m_h_hl = m_hl, KB = kb,
       delta_gap = abs(hh$Delta_m - hl$Delta_m),
       contrast = c(S2t_K2 = mean(s2t[kb == 2L]), S2t_K1 = mean(s2t[kb == 1L])),
       m_equal = sum(m_hh) == sum(m_hl))
}

# Psi-bridge statistic (secondary Psi; the vanishing-fraction cell is the ONLY
# bridge PASS cell): |Psi_N - Psi| / sigma_m,wor, Psi the generating-law parameter of P_{0,H}
# (the P0 grid). Graded at the terminal vanishing rung only, in production; reported here.
census_bridge <- function(cen, Psi) abs(cen$Psi_N - Psi) / sqrt(cen$sigma2_m_wor)

# ---- production-population gates ---------

# A population restricted to the strata `keep_h` (canonical order kept; used for per-class values).
pop_subset_strata <- function(pop, keep_h) {
  kp <- pop$psu$h %in% keep_h; ku <- pop$units$h %in% keep_h
  psu <- pop$psu[kp, , drop = FALSE]; units <- pop$units[ku, , drop = FALSE]
  rownames(psu) <- NULL; rownames(units) <- NULL
  meta <- pop$meta; meta$N <- nrow(units)
  list(psu = psu, units = units, meta = meta)
}

# The augmentation term a = H_g0 (Y - Q0(A, C)), H_g0 = A/g0 - (1 - A)/(1 - g0), at the marginal eta0.
census_aug <- function(units) {
  H <- units$A / units$g - (1 - units$A) / (1 - units$g)
  H * (units$Y - ifelse(units$A == 1L, units$Q1, units$Q0))
}

# Informativeness (a stopping gate at >= 0.30): corr(M_k, phibar_k | stratum), phibar_k the
# size-normalised census mean of phi0 over PSU k's units; both stratum-demeaned, pooled over strata.
census_informativeness <- function(pop) {
  psu <- pop$psu; idx <- rep(seq_len(nrow(psu)), psu$M)
  phibar <- as.vector(rowsum(pop$units$phi0, idx, reorder = TRUE)) / psu$M
  M <- as.numeric(psu$M); sh <- psu$h
  dm <- function(x) x - stats::ave(x, sh)
  x <- dm(M); y <- dm(phibar)
  den <- sqrt(sum(x^2) * sum(y^2))
  if (!is.finite(den) || den <= 0) NA_real_ else sum(x * y) / den
}

# DEFF_aug (definition, a REPORTED diagnostic): sigma2_m(a), the UC (no-FPC) design variance
# of the Hajek mean of a under stage 1 SRS-WOR of m_h PSUs per stratum (m_h equal across strata) and
# stage 2 SRS-WOR of n0 = 20 units (census_targets' B_h + V_h on a), over S2_a / n, the variance of the
# unweighted mean of n = H m_h n0 units drawn by SRS with replacement (S2_a the census unit variance,
# divisor N - 1). m_h cancels. Returns the pooled value and one per K_B class.
census_deff_aug <- function(pop, m_h = 2L, n0 = 20L) {
  one <- function(p) {
    p$units$phi0 <- census_aug(p$units)
    cen <- census_targets(p, m_h = m_h, n0 = n0)
    H <- length(unique(p$psu$h))
    cen$sigma2_m / (stats::var(p$units$phi0) / (H * m_h * n0))
  }
  kb_h <- tapply(pop$psu$KB, pop$psu$h, `[`, 1L)
  c(pooled = one(pop),
    K1 = one(pop_subset_strata(pop, as.integer(names(kb_h)[kb_h == 1L]))),
    K2 = one(pop_subset_strata(pop, as.integer(names(kb_h)[kb_h == 2L]))))
}

# Clip check (a stopping gate): every unit's marginal eta0 strictly inside the fitted clip.
census_clip <- function(units, q0 = 0.005, delta = 0.01) {
  r <- c(min_g = min(units$g), max_g = max(units$g),
         min_Q = min(units$Q1, units$Q0), max_Q = max(units$Q1, units$Q0))
  c(r, pass = as.numeric(r[["min_g"]] > delta && r[["max_g"]] < 1 - delta && r[["min_Q"]] > q0 && r[["max_Q"]] < 1 - q0))
}
