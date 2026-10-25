# The grading pass: the statistics, the keyed
# population-block bootstrap and the composition rules that turn saved PRODUCTION output into
# classifications. Post-processing only: nothing here draws a population, a sample or a fit.
# Requires R/rng_streams.R and R/grading.R (B_boot, the level-5 key, the percentile interval, the decision operator).
#
# Weighting. A joint (fitted) cell's statistic is computed on its B x R array with
# populations EQUALLY weighted. This is implemented as frequency weights: every row of population b carries
# f = L / R_b, with R_b that population's retained draws of the arm and L the least common multiple of the
# R_b, so every population carries total frequency L. The statistic is then the ordinary statistic of the
# frequency-expanded sample: with every replicate present f = 1 and it is exactly the plain statistic -- the
# type-7 quantile, the n - 1 sample SD, the empirical Kolmogorov distance.
# FAILED fits -- the (12 Oct 2026): any grade using an arm with a FAILED fit in
# that cell is UNGRADABLE (fitting failures); the statistic of the surviving replicates (computed as above, its
# population weights kept equal) is reported and labelled DESCRIPTIVE, never graded. Every failure is preserved;
# nothing is retried with new seeds.
#
# Bootstrap. Joint cells: population blocks only, every selected block keeps all of its rows unchanged;
# design-conditional cells (oracle, E4): the R draws at the fixed population; census-only cells: populations.
# One index floor(n U) + 1 per draw of a resample from the level-5 stream keyed (cell, statistic, grading pass);
# boot_stream() reproduces boot_index() exactly, in chunks, so 10,000 x 10,000 indices never sit in memory.

GRADE_PASS <- "v9.8-grade-1"           # the grading-pass label in the level-5 key
GRADE_LEVELS <- c(0.80, 0.90, 0.95)    # C3b's nominal levels (/C3b)

# margins, as compliant regions [lo, hi] for (one-sided margins have an infinite side).
GRADE_MARGINS <- list(
  z_location   = c(-0.10, 0.10),   # C1: |mean Z| <= 0.10, every gradeable C1 Z
  z_scale      = c(-0.05, 0.05),   # C1: |SD(Z) - 1| <= 0.05
  z_ks         = c(-Inf, 0.05),    # C1: Kolmogorov distance <= 0.05
  rvar_oracle_mean = c(-0.05, 0.05),  # C2 oracle limb: |mean R_var - 1| <= 0.05 (implementation check)
  rvar_oracle_sd   = c(-Inf, 0.15),   # C2 oracle limb: SD(R_var) within a population <= 0.15
  rvar_fit_mean    = c(-0.15, 0.15),  # C2 gradeable fitted: |mean R_var - 1| <= 0.15 at the terminal rung
  rvar_fit_sd      = c(-Inf, 0.15),   # C2 gradeable fitted: SD(R_var), joint, <= 0.15
  transfer         = c(-Inf, 0.25),   # C2 paired transfer: mean |se2_fit - se2_or| / sigma2_m <= 0.25
  coverage         = c(-0.02, Inf),   # C3b: cov - (1 - gamma) >= -0.02 at each level
  c3a              = c(-0.05, 0.05),  # C3a: |Var(S_des)/sigma2_m - (1 - Delta_m/sigma2_m)| <= 0.05 (signed)
  negligible       = c(-Inf, 0.10),   # C4a, C4b, C6 (both objects), sec. Psi: Q_0.9 <= 0.10
  wrong_clustering = c(0.85, 1.15))   # E3 C2 power controls (G6 limb (a) region)
GRADE_SPEC_REF <- list(z_location = "§12.4 :1286", z_scale = "§12.4 :1287", z_ks = "§12.4 :1288",
  rvar_oracle_mean = "§12.4 :1289", rvar_oracle_sd = "§12.4 :1290", transfer = "§12.4 :1291",
  rvar_fit_mean = "§12.4 :1292", rvar_fit_sd = "§12.4 :1293", coverage = "§12.4 :1294", c3a = "§12.4 :1295, R9",
  negligible = "§12.4 :1296, :1309-1311; A3 §A3.13", wrong_clustering = "A3 §A3.8, G6 limb (a)")

# ---- frequency-weighted statistics (exact plain statistics when every f = 1) --------------------------
gcd2 <- function(a, b) { while (b > 0) { t <- b; b <- a %% b; a <- t }; a }
lcm_int <- function(x) { x <- unique(as.integer(x)); l <- 1L; for (v in x) l <- l %/% gcd2(l, v) * v; as.integer(l) }

# Frequencies giving every population equal total weight: f = L / R_b.
pop_freq <- function(pop) {
  n_b <- as.vector(table(pop)[as.character(pop)])
  as.numeric(lcm_int(n_b) / n_b)
}

f_mean <- function(x, f) sum(f * x) / sum(f)
f_sd <- function(x, f) { n <- sum(f); m <- sum(f * x) / n; sqrt(sum(f * (x - m)^2) / (n - 1)) }
f_var <- function(x, f) f_sd(x, f)^2
# Type-7 quantile of the frequency-expanded sample, without expanding it.
f_q7 <- function(x, f, p) {
  o <- order(x); x <- x[o]; f <- f[o]
  n <- sum(f); h <- (n - 1) * p + 1
  cf <- cumsum(f)
  at <- function(k) x[which(cf >= k - 1e-9)[1]]              # the k-th order statistic of the expansion
  lo <- floor(h); hi <- ceiling(h)
  if (lo == hi) at(lo) else at(lo) + (h - lo) * (at(hi) - at(lo))
}
# Kolmogorov distance of the frequency-weighted empirical CDF from N(0, 1).
f_ks <- function(z, f) {
  o <- order(z); z <- z[o]; cw <- cumsum(f[o]) / sum(f)
  p <- stats::pnorm(z)
  max(abs(cw - p), abs(c(0, cw[-length(cw)]) - p))
}

# ---- the graded statistics: f(data rows, frequencies) -> one number --------------------------------------
# Each takes a data.frame d (the rows) and frequencies f; the column it reads is fixed by the statistic.
STAT <- list(
  mean_z      = function(col) function(d, f) f_mean(d[[col]], f),
  sd_z_minus1 = function(col) function(d, f) f_sd(d[[col]], f) - 1,
  ks_z        = function(col) function(d, f) f_ks(d[[col]], f),
  mean_minus1 = function(col) function(d, f) f_mean(d[[col]], f) - 1,
  sd          = function(col) function(d, f) f_sd(d[[col]], f),
  mean        = function(col) function(d, f) f_mean(d[[col]], f),
  mean_abs    = function(col) function(d, f) f_mean(abs(d[[col]]), f),
  q90_abs     = function(col) function(d, f) f_q7(abs(d[[col]]), f, 0.9),
  # C3b: coverage of oracle + deterministic sigma_m (interval (2), /C3b) minus the nominal level.
  cov_det_sigma_m = function(level) function(d, f) {
    z <- stats::qnorm(1 - (1 - level) / 2)
    f_mean(as.numeric(abs(d$psi_hat - d$Psi_N) <= z * d$sigma_m), f) - level
  },
  # C3a (signed): Var(S_des) / sigma2_m - (1 - Delta_m / sigma2_m); at a fixed population
  # 1 - Delta_m / sigma2_m = sigma2_m,wor / sigma2_m exactly (Delta_m = sigma2_m - sigma2_m,wor).
  c3a = function() function(d, f) f_var(d$psi_hat, f) / d$sigma_m[1]^2 - d$sigma_m_wor[1]^2 / d$sigma_m[1]^2
)

# ---- keyed bootstrap -------------------------------------------------------------------------------------
# Apply fun(index vector) to each of n_boot resamples of 1..n drawn from the level-5 stream; identical
# indices to boot_index(root, cell, statistic, grading_pass, n, n_boot) (tested), generated in chunks.
boot_stream <- function(root, cell, statistic, grading_pass, n, fun, n_boot = GRADE_B_BOOT, chunk = NULL) {
  if (!is.numeric(n) || n < 1 || n != floor(n)) stop("n must be a whole number >= 1", call. = FALSE)
  if (is.null(chunk)) chunk <- max(1L, min(n_boot, as.integer(floor(5e6 / n))))
  st <- rng_state(rng_key("postprocess", list(root = root, cell = cell, statistic = statistic, grading_pass = grading_pass)))
  out <- numeric(n_boot)
  rng_scoped({
    rng_set_state(st)
    done <- 0L
    while (done < n_boot) {
      k <- min(chunk, n_boot - done)
      idx <- matrix(as.integer(floor(n * stats::runif(n * k))) + 1L, nrow = k, ncol = n, byrow = TRUE)
      for (j in seq_len(k)) out[done + j] <- fun(idx[j, ])
      done <- done + k
    }
  })
  out
}

# Joint cell: population blocks; every selected block keeps all of its rows; frequencies
# f = L / R_b recomputed on the resample (a block drawn twice counts twice). d: rows; pop: their population.
boot_joint_stat <- function(d, pop, stat, key, n_boot = GRADE_B_BOOT) {
  pops <- sort(unique(pop))
  rows_by <- split(seq_len(nrow(d)), factor(pop, levels = pops))
  f_by <- lapply(rows_by, function(r) rep(1 / length(r), length(r)))
  L <- lcm_int(lengths(rows_by))
  est <- stat(d, pop_freq(pop))
  boot <- boot_stream(key$root, key$cell, key$statistic, key$pass, length(pops), function(i) {
    r <- unlist(rows_by[i], use.names = FALSE)
    stat(d[r, , drop = FALSE], L * unlist(f_by[i], use.names = FALSE))
  }, n_boot = n_boot)
  list(estimate = est, boot = boot, B = length(pops), rows = nrow(d))
}

# Design-conditional cell: the R draws at one fixed population are resampled (f = 1).
boot_draws_stat <- function(d, stat, key, n_boot = GRADE_B_BOOT) {
  one <- rep(1, nrow(d))
  boot <- boot_stream(key$root, key$cell, key$statistic, key$pass, nrow(d),
                      function(i) stat(d[i, , drop = FALSE], one), n_boot = n_boot)
  list(estimate = stat(d, one), boot = boot, B = 1L, rows = nrow(d))
}

# Census-only cell (R = 0): one value per population; populations are resampled (type-7 quantiles etc.).
boot_census_stat <- function(values, stat, key, n_boot = GRADE_B_BOOT) {
  boot <- boot_stream(key$root, key$cell, key$statistic, key$pass, length(values),
                      function(i) stat(values[i]), n_boot = n_boot)
  list(estimate = stat(values), boot = boot, B = length(values), rows = length(values))
}

# on a bootstrap result: the 95% percentile interval against the compliant region.
grade_boot <- function(b, region) {
  ci <- percentile_interval(b$boot)
  data.frame(estimate = b$estimate, ci_lo = ci[1], ci_hi = ci[2], boot_sd = stats::sd(b$boot),
             region_lo = region[1], region_hi = region[2], B = b$B, rows = b$rows,
             disposition = decide(ci, region[1], region[2]), basis = "graded", stringsAsFactors = FALSE)
}

# ---- composition and the special dispositions ----------------------------------------------
# Joint (fitted) cells only: PASS iff every required grade-bearing metric PASSes, RED-FLAG iff any RED-FLAGs,
# UNGRADABLE otherwise. Design-conditional cells are never composed across populations.
# Output that cannot be graded dominates: NOT GRADED (incomplete or refused output), then BLOCKED (a failed or
# missing required check, or a non-finite value in a replicate marked OK), then the operator's outcomes, then the
# UNGRADABLE reasons in order: fitting failures, premise, MC precision.
compose <- function(d) {
  d <- as.character(d)
  if (!length(d)) return("UNGRADABLE (nothing graded)")
  for (pre in c("NOT GRADED", "BLOCKED")) { x <- d[startsWith(d, pre)]; if (length(x)) return(x[1]) }
  if (all(d == "PASS")) return("PASS")
  if (any(d == "RED-FLAG")) return("RED-FLAG")
  for (pat in c("^UNGRADABLE \\(fitting failures", "^UNGRADABLE \\(unresolved T2", "^UNGRADABLE \\(premise")) { x <- d[grepl(pat, d)]; if (length(x)) return(x[1]) }
  "UNGRADABLE (MC precision)"
}

# A power control is PREDICTED to fail its row (POWER CONTROLS): a RED-FLAG on the row is the expected
# outcome; a PASS means the diagnostic has no power (a RED-FLAG of the control).
control_disposition <- function(d) switch(d,
  "PASS" = "NO POWER (the control passed its row: RED-FLAG of the diagnostic)",
  "RED-FLAG" = "FAILS AS PREDICTED",
  d)

# Rows for an item that is not graded (premise, incomplete output, boundary, stress, display).
ungraded_row <- function(disposition, estimate = NA_real_, region = c(NA_real_, NA_real_), B = NA_integer_, rows = NA_integer_,
                         basis = if (is.na(estimate)) "" else "descriptive")
  data.frame(estimate = estimate, ci_lo = NA_real_, ci_hi = NA_real_, boot_sd = NA_real_,
             region_lo = region[1], region_hi = region[2], B = B, rows = rows,
             disposition = disposition, basis = basis, stringsAsFactors = FALSE)
