# Grading primitives: the decision operator (an interval against a margin region) and the
# Monte Carlo interval.
#
#   joint (fitted) cells      population BLOCKS are the only resampling level: each selected population
#                             contributes ALL of its R draws, unchanged (never resampled within a block),
#                             and the complete predeclared statistic is recomputed on the resampled B x R
#                             array with populations equally weighted (each block's rows share weight
#                             1/(B R_b)); blocks are never replaced by population means unless the
#                             statistic itself is defined on population means.
#   design-conditional cells  the R draws at the fixed population (oracle cells, E4); one disposition per
#                             population, never composed across populations.
#   census-only cells (R = 0) the percentile bootstrap over populations only (the Psi-bridge's Q_0.9).
#   B_boot = 10,000; the fifth, post-processing RNG level keyed hash(cell, statistic, grading pass);
#   one index floor(n U) + 1 per draw of a resample; the 95% percentile interval (type-7 quantiles).
# Decision operator: PASS iff the whole interval lies in the compliant region; RED-FLAG iff it lies wholly outside;
#   UNGRADABLE (MC precision) otherwise.
# Requires R/rng_streams.R. Post-processing only: nothing here draws a population or a sample.

GRADE_B_BOOT <- 10000L
GRADE_LEVEL <- 0.95

# n_boot x n matrix of resampling indices from the level-5 stream (row b is resample b).
boot_index <- function(root, cell, statistic, grading_pass, n, n_boot = GRADE_B_BOOT) {
  if (!is.numeric(n) || n < 1 || n != floor(n)) stop("n must be a whole number >= 1", call. = FALSE)
  u <- rng_uniforms("postprocess", list(root = root, cell = cell, statistic = statistic, grading_pass = grading_pass),
                    n * n_boot)
  matrix(as.integer(floor(n * u)) + 1L, nrow = n_boot, ncol = n, byrow = TRUE)
}

# Weighted summaries (weights summing to 1).
w_mean <- function(x, w) sum(w * x)
w_sd <- function(x, w) { m <- sum(w * x); sqrt(sum(w * (x - m)^2)) }
w_quantile <- function(x, w, p) {
  o <- order(x); cw <- cumsum(w[o]) / sum(w)
  x[o][which(cw >= p - 1e-12)[1]]
}
w_ks_normal <- function(z, w) {
  o <- order(z); z <- z[o]; cw <- cumsum(w[o]) / sum(w)
  p <- stats::pnorm(z)
  max(abs(cw - p), abs(c(0, cw[-length(cw)]) - p))
}

# Population-equal weights for the rows of a joint array: every population's rows share 1/(B R_b).
pop_equal_weights <- function(pop) {
  n_b <- as.vector(table(pop)[as.character(pop)])
  w <- 1 / n_b
  w / sum(w)
}

# Joint cell. x: one value per row (or a matrix, one row per replicate); pop: the
# population of each row; stat(x_rows, w) the complete statistic. Returns the point estimate on the
# original array and the bootstrap replicates.
boot_joint <- function(x, pop, stat, idx) {
  x <- as.matrix(x)
  pops <- sort(unique(pop))
  if (ncol(idx) != length(pops)) stop("idx must have one column per population", call. = FALSE)
  rows_by <- split(seq_len(nrow(x)), factor(pop, levels = pops))
  est <- stat(x, pop_equal_weights(pop))
  reps <- vapply(seq_len(nrow(idx)), function(b) {
    blk <- rows_by[idx[b, ]]
    r <- unlist(blk, use.names = FALSE)
    w <- unlist(lapply(blk, function(v) rep(1 / length(v), length(v))), use.names = FALSE)
    stat(x[r, , drop = FALSE], w / sum(w))
  }, numeric(1))
  list(estimate = est, boot = reps)
}

# Census-only cell (R = 0): values, one per population; the percentile bootstrap over populations.
boot_census <- function(values, stat, idx) {
  if (ncol(idx) != length(values)) stop("idx must have one column per population", call. = FALSE)
  list(estimate = stat(values), boot = apply(idx, 1L, function(i) stat(values[i])))
}

# Design-conditional cell: values, one per draw at one fixed population; the draws are resampled.
boot_draws <- function(values, stat, idx) boot_census(values, stat, idx)

percentile_interval <- function(boot, level = GRADE_LEVEL)
  unname(stats::quantile(boot, c((1 - level) / 2, 1 - (1 - level) / 2), type = 7, names = FALSE))

# The decision operator against a compliant region [lo, hi] (either side may be infinite, for one-sided margins).
decide <- function(ci, lo = -Inf, hi = Inf) {
  if (length(ci) != 2L || anyNA(ci)) return("UNGRADABLE (no interval)")
  if (ci[1] >= lo && ci[2] <= hi) "PASS" else if (ci[2] < lo || ci[1] > hi) "RED-FLAG" else "UNGRADABLE (MC precision)"
}

grade <- function(b, lo = -Inf, hi = Inf, level = GRADE_LEVEL) {
  ci <- percentile_interval(b$boot, level)
  list(estimate = b$estimate, ci = ci, boot_sd = stats::sd(b$boot), disposition = decide(ci, lo, hi))
}
