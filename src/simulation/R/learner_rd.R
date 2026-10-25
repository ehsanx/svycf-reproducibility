# R-d nuisance learner: the fixed whole-PSU training throttle on top of R-a
# (G3(a); level-4 keying).
#
# R-d is R-a (R/learner_ra.R) fitted on a random subsample of
# q_m(alpha) WHOLE PSUs drawn from the fold's TRAINING set. G5 (r_m tracks m^-alpha) is
# evaluated on the dry run's full ladder m in {200, 400, 800, 1600, 3200} (results/simulation/execution-check).
#
# Fixed settings: alpha in {0.20, 0.25, 0.30, 0.35}; m in {200, 400, 800,
# 1600, 3200}; q_m(alpha) = floor(q0 * (m/m0)^(2 alpha)) with (m0, q0) = (200, 30);
# q_m is the number of training PSUs made available to the throttled learner (at V = 2
# the training pool is m/2); the subsample is data-dependent/T_v-measurable (drawn from
# the fold's training set only); NO CLIPPING: q_m > pool is a hard error and must never
# fall back to R-a behaviour; mechanical validation: q_m must equal the
# fixed table exactly; the throttle is never retuned.
#
# TECHNICAL CONVENTIONS:
# C-RD-1 Arm names: "Rd-AIPW-a020", "-a025", "-a030", "-a035"; alpha
#        enters ONLY at level 4 through the arm, so the four rungs share the
#        replicate's sample and fold structure.
# C-RD-2 Subsample draw: the level-4 learner stream keyed (root, cell, pop, rep, arm, fold)
# gives one uniform per training PSU, the PSUs taken in canonical
#        (h, k) order; the q_m PSUs with the smallest uniforms are kept (SRS-WOR by
#        ranking uniforms, as the stage-1 sample draw TC-S1 does), ties broken by
#        canonical position. Scoped stream: the ambient RNG is never touched.
# C-RD-3 Whole PSUs: every training unit of a selected PSU enters the R-a fit, and no
#        other unit does. The fit is fit_ra() on exactly those rows, so R-d inherits
#        R-a's conventions (per-class fits, nAGQ = 9, GH64 marginal, clip/truncation).
# C-RD-4 Failures: q_m > pool is a hard specification error (a failed required check). An R-a fit
#        that fails after the single Nelder-Mead refit -- or a drawn subsample that holds no PSU of a
#        K_B class, so R-a cannot be fitted for it -- raises "ra_fit_failed": the caller records this
#        rung FAILED for the replicate (counted and reported per cell, never dropped silently),
#        never retried with other settings, never replaced by the full-training R-a fit.

RD_ALPHAS <- c(0.20, 0.25, 0.30, 0.35)                 # fixed
RD_M_GRID <- c(200L, 400L, 800L, 1600L, 3200L)          # fixed
RD_M0 <- 200L; RD_Q0 <- 30L                            # anchor (m0, q0), frozen
# The table, transcribed for mechanical validation (layer (i)): rows m, columns alpha.
RD_Q_TABLE <- matrix(c(30, 30, 30, 30,
                       39, 42, 45, 48,
                       52, 60, 68, 79,
                       68, 84, 104, 128,
                       90, 120, 158, 208),
                     nrow = 5L, byrow = TRUE,
                     dimnames = list(as.character(RD_M_GRID), sprintf("%.2f", RD_ALPHAS)))

rd_arm <- function(alpha) {
  if (length(alpha) != 1L || !any(abs(alpha - RD_ALPHAS) < 1e-12))
    stop("alpha must be one of the frozen 0.20, 0.25, 0.30, 0.35", call. = FALSE)
  sprintf("Rd-AIPW-a%03d", as.integer(round(100 * alpha)))    # C-RD-1
}

# q_m(alpha) = floor(q0 (m/m0)^(2 alpha)). Defined on the frozen m-grid only.
rd_q <- function(m, alpha) {
  rd_arm(alpha)
  if (length(m) != 1L || !(m %in% RD_M_GRID))
    stop("m must lie on R-d's frozen grid {200, 400, 800, 1600, 3200}", call. = FALSE)
  as.integer(floor(RD_Q0 * (m / RD_M0)^(2 * alpha)))
}

# The subsample (C-RD-2). train: the fold's training units (columns h, k at least).
# Returns the selected PSUs (canonical order), the pool size, q and the training rows used.
rd_subsample <- function(train, q, root, cell, population_id, rep_id, alpha, fold) {
  assert_engine_root(root, "rd_subsample")
  pool <- unique(data.frame(h = as.integer(train$h), k = as.integer(train$k)))
  pool <- pool[order(pool$h, pool$k, method = "radix"), , drop = FALSE]
  rownames(pool) <- NULL
  if (!is.numeric(q) || length(q) != 1L || q != floor(q) || q < 1)
    stop("q must be a whole number >= 1", call. = FALSE)
  if (q > nrow(pool))
    stop(sprintf("R-d throttle: q_m = %d exceeds the training pool of %d PSUs -- hard invariant, no clipping, no fallback to R-a",
                 as.integer(q), nrow(pool)), call. = FALSE)
  u <- rng_uniforms("learner", list(root = root, cell = cell, population_id = population_id,
                                    rep_id = rep_id, arm = rd_arm(alpha), fold = fold), nrow(pool))
  keep <- sort(order(u, seq_len(nrow(pool)))[seq_len(q)])
  psus <- pool[keep, , drop = FALSE]
  rownames(psus) <- NULL
  rows <- which(paste(train$h, train$k) %in% paste(psus$h, psus$k))
  list(psus = psus, pool = nrow(pool), q = as.integer(q), rows = rows)
}

#' R-d fit on one fold's training set: R-a on the throttled whole-PSU subsample.
fit_rd <- function(train, alpha, m, root, cell, population_id, rep_id, fold, a_cpsu_term) {
  train <- as.data.frame(train)
  q <- rd_q(m, alpha)
  sub <- rd_subsample(train, q, root, cell, population_id, rep_id, alpha, fold)
  missing_cls <- setdiff(sort(unique(train$KB)), sort(unique(train$KB[sub$rows])))
  if (length(missing_cls))
    stop(ra_fit_failed(sprintf("R-d throttle: the %d drawn PSUs contain no K_B = %s PSU, so R-a cannot be fitted for it (fit FAILED)",
                               q, paste(missing_cls, collapse = ", "))))
  fit <- fit_ra(train[sub$rows, , drop = FALSE], a_cpsu_term = a_cpsu_term)
  fit$throttle <- list(alpha = alpha, arm = rd_arm(alpha), m = as.integer(m), q = sub$q, pool = sub$pool,
                       psus = sub$psus, n_units_fit = length(sub$rows), n_units_pool = nrow(train),
                       fold = fold, convention = "C-RD-1..4 (R/learner_rd.R), level-4 stream of the adopted RNG contract (R2)")
  fit$conventions$learner <- "R-d: R-a on a q_m(alpha) whole-PSU subsample of the training fold (§1.7a)"
  fit
}

# G3(a) by arithmetic for one cell, every frozen alpha: the q_m value against the
# frozen table (layer (i), required), the hard invariant q_m <= m(V-1)/V (required), the
# design margin q_m <= 0.5 m(V-1)/V (reported: required = FALSE), and integrality V | m_h.
rd_g3_rows <- function(m, m_h, V) {
  pool <- m * (V - 1) / V
  do.call(rbind, lapply(RD_ALPHAS, function(a) {
    q <- rd_q(m, a); qt <- RD_Q_TABLE[as.character(m), sprintf("%.2f", a)]
    data.frame(check = c(sprintf("G3(a) %s: q_m = floor(30 (m/200)^(2 alpha)) equals the frozen §1.7a table (m = %d)", rd_arm(a), m),
                         sprintf("G3(a-i) %s: q_m <= training pool m(V-1)/V (hard invariant)", rd_arm(a)),
                         sprintf("G3(a-ii) %s: q_m <= 0.5 x pool (design margin; realized fraction in value)", rd_arm(a)),
                         sprintf("G3(c) %s: V divides m_h and q_m is whole", rd_arm(a))),
               value = c(q, q, q / pool, m_h %% V), reference = c(qt, pool, 0.5, 0),
               abs_diff = c(abs(q - qt), NA_real_, NA_real_, NA_real_), tolerance = c(0, NA_real_, NA_real_, NA_real_),
               pass = c(q == qt, q <= pool, q / pool <= 0.5, m_h %% V == 0 && q == floor(q)),
               required = c(TRUE, TRUE, FALSE, TRUE), stringsAsFactors = FALSE)
  }))
}

# Per-replicate throttle checks for one R-d arm's PSU-CF fit (cf from psu_cf_aipw with fit_rd):
# per fold v, the drawn PSU count equals q_m(alpha); every drawn PSU lies in fold v's training
# set; whole PSUs (exactly the drawn PSUs' training units were fitted); no drawn PSU is scored
# in fold v; fewer units than the full training pool (no fallback); the pool is m(V-1)/V PSUs.
check_rd_throttle <- function(sample, cf, alpha, m) {
  key <- paste(sample$h, sample$k)
  q <- rd_q(m, alpha); V <- length(cf$fits)
  do.call(rbind, lapply(seq_len(V), function(v) {
    th <- cf$fits[[v]]$throttle
    drawn <- paste(th$psus$h, th$psus$k)
    train_psu <- unique(key[cf$fit_rows[[v]]])
    score_psu <- unique(key[cf$score_rows[[v]]])
    fitted_units <- sum(key %in% drawn)          # every sampled unit of the drawn PSUs, anywhere
    lab <- sprintf("Rd throttle %s fold %d", rd_arm(alpha), v)
    data.frame(check = c(paste(lab, ": drawn whole PSUs = q_m(alpha)"),
                         paste(lab, ": every drawn PSU is in the fold's training set"),
                         paste(lab, ": whole PSUs -- fitted units = every sampled unit of the drawn PSUs"),
                         paste(lab, ": no drawn PSU is scored in this fold (disjoint)"),
                         paste(lab, ": fewer units than the full training pool (no fallback to R-a)"),
                         paste(lab, ": training pool = m(V-1)/V PSUs")),
               value = c(length(unique(drawn)), sum(drawn %in% train_psu), th$n_units_fit,
                         sum(drawn %in% score_psu), th$n_units_fit, length(train_psu)),
               reference = c(q, length(drawn), fitted_units, 0, length(cf$fit_rows[[v]]), m * (V - 1) / V),
               abs_diff = NA_real_, tolerance = 0,
               pass = c(length(unique(drawn)) == q && th$q == q,
                        all(drawn %in% train_psu),
                        th$n_units_fit == fitted_units && th$n_units_pool == length(cf$fit_rows[[v]]),
                        !any(drawn %in% score_psu),
                        th$n_units_fit < length(cf$fit_rows[[v]]) && q < length(train_psu),
                        length(train_psu) == m * (V - 1) / V && th$pool == length(train_psu)),
               stringsAsFactors = FALSE)
  }))
}
