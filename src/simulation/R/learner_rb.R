# R-b nuisance learner: a single fixed, moderate random forest.
#
# Every setting is fixed; the same conventions apply to R-c below. The level-4 seed stream is the
# random-number contract's (R/rng_streams.R).
# R-b's premise gates (nested L2 / quarter-rate verdicts on its own m-ladder)
# were evaluated on the dry run's ladder (results/simulation/execution-check). No forest T1/T2 is
# computed. Nothing here is graded.
#
# Settings: ranger, num.trees = 500, min.node.size = 25 (absolute,
# constant across m), mtry = 3 of the six covariates (C1..C4, C_psu, K_B),
# replace = TRUE, sample.fraction = 1, respect.unordered.factors = "order", no honesty.
# Fitted nuisances are clipped/truncated with the fixed constants:
# clip_fit() from R/remainders.R, q0 = 0.005, delta = 0.01.
#
# TECHNICAL CONVENTIONS:
# C-RB-1 Q by treatment arm: two forests, one on the A = 1 and one on the A = 0
#        training units, each on the six covariates. This is the reading of "mtry = 3
#        of the six covariates" under which A is not a split variable.
# C-RB-2 Regression forests on the 0/1 response (splitrule "variance");
#        predictions are leaf means in [0, 1].
# C-RB-3 K_B enters as a two-level factor, so respect.unordered.factors = "order" applies
#        as set; C1..C4 and C_psu are numeric.
# C-RB-4 Unweighted fits: like R-a, C_psu carries the size driver, so stage-2 selection is
#        ignorable given C.
# C-RB-5 RNG: the level-4 learner stream keyed
#        hash(cell, pop, rep, arm, fold_index) yields three uniforms in a
#        fixed order (Q1, Q0, g); seed = floor(u * (2^31 - 1)) + 1 goes to ranger's own
#        generator, with num.threads = 1. The arm key is "Rb-AIPW".
# C-RB-6 Pairing: Rb-TMLE targets Rb-AIPW's own clipped out-of-fold fits (shares
#        the stream and the fits), which keeps learner-seed noise out of the paired
#        contrast.
# C-RB-7 Failures: a warning is recorded in $diagnostics and the fit is kept; an error
#        stops this arm for the replicate and is recorded by the caller as a failed
#        required check (a learner failure), never retried with other settings.
# C-RB-8 Prediction seed. ranger 0.18.0's predict() draws
#        runif() from R's ambient generator whenever seed = NULL. Every prediction here
#        passes the forest's own training seed, so the ambient RNG is never touched.
#        Regression-forest predictions do not depend on that seed; rb_prediction_check()
#        verifies it on every fit the runner makes (identical() to seed = NULL).

# R-c: the E3 "aggressive" foil, a full-depth interpolating forest, with these settings:
# C-RC-1 ranger, 500 trees, mtry 3 of the six covariates, min.node.size = 1, replace = FALSE,
#        sample.fraction = 1: every tree sees the whole training set and grows to pure leaves,
#        so the forest interpolates its training data; no OOB error (there are no OOB units).
# C-RC-2 everything else as R-b's conventions: Q by treatment group, regression forests,
#        K_B factor, unweighted, explicit training seed at prediction (C-RB-1/2/3/4/8).
# C-RC-3 level-4 arms "Rc-AIPW-psucf" (folds 1..V), "Rc-AIPW-nocf" (fold 0, the single
#        full-sample fit) and "Rc-AIPW-indcf" (unit-level folds 1..V).
# R-c is ungraded (it may be UNGRADABLE by design) and, under, gets no T1/T2.

RB_NUM_TREES <- 500L
RB_MIN_NODE  <- 25L
RB_MTRY      <- 3L
RB_FEATURES  <- c("C1", "C2", "C3", "C4", "Cpsu", "KB")
RB_ARM_KEY   <- "Rb-AIPW"           # C-RB-5 / C-RB-6
RB_SETTINGS <- list(name = "R-b", num.trees = RB_NUM_TREES, mtry = RB_MTRY, min.node.size = RB_MIN_NODE,
                    replace = TRUE, sample.fraction = 1, oob.error = TRUE)            # fixed
RC_SETTINGS <- list(name = "R-c", num.trees = 500L, mtry = 3L, min.node.size = 1L,
                    replace = FALSE, sample.fraction = 1, oob.error = FALSE)          # C-RC-1
ARM_RC_PSUCF <- "Rc-AIPW-psucf"; ARM_RC_NOCF <- "Rc-AIPW-nocf"; ARM_RC_INDCF <- "Rc-AIPW-indcf"   # C-RC-3

.rb_x <- function(d) {
  x <- as.data.frame(d)[, RB_FEATURES]
  x$KB <- factor(x$KB, levels = 1:2)
  x
}

# Three ranger seeds for one (cell, pop, rep, arm, fold) from the level-4 stream (adopted contract).
# arm defaults to R-b's key; R-c passes its own arm (C-RC-3); fold = 0 is the no-CF fit.
rb_seeds <- function(root, cell, population_id, rep_id, fold, arm = RB_ARM_KEY) {
  assert_engine_root(root, "rb_seeds")
  u <- rng_uniforms("learner", list(root = root, cell = cell, population_id = population_id,
                                    rep_id = rep_id, arm = arm, fold = fold), 3L)
  as.integer(floor(u * 2147483647)) + 1L
}

.rb_forest <- function(y, x, seed, st = RB_SETTINGS) {
  warn <- character(0)
  t0 <- proc.time()[["elapsed"]]
  f <- withCallingHandlers(
    ranger::ranger(x = x, y = y, num.trees = st$num.trees, mtry = st$mtry, min.node.size = st$min.node.size,
                   replace = st$replace, sample.fraction = st$sample.fraction, respect.unordered.factors = "order",
                   splitrule = "variance", num.threads = 1L, seed = seed, verbose = FALSE, oob.error = st$oob.error),
    warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(forest = f, seed = seed, seconds = proc.time()[["elapsed"]] - t0, warnings = warn)
}

# C-RB-8: predictions always pass an explicit seed (the forest's training seed).
.rb_predict <- function(forest, x, seed) {
  if (length(seed) != 1L || is.na(seed)) stop("an explicit prediction seed is required (C-RB-8)", call. = FALSE)
  stats::predict(forest, data = x, seed = seed, num.threads = 1L)$predictions
}

#' R-b fit on one training set. Returns list(Q = function(a, newdata), g = function(newdata),
#' raw, diagnostics, conventions), Q and g CLIPPED as the arm uses them.
fit_rb <- function(train, seeds, q0 = 0.005, delta = 0.01, settings = RB_SETTINGS) {
  if (length(seeds) != 3L || anyNA(seeds)) stop("fit_rb: three seeds required (Q1, Q0, g)", call. = FALSE)
  train <- as.data.frame(train)
  x <- .rb_x(train)
  i1 <- which(train$A == 1); i0 <- which(train$A == 0)
  if (!length(i1) || !length(i0)) stop("fit_rb: a treatment arm has no training units", call. = FALSE)
  f1 <- .rb_forest(as.numeric(train$Y[i1]), x[i1, , drop = FALSE], seeds[1], settings)
  f0 <- .rb_forest(as.numeric(train$Y[i0]), x[i0, , drop = FALSE], seeds[2], settings)
  fg <- .rb_forest(as.numeric(train$A), x, seeds[3], settings)
  pred <- function(f, nd) .rb_predict(f$forest, .rb_x(nd), f$seed)
  raw <- list(Q = function(a, newdata) {
                if (!length(a) %in% c(1L, nrow(newdata)) || !all(a %in% 0:1)) stop("a must be 0/1, scalar or one per row", call. = FALSE)
                a <- rep_len(a, nrow(newdata))
                out <- numeric(nrow(newdata))
                if (any(a == 1)) out[a == 1] <- pred(f1, newdata[a == 1, , drop = FALSE])
                if (any(a == 0)) out[a == 0] <- pred(f0, newdata[a == 0, , drop = FALSE])
                out
              },
              g = function(newdata) pred(fg, newdata))
  comp <- list(Q1 = f1, Q0 = f0, g = fg)
  diag <- data.frame(component = names(comp), seed = seeds, n_train = c(length(i1), length(i0), nrow(train)),
                     oob_mse = vapply(comp, function(f) { pe <- f$forest$prediction.error
                                          if (is.null(pe) || !length(pe)) NA_real_ else as.numeric(pe) }, numeric(1)),
                     seconds = vapply(comp, `[[`, numeric(1), "seconds"),
                     n_warnings = vapply(comp, function(f) length(f$warnings), integer(1)),
                     warnings = vapply(comp, function(f) paste(f$warnings, collapse = " | "), ""),
                     stringsAsFactors = FALSE, row.names = NULL)
  c(clip_fit(raw, q0 = q0, delta = delta),
    list(raw = raw, diagnostics = diag, forests = lapply(comp, `[[`, "forest"),
         conventions = list(status = "frozen (spec v9.8; R18 for R-c)", learner = settings$name, num.trees = settings$num.trees,
                            mtry = settings$mtry, min.node.size = settings$min.node.size,
                            replace = settings$replace, sample.fraction = settings$sample.fraction,
                            respect.unordered.factors = "order", splitrule = "variance",
                            q_by_arm = TRUE, weights = "none", prediction_seed = "training seed (C-RB-8)",
                            q0 = q0, delta = delta,
                            ranger = as.character(utils::packageVersion("ranger")))))
}

# OI-RB-3 verification: the raw predictions fit_rb serves (explicit training seed, C-RB-8)
# against the same forests predicted the pre-fix way (seed = NULL). The seed = NULL calls
# run inside rng_scoped(), so the check itself cannot move the ambient RNG. Q1 and Q0 are
# compared on every row of newdata, whatever its A. One row per component.
rb_prediction_check <- function(fit, newdata) {
  x <- .rb_x(newdata)
  served <- list(Q1 = fit$raw$Q(1, newdata), Q0 = fit$raw$Q(0, newdata), g = fit$raw$g(newdata))
  legacy <- rng_scoped(lapply(fit$forests, function(f)
    stats::predict(f, data = x, seed = NULL, num.threads = 1L)$predictions))
  data.frame(component = names(served), n = nrow(x),
             identical = vapply(names(served), function(k) identical(served[[k]], legacy[[k]]), logical(1)),
             max_abs_diff = vapply(names(served), function(k) max(abs(served[[k]] - legacy[[k]])), numeric(1)),
             stringsAsFactors = FALSE, row.names = NULL)
}

# R-c fit (C-RC-1..3): R-b's code path with the full-depth settings.
fit_rc <- function(train, seeds, q0 = 0.005, delta = 0.01) fit_rb(train, seeds, q0, delta, settings = RC_SETTINGS)

ARM_RB_AIPW <- "Rb-AIPW"            # RNG draft vocabulary (TC-E5)
ARM_RB_TMLE <- "Rb-TMLE"

# R-c leakage-neutral sensitivity:
# the individual-CF R-c arm with the learner's C_psu feature replaced by its bin in
# {[-1, -0.5), [-0.5, 0), [0, 0.5), [0.5, 1]}. A SENSITIVITY ANALYSIS THAT ALSO REMOVES PREDICTIVE
# INFORMATION: C_psu is a true confounder of Q0 and g0, so coarsening it removes genuine
# predictive signal together with the PSU-identifying signal that individual-CF can leak; a change
# between the two arms cannot be attributed to leakage alone.
# C-RC-4 (confirmed) bins are represented by their midpoints (-0.75, -0.25, 0.25, 0.75); a tree only
#        sees the order of values, so any increasing coding gives the same forest. Coarsening is
#        applied to the learner's training AND prediction data only; the sample, its weights and
#        the oracle eta0 keep the true C_psu.
# C-RC-5 (confirmed) pairing: arm "Rc-AIPW-indcf-coarse" shares the unit-level folds and the level-4
#        stream of "Rc-AIPW-indcf" (the same three ranger seeds per fold), so the pair differs only in
#        the C_psu feature -- the C-RB-6 pairing logic applied to a sensitivity.
RC_CPSU_BREAKS <- c(-1, -0.5, 0, 0.5, 1)
RC_CPSU_MIDS <- c(-0.75, -0.25, 0.25, 0.75)
ARM_RC_INDCF_COARSE <- "Rc-AIPW-indcf-coarse"

coarsen_cpsu <- function(d) {
  d <- as.data.frame(d)
  d$Cpsu <- RC_CPSU_MIDS[findInterval(d$Cpsu, RC_CPSU_BREAKS, rightmost.closed = TRUE, all.inside = TRUE)]
  d
}

fit_rc_coarse <- function(train, seeds, q0 = 0.005, delta = 0.01) {
  inner <- fit_rc(coarsen_cpsu(train), seeds, q0, delta)
  wrap_Q <- function(f) function(a, newdata) f(a, coarsen_cpsu(newdata))
  wrap_g <- function(f) function(newdata) f(coarsen_cpsu(newdata))
  list(Q = wrap_Q(inner$Q), g = wrap_g(inner$g),
       raw = list(Q = wrap_Q(inner$raw$Q), g = wrap_g(inner$raw$g)),
       diagnostics = inner$diagnostics, forests = inner$forests, inner = inner,
       n_cpsu_values_train = length(unique(coarsen_cpsu(train)$Cpsu)),
       conventions = c(inner$conventions, list(coarsening = "C_psu -> quartile-bin midpoint (R15; C-RC-4)",
                                               pairing = "seeds of Rc-AIPW-indcf (C-RC-5)")))
}
