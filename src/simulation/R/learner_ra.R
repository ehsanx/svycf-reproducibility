# R-a nuisance learner: the correctly specified, class-specific (K_B)-variance
# random-intercept model, Gauss-Hermite marginalised and prevalidated against the
# oracle eta0.
#
# The configuration below is fixed; a_cpsu = 1.25.
#
# What R-a is:
# - a class-specific (K_B)-variance random-effects model, Gauss-Hermite marginalised and
#   prevalidated against eta0; it contains the truth, so its L2 and quarter-rate premises hold by
#   construction and no empirical rate gate is needed.
# - The law it models: logistic A and Y with PSU random
#   intercepts u^A_k ~ N(0, sigma^2_{A,K_B}), u^Y_k ~ N(0, sigma^2_{Y,K_B}); K_B is
#   observed. The oracle plugs MARGINAL eta0, the GH integral over
#   u^A, u^Y at 64 nodes, so R-a's Q and g are marginal too.
# - The nuisance fit is unweighted: C_psu carries the size driver, so stage-2
#   selection is ignorable given C.
# - Fitted nuisances are clipped and the clipped values are what T2 uses;
#   constants q0 = 0.005, delta = 0.01, applied with
#   clip_fit() from R/remainders.R.
#
# Within-class model (the term A:Cpsu is present: a_cpsu = 1.25):
#   Q: Y ~ A + C1 + C2 + C3 + C4 + Cpsu + C1:C3 + A:Cpsu + (1 | psu)
#   g: A ~ C1 + C2 + C3 + C4 + Cpsu + C1:C3 + (1 | psu)
# Inside one K_B class this is exactly the generating model's form: 1{K_B=2} is a class
# intercept shift and 1{K_B=2}*Csum adds 0.25 to each C_j slope, both absorbed by
# fitting each class separately.
#
# THE CONFIGURATION:
# C-RA-1 Per-class fits: fixed effects AND the random-intercept variance fitted separately within
#        each K_B class. Each per-class model contains the truth and keeps one scalar random
#        effect (so adaptive GH is valid).
# C-RA-2 Likelihood: adaptive Gauss-Hermite, nAGQ = 9 (lme4), not Laplace. With PSU clusters of
#        fixed size n0 = 20 the Laplace estimator's bias does not shrink with m (measured:
#        Laplace sigma_hat 0.970 vs 0.979 under AGQ, stable from 7 to 25 nodes).
# C-RA-3 Fit acceptance: optimizer bobyqa (maxfun 1e5); lme4's convergence checks at their
#        defaults; every warning and message recorded in $diagnostics and the fit KEPT -- singular
#        fits (sigma_hat = 0, a boundary ML estimate) are kept. On a HARD fitting error the model is
#        refitted ONCE with the Nelder-Mead optimizer (RA_REFIT_OPTIMIZER); if that also errors, the
#        fit raises a condition of class "ra_fit_failed" and the caller marks the replicate's R-a
#        arm (its paired TMLE, and any R-d rung whose subsample fit failed) FAILED: counted and
#        reported per cell, never dropped silently, never replaced by another learner.
#        on_warning = "error" (tests only) turns a warning into an immediate error, without refit.
# C-RA-4 Grouping factor: psu = (h, k); k is unique within stratum h.
# C-RA-5 Marginalisation uses gh_marginal() and gauss_hermite(64) from
#        R/nuisance_truth.R and R/quadrature.R, the same code as eta0, with the
#        fitted class sigma in place of the true one.
# C-RA-6 Prevalidation (an implementation check, reported, not a gate): ra_synthetic_sample's
#        SYNTHETIC two-class sample, errors against eta0 on the p0_context() primary grid. The
#        tests' thresholds are test conventions, not a gate.
# C-RA-7 Rows are put in a canonical order (h, k, C1..C4, A, Y; radix sort) and the
#        PSU factor levels follow it, so permuting train's rows cannot change a
#        single bit of the fit (tested).
# Downstream: R-a AIPW is graded for C2 at m = 1600; R-a TMLE is the C6 backstop at m = 800;
# fixed-strata R-a feeds C4a; R-d is R-a on 30-208 whole PSUs and inherits this rule.

RA_GH_NODES <- 64L
RA_Q_CLIP   <- 0.005    # q0,
RA_G_TRUNC  <- 0.01     # delta, standing (v),
RA_NAGQ     <- 9L       # C-RA-2
RA_N0       <- 20L      # n0_h, (prevalidation sample only)

RA_OPTIMIZER <- c("bobyqa", "bobyqa")                  # C-RA-3
RA_REFIT_OPTIMIZER <- c("Nelder_Mead", "Nelder_Mead")  # the single refit after a hard error

ra_control <- function(optimizer = RA_OPTIMIZER) {
  lme4::glmerControl(optimizer = optimizer, optCtrl = list(maxfun = 1e5))
}

# The condition an R-a fit raises when the refit also fails. Callers catch it by class and
# record the arm FAILED (never a silent drop, never another learner).
ra_fit_failed <- function(msg) structure(class = c("ra_fit_failed", "error", "condition"),
                                          list(message = msg, call = NULL))

# Fixed-effect formulas. The S11 term enters only on request (no default).
ra_formulas <- function(a_cpsu_term) {
  if (missing(a_cpsu_term) || !is.logical(a_cpsu_term) || length(a_cpsu_term) != 1L || is.na(a_cpsu_term))
    stop("a_cpsu_term must be TRUE or FALSE (S11 A:Cpsu term; no default)", call. = FALSE)
  q <- if (a_cpsu_term) Y ~ A + C1 + C2 + C3 + C4 + Cpsu + C1:C3 + A:Cpsu
       else             Y ~ A + C1 + C2 + C3 + C4 + Cpsu + C1:C3
  list(Q = q, g = A ~ C1 + C2 + C3 + C4 + Cpsu + C1:C3)
}

.ra_check_train <- function(train) {
  need <- c("h", "k", "KB", "C1", "C2", "C3", "C4", "Cpsu", "A", "Y")
  miss <- setdiff(need, names(train))
  if (length(miss)) stop("train lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  if (anyNA(train[need])) stop("train has missing values", call. = FALSE)
  if (!all(train$A %in% 0:1) || !all(train$Y %in% 0:1)) stop("A and Y must be 0/1", call. = FALSE)
  if (!all(train$KB %in% 1:2)) stop("KB must be 1 or 2 (:367, :1221)", call. = FALSE)
  psu <- paste(train$h, train$k, sep = ":")
  if (any(tapply(train$KB, psu, function(x) length(unique(x))) > 1L))
    stop("KB must be constant within a PSU (it is a stratum label, :1221-1224)", call. = FALSE)
  if (any(tapply(train$Cpsu, psu, function(x) length(unique(x))) > 1L))
    stop("Cpsu must be constant within a PSU (:367, :1221)", call. = FALSE)
  invisible(psu)
}

# One glmer fit with every warning and message captured (C-RA-3). On a hard error it refits
# once with Nelder-Mead; a second error raises ra_fit_failed. A warning turned into an error
# (on_warning = "error", tests only) is re-raised at once, without refit.
.ra_glmer <- function(formula, data, nAGQ, on_warning, label) {
  attempt <- function(optimizer) {
    warns <- character(); msgs <- character()
    fit <- withCallingHandlers(
      lme4::glmer(formula, data = data, family = stats::binomial(), nAGQ = nAGQ, control = ra_control(optimizer)),
      warning = function(w) {
        if (identical(on_warning, "error"))
          stop(structure(class = c("ra_warning_as_error", "error", "condition"),
                         list(message = paste0("R-a ", label, ": ", conditionMessage(w)), call = NULL)))
        warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
      },
      message = function(m) { msgs <<- c(msgs, trimws(conditionMessage(m))); invokeRestart("muffleMessage") })
    list(fit = fit, warns = warns, msgs = msgs, optimizer = paste(optimizer, collapse = ","))
  }
  first_error <- ""
  res <- tryCatch(attempt(RA_OPTIMIZER), ra_warning_as_error = function(e) stop(e), error = function(e) e)
  if (inherits(res, "error")) {
    first_error <- conditionMessage(res)
    res <- tryCatch(attempt(RA_REFIT_OPTIMIZER), ra_warning_as_error = function(e) stop(e), error = function(e) e)
    if (inherits(res, "error"))
      stop(ra_fit_failed(sprintf("R-a %s FAILED after one Nelder-Mead refit -- first error: %s; refit error: %s",
                                 label, first_error, conditionMessage(res))))
  }
  fit <- res$fit; warns <- res$warns; msgs <- res$msgs
  conv <- fit@optinfo$conv
  lme4_msgs <- unlist(conv$lme4$messages, use.names = FALSE)
  sigma <- unname(sqrt(lme4::VarCorr(fit)$psu[1, 1]))
  # Fixed-effect SEs are kept only for prevalidation z-scores; prediction uses beta, sigma.
  se <- sqrt(diag(as.matrix(stats::vcov(fit))))
  names(se) <- names(lme4::fixef(fit))
  list(beta = lme4::fixef(fit), sigma = sigma, se = se,
       diag = data.frame(model = label, n_units = nrow(data), n_psu = nlevels(data$psu),
                         sigma = sigma, singular = lme4::isSingular(fit),
                         opt_code = if (is.null(conv$opt)) NA_integer_ else as.integer(conv$opt),
                         n_warnings = length(warns),
                         warnings = paste(warns, collapse = " | "),
                         lme4_messages = paste(lme4_msgs, collapse = " | "),
                         messages = paste(msgs, collapse = " | "),
                         optimizer = res$optimizer, refit = nzchar(first_error), first_error = first_error,
                         # lme4's gradient/Hessian checks are warnings, its
                         # singular-fit check a message (C-RA-3): a singular fit
                         # (sigma_hat = 0, a boundary ML estimate) is converged
                         # but not clean.
                         converged = length(warns) == 0L && identical(as.integer(conv$opt), 0L),
                         clean = length(warns) == 0L && length(lme4_msgs) == 0L &&
                                 identical(as.integer(conv$opt), 0L),
                         stringsAsFactors = FALSE))
}

# Fixed-effect linear predictor of one class's fit at newdata (A set to `a` for Q).
.ra_lp <- function(par, tt, newdata) {
  X <- stats::model.matrix(tt, newdata)
  if (!setequal(colnames(X), names(par$beta)))
    stop("design columns do not match the fitted fixed effects", call. = FALSE)
  drop(X[, names(par$beta), drop = FALSE] %*% par$beta)
}

#' Fit R-a on a fold's training set.
#'
#' train: one row per sampled unit, columns h, k (PSU id unique within stratum),
#'   KB, C1..C4, Cpsu, A, Y. Unweighted.
#' a_cpsu_term: REQUIRED logical; TRUE includes the S11 A:Cpsu term in Q. The
#'   frozen a_cpsu = 1.25, so the runners pass TRUE.
#' Returns list(Q = function(a, newdata), g = function(newdata), raw, par,
#'   diagnostics, conventions). Q and g are the MARGINAL, CLIPPED predictions
#'; raw holds the same marginals before clipping. newdata needs C1..C4,
#'   Cpsu, KB; a is 0/1, scalar or one per row.
fit_ra <- function(train, a_cpsu_term, nAGQ = RA_NAGQ, on_warning = c("record", "error"),
                   q0 = RA_Q_CLIP, delta = RA_G_TRUNC, gh_nodes = RA_GH_NODES) {
  if (missing(a_cpsu_term)) stop("a_cpsu_term is required (S11; no default)", call. = FALSE)
  on_warning <- match.arg(on_warning)
  fm <- ra_formulas(a_cpsu_term)
  train <- as.data.frame(train)
  .ra_check_train(train)
  # Canonical row order (C-RA-7): the fit becomes a function of the training SET,
  # not of the order rows arrive in, which U6's invariance limb relies on for R-a
  # ("a deterministic function of the training set").
  train <- train[order(train$h, train$k, train$C1, train$C2, train$C3, train$C4, train$A, train$Y,
                       method = "radix"), , drop = FALSE]
  psu <- paste(train$h, train$k, sep = ":")
  classes <- sort(unique(as.integer(train$KB)))
  par <- list(); diag <- list()
  for (cls in classes) {
    d <- train[train$KB == cls, , drop = FALSE]
    d$psu <- factor(psu[train$KB == cls], levels = unique(psu[train$KB == cls]))
    fq <- .ra_glmer(stats::update(fm$Q, . ~ . + (1 | psu)), d, nAGQ, on_warning, sprintf("Q[K_B=%d]", cls))
    fg <- .ra_glmer(stats::update(fm$g, . ~ . + (1 | psu)), d, nAGQ, on_warning, sprintf("g[K_B=%d]", cls))
    par[[as.character(cls)]] <- list(Q = fq[c("beta", "sigma", "se")], g = fg[c("beta", "sigma", "se")])
    diag[[length(diag) + 1L]] <- cbind(KB = cls, fq$diag)
    diag[[length(diag) + 1L]] <- cbind(KB = cls, fg$diag)
  }
  diag <- do.call(rbind, diag)
  tt <- list(Q = stats::delete.response(stats::terms(fm$Q)), g = stats::delete.response(stats::terms(fm$g)))
  gh <- gauss_hermite(gh_nodes)

  marginal <- function(which, newdata, a = NULL) {
    if (!all(newdata$KB %in% classes))
      stop("newdata has a K_B class R-a was not trained on", call. = FALSE)
    if (!is.null(a)) {
      if (!all(a %in% 0:1) || !(length(a) %in% c(1L, nrow(newdata)))) stop("a must be 0/1, scalar or one per row", call. = FALSE)
      newdata$A <- a
    }
    out <- rep(NA_real_, nrow(newdata))
    for (cls in classes) {
      i <- which(newdata$KB == cls)
      if (!length(i)) next
      p <- par[[as.character(cls)]][[which]]
      lp <- .ra_lp(p, tt[[which]], newdata[i, , drop = FALSE])
      out[i] <- gh_marginal(lp, rep(p$sigma, length(i)), gh)
    }
    out
  }
  raw <- list(Q = function(a, newdata) marginal("Q", newdata, a), g = function(newdata) marginal("g", newdata))
  clipped <- clip_fit(raw, q0 = q0, delta = delta)
  c(clipped, list(raw = raw, par = par, diagnostics = diag,
                  conventions = list(status = "frozen (R4, spec v9.8)", a_cpsu_term = a_cpsu_term, nAGQ = nAGQ,
                                     optimizer = "bobyqa,bobyqa; maxfun=1e5; one Nelder-Mead refit on a hard error, then FAILED (R4)", gh_nodes = gh_nodes,
                                     q0 = q0, delta = delta, on_warning = on_warning)))
}

# ---- prevalidation against eta0 -------------------------------

# SYNTHETIC sample drawn directly from the law (C-RA-6). Per-PSU entity
# streams (R/rng_streams.R, the adopted RNG contract): within a PSU the draw order is
# [C_psu, u^A, u^Y uniforms], then 6 x n0 uniforms read column-major (C1..C4, the
# A and Y thresholds). A PSU's draws depend only on (root, seed, h, k), so a
# smaller sample is a prefix (k <= n) of a larger one with the same seed.
ra_synthetic_sample <- function(co, n_psu_per_class, seed, n0 = RA_N0, root = "SYNTHETIC-RA-PREVAL") {
  if (!is.character(root) || length(root) != 1L || !startsWith(root, "SYNTHETIC"))
    stop("root must start with 'SYNTHETIC' (no production populations)", call. = FALSE)
  stopifnot(n_psu_per_class >= 1, n0 >= 1)
  hk <- data.frame(h = rep(1:2, each = n_psu_per_class), k = rep(seq_len(n_psu_per_class), 2L))
  ctx <- list(root = root, master_id = "RA-PREVAL", population_id = seed)
  states <- rng_states(rng_entity_keys("population", ctx, "psu", hk))
  draws <- rng_scoped(lapply(seq_len(nrow(hk)), function(r) {
    rng_set_state(states[r, ])
    list(u0 = stats::runif(3L), U = matrix(stats::runif(6L * n0), nrow = n0))
  }))
  KB <- ifelse(hk$h %% 2L == 1L, 1L, 2L)                       # parity rule,
  u0 <- do.call(rbind, lapply(draws, `[[`, "u0"))
  U <- do.call(rbind, lapply(draws, `[[`, "U"))
  idx <- rep(seq_len(nrow(hk)), each = n0)
  d <- data.frame(h = hk$h[idx], k = hk$k[idx], i = rep(seq_len(n0), nrow(hk)), KB = KB[idx],
                  C1 = -1 + 2 * U[, 1], C2 = -1 + 2 * U[, 2], C3 = -1 + 2 * U[, 3], C4 = -1 + 2 * U[, 4],
                  Cpsu = (-1 + 2 * u0[, 1])[idx])
  uA <- (stats::qnorm(u0[, 2]) * co$sigma_A[KB])[idx]          # u^A_k ~ N(0, sigma^2_{A,K_B}),
  uY <- (stats::qnorm(u0[, 3]) * co$sigma_Y[KB])[idx]          # u^Y_k ~ N(0, sigma^2_{Y,K_B}),
  d$A <- as.integer(U[, 5] < stats::plogis(lp_g(d, d$KB, co) + uA))
  d$Y <- as.integer(U[, 6] < stats::plogis(lp_Q(d$A, d, d$KB, co) + uY))
  attr(d, "label") <- paste0("DEV/provisional SYNTHETIC sample, root ", root, ", seed ", seed)
  d
}

# The true within-class coefficients of in the fitted parameterisation.
ra_true_coef <- function(co, a_cpsu_term) {
  out <- list()
  for (cls in 1:2) {
    k2 <- as.numeric(cls == 2L)
    g <- co$g; Q <- co$Q
    gb <- c("(Intercept)" = g[["icpt"]] + g[["kb2"]] * k2, C1 = g[["c1"]], C2 = g[["c2"]], C3 = g[["c3"]],
            C4 = g[["c4"]], Cpsu = g[["cpsu"]], "C1:C3" = g[["c1c3"]])
    s <- Q[["kb2_csum"]] * k2
    qb <- c("(Intercept)" = Q[["icpt"]] + Q[["kb2"]] * k2, A = Q[["a"]], C1 = Q[["c1"]] + s, C2 = Q[["c2"]] + s,
            C3 = Q[["c3"]] + s, C4 = Q[["c4"]] + s, Cpsu = Q[["cpsu"]], "C1:C3" = Q[["c1c3"]])
    if (a_cpsu_term) qb <- c(qb, "A:Cpsu" = co$a_cpsu)
    out[[as.character(cls)]] <- list(Q = list(beta = qb, sigma = co$sigma_Y[cls]),
                                     g = list(beta = gb, sigma = co$sigma_A[cls]))
  }
  out
}

# The p0_context() primary grid and its eta0, without building the (large)
# second grid that only U9 needs. Same rule and orders as p0_context's default.
ra_p0_primary <- function(H, co) {
  prim <- eval(formals(p0_context)$primary)
  gp <- do.call(p0_grid, c(list(H = H), prim))
  list(primary = gp, primary_truth = truth_at(gp, co, RA_GH_NODES), primary_hash = p0_grid_hash(gp))
}

# Error of a fit against eta0 on the primary grid: max, P0-mean and L2(P0) of
# |fit - eta0| for Q(1,.), Q(0,.) and g, overall and per class.
ra_errors_vs_truth <- function(fit, ctx) {
  gp <- ctx$primary; tr <- ctx$primary_truth
  est <- list(Q1 = fit$Q(1, gp), Q0 = fit$Q(0, gp), g = fit$g(gp))
  rows <- list()
  for (nm in names(est)) {
    e <- abs(est[[nm]] - tr[[nm]])
    for (cls in c(NA, 1L, 2L)) {
      i <- if (is.na(cls)) seq_along(e) else which(gp$KB == cls)
      w <- gp$w[i] / sum(gp$w[i])
      rows[[length(rows) + 1L]] <- data.frame(
        nuisance = nm, KB = if (is.na(cls)) "all" else as.character(cls),
        max_abs = max(e[i]), mean_abs_P0 = sum(w * e[i]), L2_P0 = sqrt(sum(w * e[i]^2)),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

#' Prevalidate R-a against eta0: fit on a SYNTHETIC sample
#' with n_psu_per_class PSUs per K_B class and report |Qhat - Q0|, |ghat - g0| on
#' the primary P0 grid, plus fitted vs true sigma and fixed effects.
#' ctx: a p0_context() (or ra_p0_primary()) for the same co; built if NULL.
#' The A:Cpsu term follows co$a_cpsu: the correct model has it iff a_cpsu != 0.
prevalidate_ra <- function(co, n_psu_per_class, seed, ctx = NULL, H = 2L, n0 = RA_N0, nAGQ = RA_NAGQ) {
  a_term <- co$a_cpsu != 0
  if (is.null(ctx)) ctx <- ra_p0_primary(H, co)
  t0 <- proc.time()[["elapsed"]]
  train <- ra_synthetic_sample(co, n_psu_per_class, seed, n0 = n0)
  fit <- fit_ra(train, a_cpsu_term = a_term, nAGQ = nAGQ)
  fit_sec <- proc.time()[["elapsed"]] - t0
  truth <- ra_true_coef(co, a_term)
  sig <- do.call(rbind, lapply(names(fit$par), function(cls) data.frame(
    KB = cls, nuisance = c("Q", "g"),
    sigma_hat = c(fit$par[[cls]]$Q$sigma, fit$par[[cls]]$g$sigma),
    sigma_true = c(truth[[cls]]$Q$sigma, truth[[cls]]$g$sigma))))
  coef <- do.call(rbind, lapply(names(fit$par), function(cls) do.call(rbind, lapply(c("Q", "g"), function(w) {
    b <- fit$par[[cls]][[w]]$beta; bt <- truth[[cls]][[w]]$beta[names(b)]
    se <- fit$par[[cls]][[w]]$se[names(b)]
    data.frame(KB = cls, nuisance = w, term = names(b), hat = unname(b), true = unname(bt),
               se = unname(se), z = unname((b - bt) / se))
  }))))
  list(label = "R-a prevalidation (SYNTHETIC sample; an implementation check, reported, not a gate: R4)",
       n_psu_per_class = n_psu_per_class, n0 = n0, seed = seed, a_cpsu = co$a_cpsu, nAGQ = nAGQ,
       grid_hash = ctx$primary_hash, grid_points = nrow(ctx$primary), fit_seconds = fit_sec,
       errors = ra_errors_vs_truth(fit, ctx), sigma = sig, coef = coef, diagnostics = fit$diagnostics)
}
