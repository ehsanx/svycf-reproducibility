# T2 BY DETERMINISTIC REFINEMENT: the ladder, the convergence tolerance, the trigger, the per-fit cap and the
# unresolved rule. U9's absolute 1e-6 is KEPT.
#
# Why. Where a fold fit's clipped or truncated nuisance puts a kink in the T2
# integrand, the fixed primary grid GL(6,2) can miss U9's 1e-6 (1.5e-5 on one dry-run fit). The closed form
# itself is right: it equals the brute force on identical nodes to 1.7e-18.
#
# What changes, and what does not.
#   - A fold fit that does NOT trigger keeps the frozen primary grid's P0 f_v, bit for bit: the same nodes, truth and
#     arithmetic as p0_f_closed() (tested with identical()). An arm-replicate none of whose folds trigger has the T2
#     t2_closed_form() gives, bit for bit.
#   - A fold fit that triggers climbs a fixed ladder of finer Gauss-Legendre grids (the same law, the same rule
#     family as the primary, so rung 0 IS the frozen primary) until two successive rungs agree to tol_conv.
#   - Nothing else: no learner, clip constant, population, seed, frozen grid or tolerance changes.
#
# Settings (T2_REFINE_POLICY):
#   trigger    on the frozen primary grid's nodes, per fold fit and per surface (Q1, Q0, g): the smallest distance,
#              on the logit scale, from a fitted value to its clip bound (q0, 1 - q0 for Q; delta, 1 - delta for g) is
#              at most max(1, reach), where reach = half the largest logit difference between adjacent nodes,
#              summed over the five continuous axes (C1..C4, C_psu). A clipped node has distance 0, so every fit with
#              a clipped node triggers; reach extends that to a clip boundary that a locally linear surface could
#              reach between nodes.
#   ladder     GL(6,2) [frozen primary, rung 0] -> GL(8,3) -> GL(10,4) -> GL(12,6)
#   accept     the first rung k >= 1 with step |P0f_k - P0f_{k-1}| <= tol_conv = 5e-7 (half of U9's 1e-6) and, for
#              k >= 2, a step smaller than the one before (t2r_accept); P0f_k is used
#   cap        GL(12,6): at most 150,336 + 712,704 + 2,320,000 + 7,216,128 = 10,399,168 point evaluations per fold fit
#   unresolved no rung accepted: the fold's P0 f_v is UNRESOLVED. The arm-replicate's T2 (and T1 = paired - T2) is
#              recorded NA with T2_status "unresolved"; every rung's value and step stay in t2_refinement.csv; nothing
#              is substituted, the replicate is kept, and its estimation and coverage columns are untouched. The
#              grading pass marks the remainder grades that use it UNGRADABLE (R/grade_run.R).
# LIMITATION (documented and tested, tests/testthat/test-t2-refine.R): the trigger sees a fit only at the primary
# grid's nodes. A clipped region narrower than the node spacing, on a surface that is not locally linear, can lie
# between nodes without any node being near the bound; such a fit does not trigger and keeps the primary grid's
# value. The dry run's U9 limb 2 (an independent rule family, different nodes) and a validation on finer rules
# for a sample of untriggered fits are the checks for that case.
#
# U9 limb 2 for a refined fold (dry run only; production runs no U9): (a) closed form vs brute force on the primary
# grid's IDENTICAL nodes and truth, 1e-12 -- the code check, independent of quadrature; (b) the value used vs the
# brute force on the independent Clenshaw-Curtis family at 1e-6: CC(9,5) for an untriggered fold (the unrefined row,
# unchanged), CC(13,7) for a triggered one. A fold whose P0 f_v is unresolved has an UNRESOLVED row: it does not
# pass, and it is named as unresolved (not as a disagreement).
#
# Requires R/quadrature.R, R/nuisance_truth.R, R/p0_grid.R, R/remainders.R, R/quad_rules_io.R.

T2_REFINE_POLICY <- list(
  id = "t2-refine-v1 (PROPOSED, round 3d)",
  ladder = data.frame(rule = "gauss-legendre", n_c = c(6L, 8L, 10L, 12L), n_psu = c(2L, 3L, 4L, 6L), stringsAsFactors = FALSE),
  tol_conv = 5e-7,
  margin_floor = 1,
  bf_triggered = list(rule = "clenshaw-curtis", n_c = 13L, n_psu = 7L),
  u9_tol = 1e-6,
  identical_tol = 1e-12,
  chunk_points = 2e5
)

# A small stand-in for unit tests and --mode dev only (the runner refuses it in dry and prod).
T2_REFINE_POLICY_DEV <- within(T2_REFINE_POLICY, {
  id <- "t2-refine-v1 DEV (unit tests and --mode dev only)"
  ladder <- data.frame(rule = "gauss-legendre", n_c = c(4L, 5L, 6L), n_psu = c(2L, 2L, 3L), stringsAsFactors = FALSE)
  bf_triggered <- list(rule = "clenshaw-curtis", n_c = 6L, n_psu = 4L)
})

t2r_or <- function(a, b) if (is.null(a) || !length(a)) b else a

# Point evaluations of one rung (C1..C4 tensor x 29 C_psu cells x n_psu nodes x the K_B classes present).
t2r_points <- function(n_c, n_psu, classes = 2L) as.numeric(n_c)^4 * 29 * n_psu * classes
t2_policy_max_points <- function(policy = T2_REFINE_POLICY, classes = 2L)
  sum(t2r_points(policy$ladder$n_c, policy$ladder$n_psu, classes))

# ---- the saved base rules the refined grids need (grids/REFINE-MANIFEST.json) ----------------------------------
# The frozen artifact (grids/MANIFEST.json) holds only the frozen grids' rules and is not touched. A second, separately
# hashed artifact holds the finer rules; it is merged into the active store for keys the frozen one lacks. A key both
# hold must agree bit for bit, or loading stops.
t2_refine_rules <- function(engine_dir) {
  man <- jsonlite::fromJSON(file.path(engine_dir, "grids", "REFINE-MANIFEST.json"), simplifyVector = TRUE)
  c(quad_rules_load(file.path(engine_dir, "grids", man$file), man$sha256), list(manifest = man))
}
t2_refine_rules_extend <- function(loaded) {
  if (is.null(.QUAD_STORE$rules)) stop("activate the frozen quadrature artifact before extending it", call. = FALSE)
  have <- names(.QUAD_STORE$rules)
  for (k in intersect(names(loaded$rules), have))
    if (!identical(loaded$rules[[k]], .QUAD_STORE$rules[[k]])) stop("the refinement artifact disagrees with the frozen artifact on ", k, call. = FALSE)
  for (k in setdiff(names(loaded$rules), have)) .QUAD_STORE$rules[[k]] <- loaded$rules[[k]]
  .QUAD_STORE$sha256_refine <- loaded$sha256
  invisible(setdiff(names(loaded$rules), have))
}

# ---- rung 0: the frozen primary grid, from values already computed ----------------------------------------------
# The same expressions, in the same order, as p0_f_closed() (R/remainders.R), so the result is bit-identical.
t2r_closed_from_values <- function(Q1h, Q0h, gh, grid, truth) {
  if (any(gh <= 0 | gh >= 1)) stop("fitted propensity must lie strictly inside (0,1); clip it first")
  arm1 <- (Q1h - truth$Q1) * (gh - truth$g) / gh
  arm0 <- (Q0h - truth$Q0) * ((1 - gh) - (1 - truth$g)) / (1 - gh)
  sum(grid$w * (arm1 - arm0))
}
# The brute-force integrand of p0_f_bruteforce() (phi itself, enumerated over (A, Y)), at given values and truth.
t2r_brute_integrand <- function(Q1h, Q0h, gh, tr) {
  tot <- numeric(length(gh))
  for (a in 0:1) {
    pa <- if (a == 1) tr$g else 1 - tr$g
    Qa <- if (a == 1) tr$Q1 else tr$Q0
    for (y in 0:1) {
      py <- if (y == 1) Qa else 1 - Qa
      tot <- tot + pa * py * (aipw_phi(y, a, Q1h, Q0h, gh) - aipw_phi(y, a, tr$Q1, tr$Q0, tr$g))
    }
  }
  tot
}
t2r_closed_integrand <- function(Q1h, Q0h, gh, tr)
  (Q1h - tr$Q1) * (gh - tr$g) / gh - (Q0h - tr$Q0) * ((1 - gh) - (1 - tr$g)) / (1 - gh)

# ---- the trigger -------------------------------------------------------------------------------------------------
# The primary grid as a tensor: p0_grid() builds its rows with expand.grid(i1, i2, i3, i4, ip, KB) (i1 fastest) and
# drops a class whose share is 0, so values reshape to dim c(n_c, n_c, n_c, n_c, 29 n_psu, classes). NULL if the grid
# does not have that shape (then reach is not available and only the floor applies).
t2r_as_tensor <- function(v, grid) {
  def <- attr(grid, "definition")
  if (is.null(def)) return(NULL)
  n_ip <- 29L * as.integer(def$n_psu); n_c <- as.integer(def$n_c); ncls <- length(unique(grid$KB))
  d <- c(n_c, n_c, n_c, n_c, n_ip, ncls)
  if (length(v) != prod(d)) return(NULL)
  array(v, d)
}
# The largest absolute difference between nodes adjacent along axis k of a tensor.
t2r_adj_max <- function(a, k) {
  d <- dim(a)
  if (d[k] < 2L) return(0)
  s <- prod(d[seq_len(k - 1L)])
  lo <- which(slice.index(a, k) < d[k])
  max(abs(a[lo + s] - a[lo]))
}
t2r_trigger <- function(Q1h, Q0h, gh, grid, q0, delta, floor = 1) {
  surf <- list(Q1 = list(v = Q1h, lo = q0), Q0 = list(v = Q0h, lo = q0), g = list(v = gh, lo = delta))
  rows <- lapply(names(surf), function(nm) {
    v <- surf[[nm]]$v; lo <- surf[[nm]]$lo
    lv <- stats::qlogis(v)
    dist <- min(pmin(lv - stats::qlogis(lo), stats::qlogis(1 - lo) - lv))
    a <- t2r_as_tensor(lv, grid)
    reach <- if (is.null(a)) NA_real_ else sum(vapply(1:5, function(k) t2r_adj_max(a, k), 0)) / 2
    data.frame(surface = nm, dist = dist, reach = reach, margin = max(floor, reach, na.rm = TRUE),
               at_bound = sum(v == lo | v == 1 - lo), stringsAsFactors = FALSE)
  })
  r <- do.call(rbind, rows)
  r$triggered <- r$dist <= r$margin
  list(triggered = any(r$triggered), by_surface = r, dist_min = min(r$dist), reach_max = suppressWarnings(max(r$reach, na.rm = TRUE)),
       at_bound = sum(r$at_bound))
}

# ---- finer rungs: the tensor grid in chunks, its truth cached once per run --------------------------------------
# The law and node order of p0_grid(): the C1..C4 tensor (n_c^4 points, weights w/2 per axis) crossed with the 29
# size-weighted C_psu cells (n_psu nodes each) and the K_B classes with the H-specific shares. Never materialised as
# one data frame: a chunk is the C tensor crossed with a few C_psu nodes of one class.
t2r_tensor <- function(rule, n_c, n_psu, H, chunk_points = 2e5) {
  pc <- quadrature_rule(rule, n_c)
  cl <- size_law_cells(); Z <- sum(cl$M * (cl$hi - cl$lo))
  psu <- do.call(rbind, lapply(seq_len(nrow(cl)), function(i) {
    r <- quadrature_rule(rule, n_psu, cl$lo[i], cl$hi[i])
    data.frame(Cpsu = r$x, w_psu = r$w * cl$M[i] / Z)
  }))
  ix <- expand.grid(i1 = seq_len(n_c), i2 = seq_len(n_c), i3 = seq_len(n_c), i4 = seq_len(n_c), KEEP.OUT.ATTRS = FALSE)
  wc <- pc$w / 2
  C <- data.frame(C1 = pc$x[ix$i1], C2 = pc$x[ix$i2], C3 = pc$x[ix$i3], C4 = pc$x[ix$i4])
  shares <- p0_class_shares(H)
  nC <- nrow(C); k <- max(1L, as.integer(chunk_points %/% nC))
  chunks <- list()
  for (kb in which(shares > 0)) for (s in seq(1L, nrow(psu), by = k)) chunks[[length(chunks) + 1L]] <- list(kb = kb, j = s:min(nrow(psu), s + k - 1L))
  list(rule = rule, n_c = as.integer(n_c), n_psu = as.integer(n_psu), H = H, C = C, wC = wc[ix$i1] * wc[ix$i2] * wc[ix$i3] * wc[ix$i4],
       psu = psu, shares = shares, chunks = chunks,
       hash = digest::digest(list(rule, n_c, n_psu, H, as.matrix(C), psu, shares), algo = "sha256"))
}
t2r_chunk_data <- function(tp, ch) {
  nC <- nrow(tp$C); j <- ch$j
  nd <- data.frame(C1 = rep(tp$C$C1, length(j)), C2 = rep(tp$C$C2, length(j)), C3 = rep(tp$C$C3, length(j)),
                   C4 = rep(tp$C$C4, length(j)), Cpsu = rep(tp$psu$Cpsu[j], each = nC), KB = as.integer(ch$kb))
  w <- rep(tp$wC, length(j)) * rep(tp$psu$w_psu[j], each = nC) * tp$shares[ch$kb]
  list(nd = nd, w = w)
}
# One rung's tensor and its truth (at gh_n nodes), built on first use and kept in `cache` (an environment) for the
# rest of the run: every fold fit of the run shares them.
t2r_rung <- function(cache, ctx, rule, n_c, n_psu, gh_n, chunk_points) {
  key <- sprintf("%s:%d:%d:H%d:gh%d", rule, n_c, n_psu, ctx$H, gh_n)
  if (!is.null(cache[[key]])) return(cache[[key]])
  t0 <- proc.time()[["elapsed"]]
  tp <- t2r_tensor(rule, n_c, n_psu, ctx$H, chunk_points)
  truth_fun <- t2r_or(ctx$truth_fun, function(nd) truth_at(nd, ctx$co, gh_n))
  tp$truth <- lapply(tp$chunks, function(ch) truth_fun(t2r_chunk_data(tp, ch)$nd))
  tp$build_seconds <- proc.time()[["elapsed"]] - t0
  assign(key, tp, envir = cache)
  tp
}
# P0 f_v (closed form, or brute force) of one fit on one rung, summed over the chunks.
t2r_rung_value <- function(fit, tp, kind = c("closed", "brute")) {
  kind <- match.arg(kind)
  tot <- 0
  for (i in seq_along(tp$chunks)) {
    cd <- t2r_chunk_data(tp, tp$chunks[[i]]); tr <- tp$truth[[i]]
    Q1h <- fit$Q(1, cd$nd); Q0h <- fit$Q(0, cd$nd); gh <- fit$g(cd$nd)
    if (any(gh <= 0 | gh >= 1)) stop("fitted propensity must lie strictly inside (0,1); clip it first")
    tot <- tot + sum(cd$w * if (kind == "closed") t2r_closed_integrand(Q1h, Q0h, gh, tr) else t2r_brute_integrand(Q1h, Q0h, gh, tr))
  }
  tot
}

# Acceptance of the latest rung: its step is within tol_conv AND, from the second refined rung on, smaller than the
# step before it. On a kink a Gauss rule's error oscillates with the node count (tests/testthat/test-t2-refine.R: a
# one-dimensional clipped logistic gives errors +1.1e-6, +4.8e-6, +7e-7, -5e-7, +1.4e-6, -2.0e-6 on GL(4..12)), so two
# rungs can agree by coincidence; a non-decreasing step is read as not converged. A step is an error ESTIMATE, not a
# bound: the dry run's U9 limb 2 (an independent rule family) and the validation on finer rules are the checks on it.
t2r_accept <- function(steps, tol) {
  n <- length(steps)
  n >= 1L && steps[n] <= tol && (n == 1L || steps[n] < steps[n - 1L])
}

# ---- one fold fit ------------------------------------------------------------------------------------------------
# trigger: a trigger already computed for this fold (the paired TMLE's targeted fits use their AIPW fold's trigger:
# the tilt is smooth, so the kinks are the initial fit's). Returns the value used (NA if unresolved), the status, the
# accepted rung, every rung's value and step, the trigger and the refinement time.
t2_fold_refined <- function(fit, ctx, policy = T2_REFINE_POLICY, cache = new.env(parent = emptyenv()), trigger = NULL) {
  lad <- policy$ladder
  def <- attr(ctx$primary, "definition")
  if (!identical(c(def$rule, as.integer(def$n_c), as.integer(def$n_psu)), c(lad$rule[1], lad$n_c[1], lad$n_psu[1])))
    stop("the refinement ladder's rung 0 must be the context's primary grid", call. = FALSE)
  grid <- ctx$primary; tr <- ctx$primary_truth
  Q1h <- fit$Q(1, grid); Q0h <- fit$Q(0, grid); gh <- fit$g(grid)
  v0 <- t2r_closed_from_values(Q1h, Q0h, gh, grid, tr)
  if (is.null(trigger)) {                     # the fit's own clip constants; else the frozen q0 = 0.005, delta = 0.01
    q0 <- t2r_or(fit$conventions$q0, 0.005); dl <- t2r_or(fit$conventions$delta, 0.01)
    trigger <- t2r_trigger(Q1h, Q0h, gh, grid, q0, dl, policy$margin_floor)
  }
  out <- list(value = v0, status = "primary", rung = 0L, values = v0, steps = numeric(0), trigger = trigger, seconds = 0,
              brute_inputs = list(Q1h = Q1h, Q0h = Q0h, gh = gh))
  if (!isTRUE(trigger$triggered)) return(out)
  t0 <- proc.time()[["elapsed"]]
  vals <- v0; steps <- numeric(0); acc <- NA_integer_
  for (k in 2:nrow(lad)) {
    tp <- t2r_rung(cache, ctx, lad$rule[k], lad$n_c[k], lad$n_psu[k], ctx$gh, policy$chunk_points)
    vk <- t2r_rung_value(fit, tp, "closed")
    vals <- c(vals, vk); steps <- c(steps, abs(vk - vals[k - 1L]))
    if (t2r_accept(steps, policy$tol_conv)) { acc <- k - 1L; break }
  }
  out$values <- vals; out$steps <- steps; out$seconds <- proc.time()[["elapsed"]] - t0
  if (is.na(acc)) { out$value <- NA_real_; out$status <- "unresolved"; out$rung <- NA_integer_ }
  else { out$value <- vals[acc + 1L]; out$status <- "refined"; out$rung <- acc }
  out
}

# ---- one arm-replicate: T2 = sum_v omega_v P0 f_v -----------------------------------------------------------------
# All folds primary: T2 equals t2_closed_form()'s, bit for bit. Any fold unresolved: T2 is NA (status "unresolved"),
# and T2_last (the finest rung's sum) is kept only as a diagnostic.
t2_refined <- function(ctx, fits, omega, policy = T2_REFINE_POLICY, cache = new.env(parent = emptyenv()), triggers = NULL) {
  check_omega(omega, length(fits))
  fr <- lapply(seq_along(fits), function(v) t2_fold_refined(fits[[v]], ctx, policy, cache, trigger = if (is.null(triggers)) NULL else triggers[[v]]))
  st <- vapply(fr, `[[`, "", "status")
  per_fold <- vapply(fr, `[[`, 0, "value")
  last <- vapply(fr, function(x) x$values[length(x$values)], 0)
  status <- if (any(st == "unresolved")) "unresolved" else if (any(st == "refined")) "refined" else "primary"
  list(T2 = if (status == "unresolved") NA_real_ else sum(omega * per_fold), status = status, per_fold = per_fold, folds = fr,
       T2_last = sum(omega * last), n_refined = sum(st == "refined"), n_unresolved = sum(st == "unresolved"),
       n_triggered = sum(vapply(fr, function(x) isTRUE(x$trigger$triggered), TRUE)),
       rung_max = suppressWarnings(max(vapply(fr, function(x) if (is.na(x$rung)) length(x$values) - 1L else x$rung, 0L))),
       seconds = sum(vapply(fr, `[[`, 0, "seconds")), grid_hash = ctx$primary_hash)
}

# Results-row columns for an arm-replicate (NA T1/T2 when unresolved; the diagnostic sum kept apart).
t2_refined_cols <- function(t2r, psi_cf, psi_oracle, sig_wor) {
  T1 <- if (is.finite(t2r$T2)) t1_by_subtraction(psi_cf, psi_oracle, t2r$T2) else NA_real_
  data.frame(T2 = t2r$T2, T1 = T1, T2_over_sigma_wor = t2r$T2 / sig_wor, T1_over_sigma_wor = T1 / sig_wor,
             T2_status = t2r$status, T2_folds_triggered = t2r$n_triggered, T2_folds_refined = t2r$n_refined,
             T2_folds_unresolved = t2r$n_unresolved, T2_rung_max = t2r$rung_max,
             T2_unresolved_last = if (t2r$status == "unresolved") t2r$T2_last else NA_real_,
             T2_refine_seconds = t2r$seconds, stringsAsFactors = FALSE)
}

# The required "finite outputs" row of an arm with a T2: the unrefined label and values when T2 is resolved; when it is
# unresolved, T1 and T2 are NA by design (not a code failure), so the row checks the other outputs and says so.
# label: the unrefined label, which names T1 and T2 (", T1, T2" is removed for the unresolved row). Needs check_finite()
# (R/estimators.R).
t2_finite_check <- function(label, t2r, values, T1) {
  if (identical(t2r$status, "unresolved"))
    check_finite(paste0(sub(", T1, T2", "", label, fixed = TRUE), " (T1, T2 UNRESOLVED: recorded NA; rungs in t2_refinement.csv)"), values)
  else check_finite(label, c(values, list(T1, t2r$T2)))
}

# Per-fold diagnostics (t2_refinement.csv): every rung's value and step, the trigger, the status.
t2_refined_diag <- function(t2r, pop_id, rep_id, arm, max_rungs = 4L) {
  do.call(rbind, lapply(seq_along(t2r$folds), function(v) {
    x <- t2r$folds[[v]]; tg <- x$trigger
    vals <- c(x$values, rep(NA_real_, max_rungs - length(x$values)))[seq_len(max_rungs)]
    stp <- c(x$steps, rep(NA_real_, max_rungs - 1L - length(x$steps)))[seq_len(max_rungs - 1L)]
    cbind(data.frame(pop_id = pop_id, rep_id = rep_id, arm = arm, fold = v, triggered = isTRUE(tg$triggered),
                     dist_min = tg$dist_min, reach_max = tg$reach_max, nodes_at_bound = tg$at_bound,
                     status = x$status, rung = x$rung, value = x$value, seconds = x$seconds, stringsAsFactors = FALSE),
          stats::setNames(as.list(vals), paste0("rung", seq_len(max_rungs) - 1L)),
          stats::setNames(as.list(stp), paste0("step", seq_len(max_rungs - 1L))))
  }))
}

# ---- U9 limb 2 with refinement (dry run) --------------------------------------------------------------------------
# For each selected fold: (a) the identical-node code check (closed vs brute on the primary grid's nodes and truth,
# 1e-12); (b) the 1e-6 limb: an untriggered fold is the unrefined row exactly (closed on the primary vs brute force on
# the CC(9,5) second grid, GH128); a triggered fold compares the value used with the brute force on CC(13,7) (GH128);
# an unresolved fold's row is UNRESOLVED (does not pass).
u9_cross_check_refined <- function(ctx, fits, t2r, sel, policy = T2_REFINE_POLICY, cache = new.env(parent = emptyenv())) {
  rows <- list()
  for (v in sel) {
    fit <- fits[[v]]; x <- t2r$folds[[v]]; bi <- x$brute_inputs
    b_id <- sum(ctx$primary$w * t2r_brute_integrand(bi$Q1h, bi$Q0h, bi$gh, ctx$primary_truth))
    rows[[length(rows) + 1L]] <- data.frame(check = sprintf("U9 identical nodes (primary grid): closed form = brute force (fold %d, 1e-12)", v),
      value = x$values[1], reference = b_id, abs_diff = abs(x$values[1] - b_id), tolerance = policy$identical_tol,
      pass = check_pass(x$values[1], b_id, policy$identical_tol), stringsAsFactors = FALSE)
    if (!isTRUE(x$trigger$triggered)) {
      bf <- p0_f_bruteforce(fit, ctx)
      rows[[length(rows) + 1L]] <- data.frame(check = sprintf("U9 brute-force P0 f_v vs closed form (fold %d, 1e-6)", v),
        value = x$value, reference = bf, abs_diff = abs(x$value - bf), tolerance = policy$u9_tol, pass = check_pass(x$value, bf, policy$u9_tol),
        stringsAsFactors = FALSE)
    } else {
      b <- policy$bf_triggered
      tp <- t2r_rung(cache, ctx, b$rule, b$n_c, b$n_psu, ctx$gh_bf, policy$chunk_points)
      bf <- t2r_rung_value(fit, tp, "brute")
      lab <- sprintf("%s(%d,%d)", if (b$rule == "clenshaw-curtis") "CC" else "GL", b$n_c, b$n_psu)
      if (x$status == "unresolved") {
        rows[[length(rows) + 1L]] <- data.frame(check = sprintf("U9 brute-force P0 f_v vs closed form (fold %d, 1e-6; refined fold): UNRESOLVED -- the closed form did not converge by the ladder's cap; brute force on %s", v, lab),
          value = x$values[length(x$values)], reference = bf, abs_diff = abs(x$values[length(x$values)] - bf), tolerance = policy$u9_tol, pass = FALSE,
          stringsAsFactors = FALSE)
      } else {
        rows[[length(rows) + 1L]] <- data.frame(check = sprintf("U9 brute-force P0 f_v vs closed form (fold %d, 1e-6; refined fold, rung %d; brute force on %s)", v, x$rung, lab),
          value = x$value, reference = bf, abs_diff = abs(x$value - bf), tolerance = policy$u9_tol, pass = check_pass(x$value, bf, policy$u9_tol),
          stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

# The policy as recorded in a run's manifest.
t2_policy_record <- function(policy = T2_REFINE_POLICY, rules = NULL) {
  list(id = policy$id, status = "PROPOSED (author review pending; docs/ROUND3D.md)",
       ladder = sprintf("%s(%d,%d)", ifelse(policy$ladder$rule == "gauss-legendre", "GL", "CC"), policy$ladder$n_c, policy$ladder$n_psu),
       tol_conv = policy$tol_conv, margin_floor = policy$margin_floor,
       trigger = "per surface: min logit distance to the clip bound <= max(margin_floor, half the largest adjacent-node logit difference summed over C1..C4 and C_psu)",
       cap_points_per_fold_fit = t2_policy_max_points(policy),
       u9_triggered_brute_force = sprintf("%s(%d,%d), GH128", if (policy$bf_triggered$rule == "clenshaw-curtis") "CC" else "GL", policy$bf_triggered$n_c, policy$bf_triggered$n_psu),
       unresolved = "T2 and T1 recorded NA (T2_status unresolved); every rung kept in t2_refinement.csv; remainder grades UNGRADABLE",
       refine_rules_sha256 = if (is.null(rules)) NA_character_ else rules$sha256)
}
