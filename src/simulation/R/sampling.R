# Two-stage SRS-WOR sampler and the U1 weight check.
#
# Every draw comes from the entity-addressed level-2 streams of R/rng_streams.R,
# under a SYNTHETIC root (tests) or the production root (dry and production runs).
#
# Requires, sourced beforehand: R/rng_streams.R (rng_key, rng_entity_keys,
# rng_states, rng_scoped, rng_set_state, RNG_DRAFT_VERSION).
#
# DESIGN:
#   frame   the cell's first H strata and, in each, its first J PSUs (prefix rule
#           on both axes); M_pop,h = J PSUs per stratum.
#   stage 1 SRS-WOR of m_h PSUs per stratum (SRS-WOR first stage).
#   stage 2 SRS-WOR of a fixed n0_h units per sampled PSU, pi_{i|k} = n0_h/M_k
#           (n0_h = 20 in every stratum, pi_{i|k} in [0.417, 1]).
#   weights w_i = 1/(pi_k pi_{i|k}) = (M_pop,h/m_h)(M_k/n0_h) (U1); never
#           calibrated or post-stratified to N, so N_hat = sum(w) stays an
#           estimator (U2).
#   RNG     only the level-2 sample stream hash(cell, population id, rep id), which
#           is arm-free and runs in a scope that restores the caller's RNG
#           (cell := (regime, H, m_h, J, m, f, E2-construction)), V and
#           alpha excluded. Any other key field (arm, fold, ...) is
#           ignored, so every arm of a replicate sees the identical sample.
#
# TECHNICAL CONVENTIONS:
#   TC-S1 SRS-WOR by ranking uniforms: entity stage1(h) draws
#         J uniforms and keeps the m_h smallest, ties broken by frame index k;
#         entity stage2(h,k) draws M_k uniforms and keeps the n0 smallest, ties by
#         unit index i. Identical to the prototype toy_sample() (tested).
#   TC-S2 weight arithmetic: w = (M_pop,h / m_h) * (M_k / n0), two binary64
#         divisions then one product, in that order; pi_psu = m_h/M_pop,h and
#         pi_unit = n0/M_k. 1/(pi_psu*pi_unit) is the same number algebraically but
#         not bit-for-bit (agrees to <= 16 eps relative, checked). U1's bit-exact
#         limb is therefore the integer-ratio form, recomputed from the census and
#         the recorded design on a separate code path (check_u1).
#   TC-S3 a sampled PSU with M_k < n0 would need pi_{i|k} > 1; the feasibility chain
# and the frozen law exclude it, so it is a hard
#         error. allow_take_all = TRUE (synthetic fixtures only) takes the whole PSU,
#         pi_unit = 1, recorded in the design attribute.
#   TC-S4 n0 is one number for all strata (the frozen n0_h is 20 everywhere).
#   TC-S5 m_h is a scalar, a length-H vector by stratum, or class-indexed
#         c(K1 = ., K2 = .) resolved through the stratum's K_B (unequal-allocation
#         pair); it must equal key$cell$m_h, so the sample stream's key
#         always names the allocation actually drawn.
#   TC-S6 output: rows in (h, k, i) order; columns h, k, i, w, pi_psu, pi_unit,
#         M_pop_h, M_k, then every unit column, then the PSU covariates psu_cols
#         (default K_B, C_psu: observed PSU-level confounders) that
#         the unit table lacks. The recorded design is attr(., "design").

SAMPLING_DEV_LABEL <- paste0(
  "sample: adopted RNG contract (sample stream, docs/RNG-CONVENTIONS-DRAFT.md §7, R2); ",
  "the runner's run mode says whether it is a dry-run or production draw")

.SAMPLING_ADDED_COLS <- c("w", "pi_psu", "pi_unit", "M_pop_h", "M_k")

.assert_sampling_root <- function(root) assert_engine_root(root, "sampler")

.whole_vec <- function(x, what, lo = 1) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x != floor(x)) || any(x < lo))
    stop(what, " must be whole number(s) >= ", lo, call. = FALSE)
  invisible(x)
}

# Frame of a cell: first H strata x first J PSUs, units sorted by
# (h, k, i) with a per-PSU row offset, so each draw indexes rows in O(sample).
# Built once per (population, H, J) and reusable across every rep.
sampling_frame <- function(pop, H, J, psu_cols = c("KB", "Cpsu")) {
  if (inherits(pop, "sampling_frame")) {
    if (pop$H != H || pop$J != J) stop("sampling_frame was built for another (H, J)", call. = FALSE)
    return(pop)
  }
  .whole_vec(H, "H"); .whole_vec(J, "J")
  if (!is.list(pop) || !is.data.frame(pop$psu) || !is.data.frame(pop$units))
    stop("pop must be a list with data frames $psu (h, k, M, ...) and $units (h, k, i, ...)",
         call. = FALSE)
  ps <- pop$psu
  un <- pop$units
  for (cl in c("h", "k", "M")) if (is.null(ps[[cl]])) stop("pop$psu lacks column ", cl, call. = FALSE)
  for (cl in c("h", "k", "i")) if (is.null(un[[cl]])) stop("pop$units lacks column ", cl, call. = FALSE)
  clash <- intersect(.SAMPLING_ADDED_COLS, names(un))
  if (length(clash)) stop("pop$units already has column(s) ", paste(clash, collapse = ", "), call. = FALSE)

  ps <- ps[ps$h <= H & ps$k <= J, , drop = FALSE]
  .whole_vec(ps$h, "pop$psu$h"); .whole_vec(ps$k, "pop$psu$k"); .whole_vec(ps$M, "pop$psu$M")
  ps <- ps[order(ps$h, ps$k), , drop = FALSE]
  # Prefix frame must be complete: exactly k = 1..J in every h = 1..H.
  if (nrow(ps) != H * J || any(ps$h != rep(seq_len(H), each = J)) || any(ps$k != rep(seq_len(J), H)))
    stop("pop$psu does not contain exactly PSUs k = 1..", J, " in each stratum h = 1..", H,
         " (prefix frame)", call. = FALSE)
  rownames(ps) <- NULL
  M <- as.numeric(ps$M)

  un <- un[un$h <= H & un$k <= J, , drop = FALSE]
  .whole_vec(un$h, "pop$units$h"); .whole_vec(un$k, "pop$units$k"); .whole_vec(un$i, "pop$units$i")
  p_un <- (un$h - 1) * J + un$k                      # frame PSU index, 1..H*J
  # Counts equal to M_k, every i in 1..M_k and no duplicate together mean each
  # PSU holds exactly i = 1..M_k, which the row offsets below rely on.
  if (any(un$i > M[p_un]) || any(tabulate(p_un, nbins = H * J) != M) ||
      anyDuplicated((p_un - 1) * max(M) + un$i))
    stop("pop$units does not hold exactly the units i = 1..M_k of every frame PSU", call. = FALSE)
  un <- un[order(p_un, un$i), , drop = FALSE]
  rownames(un) <- NULL

  M_pop_h <- tabulate(ps$h, nbins = H)               # = J by the completeness check
  un$M_pop_h <- M_pop_h[un$h]
  un$M_k <- M[(un$h - 1) * J + un$k]
  for (cl in setdiff(intersect(psu_cols, names(ps)), names(un)))
    un[[cl]] <- ps[[cl]][(un$h - 1) * J + un$k]

  structure(list(H = H, J = J, psu = ps, M = M, M_pop_h = M_pop_h, units = un,
                 offset = c(0, cumsum(M))[seq_len(H * J)], label = SAMPLING_DEV_LABEL),
            class = "sampling_frame")
}

# Resolve an allocation (TC-S5) to an integer vector of length H.
.resolve_m_h <- function(m_h, frame, what) {
  H <- frame$H
  if (length(m_h) == 1L && is.null(names(m_h))) return(rep(m_h, H))
  nm <- names(m_h)
  if (!is.null(nm) && all(grepl("^K[0-9]+$", nm))) {
    kb <- frame$psu$KB
    if (is.null(kb)) stop(what, " is class-indexed but pop$psu has no KB column", call. = FALSE)
    kb_h <- kb[!duplicated(frame$psu$h)]
    if (any(tapply(kb, frame$psu$h, function(x) length(unique(x))) != 1L))
      stop("K_B is not constant within a stratum", call. = FALSE)
    out <- unname(m_h[paste0("K", kb_h)])
    if (anyNA(out)) stop(what, " lacks a class present in the frame", call. = FALSE)
    return(out)
  }
  if (length(m_h) == H && is.null(nm)) return(m_h)
  stop(what, " must be a scalar, a length-H vector, or class-indexed c(K1 = ., K2 = .)",
       call. = FALSE)
}

# draw_sample(pop, m_h, n0, key, root): one two-stage SRS-WOR sample.
#   pop  list(psu = data.frame(h, k, M, KB, Cpsu, ...), units = data.frame(h, k, i, ..., phi0)),
#        or a sampling_frame() of it (faster across reps).
#   key  list(cell = <rng_streams cell list>, population_id, rep_id); other fields ignored.
draw_sample <- function(pop, m_h, n0 = 20, key, root = "SYNTHETIC-DEV-E2E",
                        allow_take_all = FALSE, psu_cols = c("KB", "Cpsu")) {
  .assert_sampling_root(root)
  if (!is.list(key) || is.null(key$cell)) stop("key must be list(cell, population_id, rep_id)", call. = FALSE)
  if (!is.null(key$root) && !identical(key$root, root))
    stop("key$root differs from root", call. = FALSE)
  ctx <- key
  ctx$root <- root
  sample_key <- rng_key("sample", ctx)               # validates cell, population_id, rep_id
  cell <- key$cell
  H <- cell$H
  J <- cell$J
  frame <- sampling_frame(pop, H, J, psu_cols)

  .whole_vec(n0, "n0")
  if (length(n0) != 1L) stop("n0 must be one number (frozen n0_h = 20 in every stratum, :1242)", call. = FALSE)
  mh <- .resolve_m_h(m_h, frame, "m_h")
  mh_cell <- .resolve_m_h(cell$m_h, frame, "key$cell$m_h")
  .whole_vec(mh, "m_h")
  if (!identical(as.numeric(mh), as.numeric(mh_cell)))
    stop("m_h differs from key$cell$m_h: the sample stream key must name the allocation drawn",
         call. = FALSE)
  if (any(mh > J)) stop("m_h exceeds J (SRS-WOR of m_h of J PSUs)", call. = FALSE)

  # Stage 1, entity stage1(h) (TC-S1).
  st1 <- rng_states(rng_entity_keys("sample", ctx, "stage1", list(h = seq_len(H))))
  sel_k <- rng_scoped(lapply(seq_len(H), function(h) {
    rng_set_state(st1[h, ])
    u <- stats::runif(J)
    sort(order(u, seq_len(J))[seq_len(mh[h])])
  }))
  s_h <- rep(seq_len(H), times = mh)
  s_k <- unlist(sel_k, use.names = FALSE)
  p_idx <- (s_h - 1) * J + s_k
  Mk <- frame$M[p_idx]
  if (any(Mk < n0) && !allow_take_all)
    stop("a sampled PSU has M_k < n0 = ", n0, ": pi_{i|k} > 1 is infeasible ",
         "(each sampled PSU needs at least n0 units); allow_take_all is for synthetic fixtures only",
         call. = FALSE)
  nk <- pmin(n0, Mk)

  # Stage 2, entity stage2(h,k) (TC-S1).
  st2 <- rng_states(rng_entity_keys("sample", ctx, "stage2", list(h = s_h, k = s_k)))
  sel_i <- rng_scoped(lapply(seq_along(p_idx), function(r) {
    rng_set_state(st2[r, ])
    u2 <- stats::runif(Mk[r])
    sort(order(u2, seq_len(Mk[r]))[seq_len(nk[r])])
  }))

  rows <- rep(frame$offset[p_idx], times = nk) + unlist(sel_i, use.names = FALSE)
  out <- frame$units[rows, , drop = FALSE]
  rownames(out) <- NULL
  mh_u <- mh[out$h]
  n_u <- rep(nk, times = nk)
  out$pi_psu <- mh_u / out$M_pop_h                   # pi_k (TC-S2)
  out$pi_unit <- n_u / out$M_k                       # pi_{i|k}
  out$w <- (out$M_pop_h / mh_u) * (out$M_k / n_u)    # U1 form, (TC-S2)
  lead <- c("h", "k", "i", .SAMPLING_ADDED_COLS)
  out <- out[, c(lead, setdiff(names(out), lead)), drop = FALSE]

  attr(out, "design") <- list(
    label = SAMPLING_DEV_LABEL, rng_version = RNG_DRAFT_VERSION, root = root,
    sample_key = sample_key, H = H, J = J, m_h = stats::setNames(as.integer(mh), seq_len(H)),
    n0 = n0, M_pop_h = stats::setNames(frame$M_pop_h, seq_len(H)),
    allow_take_all = allow_take_all, n_take_all = sum(Mk < n0),
    min_pi_unit = min(out$pi_unit),
    spec = "letter/SIM-DESIGN-SPEC.md :339 (U1), :385-386, :514-519, :535-538, :565-568, :1241-1243")
  out
}

# U1 (bit-exact): every unit's w equals (M_pop,h/m_h)(M_k/n0_h) recomputed
# from the CENSUS (pop$psu, looked up by (h, k), PSU counts tabulated afresh) and
# the RECORDED design (attr(sample, "design")), on a code path that shares no
# lookup with sampling_frame(). Also asserts the realised design matches the
# recorded one. Returns list(pass, checks).
check_u1 <- function(sample, pop) {
  d <- attr(sample, "design")
  if (is.null(d)) stop("sample has no recorded design (attr 'design'); call check_u1 on draw_sample() output",
                       call. = FALSE)
  ps <- pop$psu
  ps <- ps[ps$h <= d$H & ps$k <= d$J, c("h", "k", "M")]
  base <- d$J + 1
  Mk <- as.numeric(ps$M[match(sample$h * base + sample$k, ps$h * base + ps$k)])
  Mpop <- as.numeric(table(factor(ps$h, levels = seq_len(d$H))))[sample$h]
  mh <- as.numeric(d$m_h)[sample$h]
  nk <- pmin(d$n0, Mk)
  w_ref <- (Mpop / mh) * (Mk / nk)

  psu_id <- sample$h * base + sample$k
  per_h <- tapply(psu_id, factor(sample$h, levels = seq_len(d$H)), function(x) length(unique(x)))
  per_h[is.na(per_h)] <- 0
  per_psu <- table(psu_id)
  nk_psu <- pmin(d$n0, Mk[match(as.numeric(names(per_psu)), psu_id)])
  inv_prod <- 1 / (sample$pi_psu * sample$pi_unit)
  stage2_sum <- tapply(1 / sample$pi_unit, psu_id, sum)   # HT size of each sampled PSU
  Mk_psu <- Mk[match(as.numeric(names(stage2_sum)), psu_id)]

  rel <- function(a, b) max(abs(a - b) / abs(b))
  checks <- data.frame(
    check = c("census lookup complete",
              "w bit-exact = (M_pop,h/m_h)(M_k/n0) (U1)",
              "pi_psu bit-exact = m_h/M_pop,h",
              "pi_unit bit-exact = n0/M_k",
              "w vs 1/(pi_psu*pi_unit), rel <= 16 eps",
              "m_h PSUs in every stratum",
              "min(n0, M_k) units in every sampled PSU",
              "no duplicated unit",
              "sum_i 1/pi_{i|k} = M_k per PSU, rel <= 64 eps (U2 :340)",
              "recorded M_pop,h matches census"),
    pass = c(!anyNA(Mk),
             !anyNA(Mk) && identical(sample$w, w_ref),
             !anyNA(Mk) && identical(sample$pi_psu, mh / Mpop),
             !anyNA(Mk) && identical(sample$pi_unit, nk / Mk),
             !anyNA(Mk) && rel(sample$w, inv_prod) <= 16 * .Machine$double.eps,
             identical(as.numeric(per_h), as.numeric(d$m_h)),
             !anyNA(Mk) && identical(as.numeric(per_psu), as.numeric(nk_psu)),
             !anyDuplicated(data.frame(sample$h, sample$k, sample$i)),
             !anyNA(Mk) && rel(as.numeric(stage2_sum), Mk_psu) <= 64 * .Machine$double.eps,
             identical(as.numeric(d$M_pop_h), as.numeric(table(factor(ps$h, levels = seq_len(d$H)))))),
    stringsAsFactors = FALSE)
  list(pass = all(checks$pass), checks = checks, label = SAMPLING_DEV_LABEL)
}

# ---- E2 WR code-check cell ------------------------
#
# draw_sample_wr(pop, m_h, n0, key, root): the with-replacement first stage of the E2 WR
# code-check (H = 50, m_h = 2, J = 4, prefixes of the FPC master).
# Stage 1, entity stage1(h): m_h with-replacement selections k_j = floor(J U_j) + 1,
# j = 1..m_h. Stage 2 is keyed by the SELECTION, entity stage2wr(h, j): each selection draws
# M_{k_j} uniforms and keeps the n0 smallest (ties by unit index), so a PSU selected twice gets
# two independent subsamples and every selection is its own ultimate cluster (the
# prototype is toy_sample_wr). Weights are the Hansen-Hurwitz weights
# w = (M_pop,h / m_h)(M_k / n0), in the U1 arithmetic order (TC-S2). Column j is the selection
# index within the stratum: (h, j) is the ultimate cluster that every WR variance estimator
# clusters on. Under this design the UC estimator of the census-centred known-N linearised total
# has expectation EXACTLY sigma2_m,WR = sum_h [(1 - 1/M_pop,h) B_h + V_h] (exact-enumeration
# test in tests/testthat/test-sampling-wr.R).
draw_sample_wr <- function(pop, m_h, n0 = 20, key, root = "SYNTHETIC-DEV-E2E", psu_cols = c("KB", "Cpsu")) {
  .assert_sampling_root(root)
  if (!is.list(key) || is.null(key$cell)) stop("key must be list(cell, population_id, rep_id)", call. = FALSE)
  if (!is.null(key$root) && !identical(key$root, root)) stop("key$root differs from root", call. = FALSE)
  ctx <- key
  ctx$root <- root
  sample_key <- rng_key("sample", ctx)
  cell <- key$cell
  H <- cell$H; J <- cell$J
  frame <- sampling_frame(pop, H, J, psu_cols)
  .whole_vec(n0, "n0")
  if (length(n0) != 1L) stop("n0 must be one number", call. = FALSE)
  mh <- .resolve_m_h(m_h, frame, "m_h")
  if (!identical(as.numeric(mh), as.numeric(.resolve_m_h(cell$m_h, frame, "key$cell$m_h"))))
    stop("m_h differs from key$cell$m_h: the sample stream key must name the allocation drawn", call. = FALSE)
  .whole_vec(mh, "m_h")
  st1 <- rng_states(rng_entity_keys("sample", ctx, "stage1", list(h = seq_len(H))))
  sel_k <- rng_scoped(lapply(seq_len(H), function(h) {
    rng_set_state(st1[h, ])
    floor(J * stats::runif(mh[h])) + 1
  }))
  s_h <- rep(seq_len(H), times = mh)
  s_j <- unlist(lapply(mh, seq_len), use.names = FALSE)
  s_k <- unlist(sel_k, use.names = FALSE)
  if (any(s_k < 1 | s_k > J)) stop("WR selection outside 1..J", call. = FALSE)
  p_idx <- (s_h - 1) * J + s_k
  Mk <- frame$M[p_idx]
  if (any(Mk < n0)) stop("a selected PSU has M_k < n0 = ", n0, ": pi_{i|k} > 1 is infeasible", call. = FALSE)
  st2 <- rng_states(rng_entity_keys("sample", ctx, "stage2wr", list(h = s_h, j = s_j)))
  sel_i <- rng_scoped(lapply(seq_along(p_idx), function(r) {
    rng_set_state(st2[r, ])
    u2 <- stats::runif(Mk[r])
    sort(order(u2, seq_len(Mk[r]))[seq_len(n0)])
  }))
  rows <- rep(frame$offset[p_idx], each = n0) + unlist(sel_i, use.names = FALSE)
  out <- frame$units[rows, , drop = FALSE]
  rownames(out) <- NULL
  out$j <- rep(s_j, each = n0)
  mh_u <- mh[out$h]
  out$pi_psu <- 1 / out$M_pop_h                       # per-draw selection probability p_k (WR)
  out$pi_unit <- n0 / out$M_k
  out$w <- (out$M_pop_h / mh_u) * (out$M_k / n0)      # Hansen-Hurwitz weight, U1 order (TC-S2)
  lead <- c("h", "j", "k", "i", .SAMPLING_ADDED_COLS)
  out <- out[, c(lead, setdiff(names(out), lead)), drop = FALSE]
  attr(out, "design") <- list(
    label = paste(SAMPLING_DEV_LABEL, "| WR first stage, selection-keyed stage 2"), rng_version = RNG_DRAFT_VERSION,
    root = root, sample_key = sample_key, H = H, J = J, m_h = stats::setNames(as.integer(mh), seq_len(H)),
    n0 = n0, M_pop_h = stats::setNames(frame$M_pop_h, seq_len(H)), wr = TRUE,
    n_repeat_selections = sum(duplicated(data.frame(s_h, s_k))), min_pi_unit = min(out$pi_unit),
    spec = "letter/SIM-DESIGN-SPEC.md :655-656, :1143; RNG draft §7 (stage2wr); rulings R13, R14")
  out
}

# U1 for the WR sample (bit-exact class): w = (M_pop,h/m_h)(M_k/n0) recomputed from the census and
# the recorded design on a separate path; m_h selections per stratum, n0 distinct units per
# selection, every unit inside its selected PSU, sum_i 1/pi_{i|k} = M_k per selection. Units may
# repeat ACROSS selections of the same PSU (independent subsamples), never within one.
check_u1_wr <- function(sample, pop) {
  d <- attr(sample, "design")
  if (is.null(d) || !isTRUE(d$wr)) stop("check_u1_wr needs draw_sample_wr() output", call. = FALSE)
  ps <- pop$psu[pop$psu$h <= d$H & pop$psu$k <= d$J, c("h", "k", "M")]
  base <- d$J + 1
  Mk <- as.numeric(ps$M[match(sample$h * base + sample$k, ps$h * base + ps$k)])
  Mpop <- as.numeric(table(factor(ps$h, levels = seq_len(d$H))))[sample$h]
  mh <- as.numeric(d$m_h)[sample$h]
  w_ref <- (Mpop / mh) * (Mk / d$n0)
  sel <- paste(sample$h, sample$j)
  per_h <- tapply(sample$j, factor(sample$h, levels = seq_len(d$H)), function(x) length(unique(x)))
  per_h[is.na(per_h)] <- 0
  per_sel <- as.vector(table(sel))
  k_const <- all(tapply(sample$k, sel, function(x) length(unique(x))) == 1L)
  stage2_sum <- tapply(1 / sample$pi_unit, sel, sum)
  Mk_sel <- tapply(Mk, sel, `[`, 1)
  rel <- function(a, b) max(abs(a - b) / abs(b))
  checks <- data.frame(
    check = c("census lookup complete",
              "w bit-exact = (M_pop,h/m_h)(M_k/n0) (U1, WR)",
              "pi_unit bit-exact = n0/M_k",
              "m_h selections in every stratum",
              "n0 units in every selection",
              "one PSU per selection, inside the frame (1 <= k <= J)",
              "no duplicated unit within a selection",
              "sum_i 1/pi_{i|k} = M_k per selection, rel <= 64 eps",
              "recorded M_pop,h matches census"),
    pass = c(!anyNA(Mk),
             !anyNA(Mk) && identical(sample$w, w_ref),
             !anyNA(Mk) && identical(sample$pi_unit, d$n0 / Mk),
             identical(as.numeric(per_h), as.numeric(d$m_h)),
             all(per_sel == d$n0),
             k_const && all(sample$k >= 1 & sample$k <= d$J),
             !anyDuplicated(data.frame(sample$h, sample$j, sample$i)),
             !anyNA(Mk) && rel(as.numeric(stage2_sum), as.numeric(Mk_sel)) <= 64 * .Machine$double.eps,
             identical(as.numeric(d$M_pop_h), as.numeric(table(factor(ps$h, levels = seq_len(d$H)))))),
    stringsAsFactors = FALSE)
  list(pass = all(checks$pass), checks = checks, label = SAMPLING_DEV_LABEL)
}
