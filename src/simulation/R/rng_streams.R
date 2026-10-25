# The random-number contract of the simulation: key, hash and fill conventions for every
# random draw of the dry and production runs and of the tests.
#
# Production populations are built by R/dgp.R under the adopted production root SIM_PROD_ROOT; the
# toy generator below stays SYNTHETIC-only, and its variables are placeholders that are
# deliberately NOT the DGP (no Q0/g0 coefficients, no K_B class variances, no GH
# marginalisation). assert_engine_root() accepts exactly two root families: "SYNTHETIC..." (tests,
# DEV) and SIM_PROD_ROOT (dry and production runs).
#
# Mechanism: canonical ASCII key -> SHA-256 -> 6-word MRG32k3a state -> R's
# L'Ecuyer-CMRG runif, inside a scope that restores the caller's RNG. Population
# draws are addressed per PSU, so a PSU's values never depend on how many other
# strata or PSUs are generated: that is what makes every prefix view of a master
# identical to direct generation at the smaller size (prefix rule).

RNG_DRAFT_VERSION <- "rng-draft-2"            # the contract version recorded in every run
SIM_PROD_ROOT <- "svytmle-eif-clt-sim-v9.7-rng1"   # the production root

# The two root families the engine accepts: SYNTHETIC... for tests and DEV runs, and the
# production root for dry and production runs.
assert_engine_root <- function(root, what = "engine") {
  if (!is.character(root) || length(root) != 1L || is.na(root) ||
      !(startsWith(root, "SYNTHETIC") || identical(root, SIM_PROD_ROOT)))
    stop(what, ": root must start with 'SYNTHETIC' or be the adopted production root '", SIM_PROD_ROOT, "'",
         call. = FALSE)
  invisible(root)
}
is_prod_root <- function(root) identical(root, SIM_PROD_ROOT)

# MRG32k3a moduli, the generator behind R's "L'Ecuyer-CMRG" (?RNGkind).
.MRG_M1 <- 4294967087
.MRG_M2 <- 4294944443
# R's normalisation constant for that generator; used only by the reference
# emulation that pins runif() bit-for-bit (rng_mrg32k3a_emulate).
.MRG_NORM <- 2.328306549295727688e-10

# The fields each level may depend on: exactly those named by the hierarchical
# contract plus the fifth, post-processing level. Level 4
# keys the fold INDEX, not fold content. Context fields not listed
# for a level are ignored by that level, which is what the separation tests vary.
RNG_LEVEL_FIELDS <- list(
  population  = c("master_id", "population_id"),
  sample      = c("cell", "population_id", "rep_id"),
  fold        = c("cell", "population_id", "rep_id", "fold_structure_id"),
  learner     = c("cell", "population_id", "rep_id", "arm", "fold"),
  postprocess = c("cell", "statistic", "grading_pass")
)

# cell := (regime, H, m_h, J, m, f, E2-construction). V and alpha are not
# cell fields, so a cell carrying either is rejected, not ignored.
RNG_CELL_FIELDS <- c("regime", "H", "m_h", "J", "m", "f", "e2c")

# Entity sub-addresses below a level key. They are coordinates inside the
# keyed object (a master's frame, a cell's sample), never cell or rep factors.
# stage2wr(h,j) addresses the j-th with-replacement SELECTION of the WR cell, not
# a PSU, so a PSU selected twice gets two independent subsamples (docs).
RNG_ENTITY_COORDS <- list(psu = c("h", "k"), stage1 = "h", stage2 = c("h", "k"),
                          stage2wr = c("h", "j"), stratum = "h")
RNG_LEVEL_ENTITIES <- list(population = "psu",
                           sample = c("stage1", "stage2", "stage2wr"),
                           fold = "stratum", learner = character(0),
                           postprocess = character(0))

.LABEL_RE <- "^[A-Za-z0-9._+-]+$"

# ---- canonical key encoding ---------------------------------------------------

# Whole numbers print in plain decimal; floating values never enter a key, because
# their decimal printing is not a canonical byte string across platforms.
.whole_msg <- function(what)
  paste0("key field '", what, "' must be a non-negative whole number ",
         "(floating values never enter a key)")

.enc_whole_vec <- function(x, what) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x != floor(x)) ||
      any(x < 0) || any(x >= 2^53))
    stop(.whole_msg(what), call. = FALSE)
  sprintf("%.0f", x)
}

.enc_whole <- function(x, what) {
  if (length(x) != 1L) stop(.whole_msg(what), call. = FALSE)
  .enc_whole_vec(x, what)
}

.enc_label <- function(x, what) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !grepl(.LABEL_RE, x))
    stop("key field '", what, "' must be one ASCII label matching ", .LABEL_RE,
         call. = FALSE)
  x
}

.enc_scalar <- function(x, what) {
  if (is.character(x)) .enc_label(x, what) else .enc_whole(x, what)
}

.gcd <- function(a, b) if (b == 0) a else .gcd(b, a %% b)

# f is an exact reduced rational "p/q", never a decimal (0.1 has no exact binary
# form); it must equal m_h/J where both are known.
.enc_rational <- function(x, what) {
  if (!is.character(x) || length(x) != 1L || !grepl("^[0-9]+/[0-9]+$", x))
    stop("key field '", what, "' must be a rational string 'p/q'", call. = FALSE)
  pq <- as.numeric(strsplit(x, "/", fixed = TRUE)[[1]])
  if (pq[2] == 0 || .gcd(pq[1], pq[2]) != 1)
    stop("key field '", what, "' must be a reduced rational with q > 0", call. = FALSE)
  x
}

# Class-indexed values (the unequal-allocation pair has m_h and f by K_B class)
# encode as "(K1=..,K2=..)" in C-locale name order.
.enc_maybe_classed <- function(x, what, enc) {
  if (length(x) == 1L && is.null(names(x))) return(enc(x, what))
  nm <- names(x)
  if (is.null(nm) || anyNA(nm) || !all(grepl("^K[0-9]+$", nm)) || anyDuplicated(nm))
    stop("class-indexed '", what, "' needs unique names K1, K2, ...", call. = FALSE)
  o <- order(nm, method = "radix")
  parts <- vapply(o, function(j) paste0(nm[j], "=", enc(x[[j]], what)), "")
  paste0("(", paste(parts, collapse = ","), ")")
}

rng_encode_cell <- function(cell) {
  if (!is.list(cell) || is.null(names(cell)))
    stop("cell must be a named list of the :565 fields", call. = FALSE)
  extra <- setdiff(names(cell), RNG_CELL_FIELDS)
  if (length(extra))
    stop("cell carries field(s) outside :565: ", paste(extra, collapse = ", "),
         " (V and alpha never enter cell, :565-568)", call. = FALSE)
  miss <- setdiff(RNG_CELL_FIELDS, names(cell))
  if (length(miss))
    stop("cell lacks field(s): ", paste(miss, collapse = ", "), call. = FALSE)
  .check_f(cell)
  enc <- c(
    regime = .enc_label(cell$regime, "regime"),
    H      = .enc_whole(cell$H, "H"),
    m_h    = .enc_maybe_classed(cell$m_h, "m_h", .enc_whole),
    J      = .enc_whole(cell$J, "J"),
    m      = .enc_whole(cell$m, "m"),
    f      = .enc_maybe_classed(cell$f, "f", .enc_rational),
    e2c    = .enc_label(cell$e2c, "e2c"))
  paste0("(", paste0(names(enc), "=", enc, collapse = ","), ")")
}

# f = m_h/J holds in every cell; checked exactly in integers.
.check_f <- function(cell) {
  mh <- cell$m_h
  f <- cell$f
  if (length(mh) != length(f)) stop("m_h and f must be indexed alike", call. = FALSE)
  if (length(mh) > 1L) {
    if (!setequal(names(mh), names(f))) stop("m_h and f class names differ", call. = FALSE)
    f <- f[names(mh)]
  }
  for (j in seq_along(mh)) {
    pq <- as.numeric(strsplit(.enc_rational(f[[j]], "f"), "/", fixed = TRUE)[[1]])
    if (pq[1] * cell$J != pq[2] * mh[[j]])
      stop("f = ", f[[j]], " is not m_h/J = ", mh[[j]], "/", cell$J, call. = FALSE)
  }
  invisible(TRUE)
}

.enc_entity <- function(level, entity) {
  if (is.null(entity)) return(NULL)
  type <- entity$type
  if (!is.character(type) || !type %in% RNG_LEVEL_ENTITIES[[level]])
    stop("entity type '", format(type), "' is not addressable at level ", level,
         call. = FALSE)
  co <- RNG_ENTITY_COORDS[[type]]
  if (!setequal(setdiff(names(entity), "type"), co))
    stop("entity '", type, "' needs exactly coordinates ", paste(co, collapse = ","),
         call. = FALSE)
  vals <- vapply(co, function(nm) .enc_whole(entity[[nm]], nm), "")
  paste0(type, "(", paste0(co, "=", vals, collapse = ","), ")")
}

# Canonical key: "<root>;L=<level>;<field>=<value>;...[;ent=<entity>]", fields
# in the fixed RNG_LEVEL_FIELDS order. Pure ASCII, so its bytes are the same on
# every platform and locale. The level name is inside the key, so keys of
# different levels can never coincide.
rng_key <- function(level, ctx, entity = NULL) {
  level <- match.arg(level, names(RNG_LEVEL_FIELDS))
  root <- .enc_label(ctx$root, "root")
  fields <- RNG_LEVEL_FIELDS[[level]]
  miss <- fields[vapply(fields, function(f) is.null(ctx[[f]]), TRUE)]
  if (length(miss))
    stop("level ", level, " needs field(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  vals <- vapply(fields, function(f)
    if (f == "cell") rng_encode_cell(ctx$cell) else .enc_scalar(ctx[[f]], f), "")
  parts <- c(root, paste0("L=", level), paste0(fields, "=", vals))
  ent <- .enc_entity(level, entity)
  if (!is.null(ent)) parts <- c(parts, paste0("ent=", ent))
  paste(parts, collapse = ";")
}

# Vectorised rng_key() for many entities of one level key: the level prefix is
# encoded once, the coordinates with the same rules. Must equal rng_key() per
# entity (tested).
rng_entity_keys <- function(level, ctx, type, coords) {
  level <- match.arg(level, names(RNG_LEVEL_FIELDS))
  if (!is.character(type) || length(type) != 1L || !type %in% RNG_LEVEL_ENTITIES[[level]])
    stop("entity type '", format(type), "' is not addressable at level ", level,
         call. = FALSE)
  co <- RNG_ENTITY_COORDS[[type]]
  if (!setequal(names(coords), co))
    stop("entity '", type, "' needs exactly coordinates ", paste(co, collapse = ","),
         call. = FALSE)
  inner <- paste0(co[1], "=", .enc_whole_vec(coords[[co[1]]], co[1]))
  for (nm in co[-1]) inner <- paste0(inner, ",", nm, "=", .enc_whole_vec(coords[[nm]], nm))
  paste0(rng_key(level, ctx), ";ent=", type, "(", inner, ")")
}

# ---- hash -> MRG32k3a state -----------------------------------------------------

# SHA-256 of the key's bytes, never of R's serialize() output (whose format is
# R-version dependent): serialize = FALSE hashes the string itself.
.sha256_keys <- local({
  vd <- NULL
  function(keys) {
    if (is.null(vd)) vd <<- digest::getVDigest(algo = "sha256")
    vd(keys, serialize = FALSE)
  }
})

# 64 hex digits -> 8 big-endian 32-bit words, held exactly as doubles.
.hex_words <- function(hex) {
  n <- length(hex)
  starts <- seq(1L, 61L, by = 4L)
  hw <- strtoi(substring(rep(hex, each = 16L), starts, starts + 3L), base = 16L)
  hw <- matrix(hw, nrow = n, byrow = TRUE)
  hw[, seq(1L, 15L, 2L), drop = FALSE] * 65536 + hw[, seq(2L, 16L, 2L), drop = FALSE]
}

# Words 1-3 mod m1, 4-6 mod m2 (words 7-8 unused). A state R would not accept is
# flagged: an all-zero component triple, or a word equal to 2^31, which has no
# non-NA signed-integer form in .Random.seed. R silently re-seeds from the clock
# on an out-of-range L'Ecuyer state, so validity must be enforced, not assumed.
rng_state_from_words <- function(w) {
  w <- matrix(w, ncol = 8L)
  s <- cbind(w[, 1:3, drop = FALSE] %% .MRG_M1, w[, 4:6, drop = FALSE] %% .MRG_M2)
  valid <- rowSums(s[, 1:3, drop = FALSE] != 0) > 0 &
    rowSums(s[, 4:6, drop = FALSE] != 0) > 0 & rowSums(s == 2^31) == 0
  list(state = s, valid = valid)
}

# Keys -> n x 6 state matrix. Invalid states (probability ~ 1e-9 per key) are
# re-derived from "<key>#retry=<r>", r = 1, 2, ..., which depends on that key alone.
rng_states <- function(keys) {
  out <- rng_state_from_words(.hex_words(.sha256_keys(keys)))
  for (j in which(!out$valid)) {
    r <- 0L
    repeat {
      r <- r + 1L
      one <- rng_state_from_words(.hex_words(.sha256_keys(paste0(keys[j], "#retry=", r))))
      if (one$valid) break
    }
    out$state[j, ] <- one$state
  }
  out$state
}

rng_state <- function(key) drop(rng_states(key))

.to_seed_int <- function(s) as.integer(ifelse(s >= 2^31, s - 2^32, s))

# ---- scoped streams -------------------------------------------------------------

.rng_scope <- new.env(parent = emptyenv())

# Runs `expr` under L'Ecuyer-CMRG and restores the caller's RNG kind and
# .Random.seed afterwards, including on error, and including the case where the
# caller had no .Random.seed at all (the "scoped stream that is restored
# afterwards"). Not restored: a Box-Muller caller's cached second normal,
# which any RNGkind switch clears; production uses Inversion, where this is exact.
rng_scoped <- function(expr) {
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = genv, inherits = FALSE)
  old_kind <- RNGkind()
  old_depth <- if (is.null(.rng_scope$depth)) 0L else .rng_scope$depth
  old_code <- .rng_scope$code
  on.exit({
    suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
    .rng_scope$depth <- old_depth
    .rng_scope$code <- old_code
  }, add = TRUE)
  suppressWarnings(RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"))
  .rng_scope$code <- get(".Random.seed", envir = genv)[1L]
  .rng_scope$depth <- old_depth + 1L
  force(expr)
}

# Positions the scoped generator at an entity's state. Only legal inside
# rng_scoped(), so a stream can never leak into the caller's RNG.
rng_set_state <- function(state) {
  if (is.null(.rng_scope$depth) || .rng_scope$depth < 1L)
    stop("rng_set_state() must be called inside rng_scoped()", call. = FALSE)
  assign(".Random.seed", c(.rng_scope$code, .to_seed_int(state)), envir = globalenv())
  invisible(NULL)
}

# n uniforms from one level stream (optionally one entity of it).
rng_uniforms <- function(level, ctx, n, entity = NULL) {
  st <- rng_state(rng_key(level, ctx, entity))
  rng_scoped({
    rng_set_state(st)
    stats::runif(n)
  })
}

# Plain-double reference of MRG32k3a x normc. Every product is below 2^53, so the
# arithmetic is exact; it reproduces runif() bit-for-bit, which pins the
# generator independently of R's internals and ports directly to other languages.
rng_mrg32k3a_emulate <- function(state, n) {
  s <- as.numeric(state)
  out <- numeric(n)
  for (i in seq_len(n)) {
    p1 <- (1403580 * s[2] - 810728 * s[1]) %% .MRG_M1
    s[1:3] <- c(s[2], s[3], p1)
    p2 <- (527612 * s[6] - 1370589 * s[4]) %% .MRG_M2
    s[4:6] <- c(s[5], s[6], p2)
    out[i] <- (if (p1 > p2) p1 - p2 else p1 - p2 + .MRG_M1) * .MRG_NORM
  }
  out
}

# ---- rounding -------------------------------------------------------------------

# Round half to even on the binary64 value, stated explicitly rather than relying
# on round(). Exact for x >= 0: x - floor(x) is exact there (Sterbenz).
round_half_even <- function(x) {
  if (any(!is.finite(x)) || any(x < 0))
    stop("round_half_even() is defined here for finite x >= 0", call. = FALSE)
  r <- floor(x)
  d <- x - r
  out <- r + (d > 0.5)
  tie <- d == 0.5
  out[tie] <- r[tie] + (r[tie] %% 2 != 0)
  out
}

# Size from a driver: centre + slope * c evaluated as two binary64 operations in
# that order (slope * c first), then round_half_even. The law is
# the case centre = 34, slope = 14; the toy below uses a placeholder law.
size_from_driver <- function(c, centre, slope) round_half_even(centre + slope * c)

# ---- SYNTHETIC toy generator (NOT the DGP) ------------------------------

# Placeholder constants, chosen only to exercise variable PSU sizes and binary
# thresholds; they correspond to nothing in.
TOY_SIZE <- c(centre = 6, slope = 4)   # sizes in {2, ..., 10}
TOY_COEF <- list(a = c(0.2, 0.5, -0.4, 0.3), y = c(-0.1, 0.6, 0.4, 0.5))

TOY_IDENTITY_COLS <- c("h", "k", "i", "M", "c_psu", "u_za", "u_zy",
                       "x1", "x2", "x3", "x4", "va", "vy")
TOY_DERIVED_COLS <- c("za", "zy", "pa", "a", "py", "y", "nu1", "nu0")

.assert_synthetic <- function(root) {
  if (!is.character(root) || length(root) != 1L || !startsWith(root, "SYNTHETIC"))
    stop("toy generator: root must start with 'SYNTHETIC' (no production populations)",
         call. = FALSE)
}

# One PSU's coefficient-free draws, in this within-entity order:
# [driver, PSU-effect uniform A, PSU-effect uniform Y], then 6 x M uniforms read
# column-major (variable blocks x1..x4, then the A and Y threshold uniforms).
# The draw count depends on this PSU's own M only.
.toy_psu <- function(h, k) {
  u0 <- stats::runif(3L)
  c_psu <- -1 + 2 * u0[1]
  M <- as.integer(size_from_driver(c_psu, TOY_SIZE[["centre"]], TOY_SIZE[["slope"]]))
  U <- matrix(stats::runif(6L * M), nrow = M)
  list(h = rep(as.integer(h), M), k = rep(as.integer(k), M), i = seq_len(M),
       M = rep(M, M), c_psu = rep(c_psu, M), u_za = rep(u0[2], M),
       u_zy = rep(u0[3], M), x1 = -1 + 2 * U[, 1], x2 = -1 + 2 * U[, 2],
       x3 = -1 + 2 * U[, 3], x4 = -1 + 2 * U[, 4], va = U[, 5], vy = U[, 6])
}

# SYNTHETIC toy population of a H x J frame. `traversal` permutes the order in
# which PSUs are generated; with entity streams the result cannot depend on it.
toy_population <- function(ctx, H, J, coef = TOY_COEF, traversal = NULL) {
  .assert_synthetic(ctx$root)
  hk <- data.frame(h = rep(seq_len(H), each = J), k = rep(seq_len(J), times = H))
  if (!is.null(traversal)) hk <- hk[traversal, , drop = FALSE]
  states <- rng_states(rng_entity_keys("population", ctx, "psu", hk))
  rows <- rng_scoped(lapply(seq_len(nrow(hk)), function(r) {
    rng_set_state(states[r, ])
    .toy_psu(hk$h[r], hk$k[r])
  }))
  pop <- as.data.frame(lapply(stats::setNames(nm = names(rows[[1]])), function(cl)
    unlist(lapply(rows, `[[`, cl), use.names = FALSE)))
  pop <- pop[order(pop$h, pop$k, pop$i), , drop = FALSE]
  rownames(pop) <- NULL
  toy_derive(pop, coef)
}

# Coefficient-dependent layer: thresholds the stored uniforms (inversion), so a
# coefficient change can flip A or Y but never moves a draw. nu1/nu0 stand in for
# "oracle nuisances": deterministic placeholder functions, not GH marginals.
toy_derive <- function(pop, coef = TOY_COEF) {
  pop$za <- stats::qnorm(pop$u_za)
  pop$zy <- stats::qnorm(pop$u_zy)
  ca <- coef$a
  cy <- coef$y
  pop$pa <- stats::plogis(ca[1] + ca[2] * pop$x1 + ca[3] * pop$x2 + ca[4] * pop$za)
  pop$a <- as.integer(pop$va < pop$pa)
  pop$py <- stats::plogis(cy[1] + cy[2] * pop$a + cy[3] * pop$x3 + cy[4] * pop$zy)
  pop$y <- as.integer(pop$vy < pop$py)
  pop$nu1 <- stats::plogis(cy[1] + cy[2] + cy[3] * pop$x3)
  pop$nu0 <- stats::plogis(cy[1] + cy[3] * pop$x3)
  pop
}

# A prefix view: the first H strata and the first J PSUs of each.
toy_prefix_view <- function(pop, H, J) {
  v <- pop[pop$h <= H & pop$k <= J, , drop = FALSE]
  rownames(v) <- NULL
  v
}

# SYNTHETIC sample-stream use (level 2): per-stratum stage-1 SRS-WOR of m_h of
# the first J PSUs by ranking J uniforms (ties broken by frame index k), then
# per-PSU stage-2 SRS-WOR of n0 units the same way. Explicit ranking keeps the
# draw count fixed and avoids sample()'s version-dependent algorithm. `pop` may be
# the whole master: only the cell's frame (first H strata, first J PSUs) is used.
toy_sample <- function(pop, ctx, n0) {
  .assert_synthetic(ctx$root)
  H <- ctx$cell$H
  J <- ctx$cell$J
  m_h <- ctx$cell$m_h
  if (length(m_h) != 1L) stop("toy_sample() takes a scalar m_h", call. = FALSE)
  sel <- vector("list", H)
  rng_scoped(for (h in seq_len(H)) {
    rng_set_state(rng_state(rng_key("sample", ctx, list(type = "stage1", h = h))))
    u <- stats::runif(J)
    ks <- sort(order(u, seq_len(J))[seq_len(m_h)])
    units <- lapply(ks, function(k) {
      Mk <- pop$M[pop$h == h & pop$k == k][1]
      rng_set_state(rng_state(rng_key("sample", ctx, list(type = "stage2", h = h, k = k))))
      u2 <- stats::runif(Mk)
      data.frame(h = h, k = k, i = sort(order(u2, seq_len(Mk))[seq_len(min(n0, Mk))]))
    })
    sel[[h]] <- do.call(rbind, units)
  })
  out <- do.call(rbind, sel)
  rownames(out) <- NULL
  out
}

# SYNTHETIC level-2 use for the WR code-check cell. Stage 1:
# m_h PSU indices floor(J*U)+1, with replacement, from entity stage1(h). Stage 2
# is keyed by the selection j (entity stage2wr(h,j)), never by the PSU: keyed by
# (h,k), a PSU selected twice would give both selections the same units. Selection
# keying makes se^2 unbiased for the literal WR variance. Whether that is the
# check's target: census_wr_target() in R/census.R gives that literal WR variance, which
# differs from sigma_m^2 (Supplement S1) by sum_h B_h/M_pop,h.
toy_sample_wr <- function(pop, ctx, n0) {
  .assert_synthetic(ctx$root)
  H <- ctx$cell$H
  J <- ctx$cell$J
  m_h <- ctx$cell$m_h
  if (length(m_h) != 1L) stop("toy_sample_wr() takes a scalar m_h", call. = FALSE)
  sel <- vector("list", H)
  rng_scoped(for (h in seq_len(H)) {
    rng_set_state(rng_state(rng_key("sample", ctx, list(type = "stage1", h = h))))
    ks <- floor(J * stats::runif(m_h)) + 1
    sel[[h]] <- do.call(rbind, lapply(seq_len(m_h), function(j) {
      Mk <- pop$M[pop$h == h & pop$k == ks[j]][1]
      rng_set_state(rng_state(rng_key("sample", ctx, list(type = "stage2wr", h = h, j = j))))
      u2 <- stats::runif(Mk)
      data.frame(h = h, j = j, k = ks[j],
                 i = sort(order(u2, seq_len(Mk))[seq_len(min(n0, Mk))]))
    }))
  })
  out <- do.call(rbind, sel)
  rownames(out) <- NULL
  out
}

# ---- population-store bytes and hash ------------------------------------------

# Canonical bytes: an ASCII header naming tag, columns, types and row count, then
# each column in canonical (h, k, i) row order as little-endian int32 ("i4") or
# IEEE binary64 ("f8"). No saveRDS: its bytes depend on R's serialization format
# and compression. NA/NaN are refused, so no payload ambiguity can arise.
rng_store_bytes <- function(df, cols, tag) {
  ord <- order(df$h, df$k, df$i)
  types <- vapply(cols, function(cl) if (is.integer(df[[cl]])) "i4" else "f8", "")
  header <- paste0(tag, ";cols=", paste(cols, collapse = ","), ";types=",
                   paste(types, collapse = ","), ";nrow=", nrow(df), "\n")
  body <- lapply(cols, function(cl) {
    x <- df[[cl]][ord]
    if (anyNA(x)) stop("store column '", cl, "' contains NA/NaN", call. = FALSE)
    # "+ 0" maps -0 to +0 so the two zeros cannot give different bytes.
    if (is.integer(x)) writeBin(x, raw(), size = 4L, endian = "little")
    else writeBin(as.double(x) + 0, raw(), size = 8L, endian = "little")
  })
  c(charToRaw(header), unlist(body, use.names = FALSE))
}

rng_store_hash <- function(df, cols, tag) {
  digest::digest(rng_store_bytes(df, cols, tag), algo = "sha256", serialize = FALSE)
}
