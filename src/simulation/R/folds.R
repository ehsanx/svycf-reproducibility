# Strict whole-PSU cross-fitting folds for this simulation only.
#
# Contract: G8/U4 -- V | m_h and m_{v,h} = m_h/V
# asserted per stratum, an incompatible request is a hard error, never a reduced V;
# G8/U5 -- training/scoring PSU disjointness, PSU-CF arms only;
# V fixed, whole-PSU balanced folds; G3(c)
# integrality; level 3 -- the fold stream is keyed on exactly
# (cell, population id, rep id, fold-structure id).
#
# This is NOT the shared svytmle make_cf_folds(): its automatic V reduction and
# approximate rep_len balance are intentional for the package and SPA, and are
# exactly what U4/U5 forbid here. Nothing below is
# inherited from it.
#
# Randomness enters only through the caller-supplied `draw`; the fold-stream
# keying and traversal are supplied by the caller (R/estimators.R with
# R/rng_streams.R), so no convention
# is fixed in this file.

# The level-3 key factors. V is deliberately absent: the V-
# sensitivity cell must see the same fold stream at V=2 and V=5.
.fold_key_fields <- c("cell", "population_id", "rep_id", "fold_structure_id")

# Factors that may enter only at level 4 (arm, alpha, fold) or
# never enter `cell` at all (V).
.fold_key_forbidden <- c("v", "arm", "alpha", "fold")

.check_fold_key <- function(fold_key) {
  nm <- names(fold_key)
  if (!is.list(fold_key) || is.null(nm) || anyDuplicated(nm) ||
      !setequal(nm, .fold_key_fields) || length(fold_key) != length(.fold_key_fields))
    stop("fold_key must be a list with exactly the fields ",
         paste(.fold_key_fields, collapse = ", "),
         " (fold stream = hash(cell, population id, rep id, fold-structure id), ",
         "the fold-assignment random stream)", call. = FALSE)
  for (f in .fold_key_fields) {
    x <- fold_key[[f]]
    if (!length(x) || anyNA(unlist(x, use.names = FALSE)))
      stop("fold_key$", f, " is empty or contains NA", call. = FALSE)
  }
  for (f in c("population_id", "rep_id", "fold_structure_id"))
    if (length(fold_key[[f]]) != 1L)
      stop("fold_key$", f, " must be a single id", call. = FALSE)
  bad <- intersect(tolower(names(fold_key$cell)), .fold_key_forbidden)
  if (length(bad))
    stop("fold_key$cell must not contain ", paste(bad, collapse = ", "),
         ": V is not part of cell, and arm/alpha/fold ",
         "enter only the level-4 learner stream", call. = FALSE)
  invisible(fold_key)
}

.check_V <- function(V) {
  if (!is.numeric(V) || length(V) != 1L || !is.finite(V) || V != round(V) || V < 2)
    stop("V must be a single whole number >= 2; it is fixed by the simulation cell ",
         "and is never adapted", call. = FALSE)
  as.integer(V)
}

# PSUs are keyed on the (stratum, psu) pair so ids that repeat across strata are
# distinct PSUs. Integer codes, not pasted strings, so no separator can collide.
# Canonical PSU order is a locale-independent radix sort on (stratum, psu): the
# assignment then cannot depend on the row order of the sample.
.psu_index <- function(strata, psu) {
  s_lev <- unique(strata)
  p_lev <- unique(psu)
  pair <- match(strata, s_lev) * (length(p_lev) + 1) + match(psu, p_lev)
  first <- which(!duplicated(pair))
  tab <- data.frame(stratum = strata[first], psu = psu[first],
                    stringsAsFactors = FALSE)
  ord <- order(tab$stratum, tab$psu, method = "radix")
  tab <- tab[ord, , drop = FALSE]
  rownames(tab) <- NULL
  list(tab = tab, unit_psu = match(pair, pair[first][ord]))
}

# Restores the caller's RNG exactly (seed vector, or the uninitialised state plus
# kinds), even on error. The fold stream must never consume or reseed the
# ambient stream that the population/sample layers use.
.with_preserved_rng <- function(expr) {
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had_seed) seed <- get(".Random.seed", envir = genv, inherits = FALSE)
  kinds <- RNGkind()
  on.exit({
    if (had_seed) {
      assign(".Random.seed", seed, envir = genv)
    } else {
      suppressWarnings(RNGkind(kinds[1L], kinds[2L], kinds[3L]))
      if (exists(".Random.seed", envir = genv, inherits = FALSE))
        rm(".Random.seed", envir = genv)
    }
  }, add = TRUE)
  expr
}

.valid_scores <- function(s, n) {
  if (!is.numeric(s) || length(s) != n || any(!is.finite(s)))
    stop("draw() must return ", n, " finite numeric scores, one per PSU of the ",
         "canonical PSU table it was given", call. = FALSE)
  as.numeric(unname(s))
}

# The purity check runs the draw under two different ambient states: a draw that
# reads the global stream instead of its key gives different scores and is
# rejected. The two set.seed() calls only perturb the ambient state inside the
# restored scope; they cannot reach the result of a pure draw.
.draw_scores <- function(draw, fold_key, psus, check_purity) {
  n <- nrow(psus)
  .with_preserved_rng({
    if (check_purity) {
      set.seed(1L)
      s1 <- .valid_scores(draw(fold_key, psus), n)
      set.seed(2L)
      s2 <- .valid_scores(draw(fold_key, psus), n)
      if (!identical(s1, s2))
        stop("draw() depends on the ambient RNG state; fold randomness must come ",
             "only from fold_key", call. = FALSE)
      s1
    } else {
      .valid_scores(draw(fold_key, psus), n)
    }
  })
}

.report_strata <- function(m_h, which_bad, max_show = 10L) {
  b <- which(which_bad)
  shown <- utils::head(b, max_show)
  txt <- paste0("'", names(m_h)[shown], "' (m_h = ", m_h[shown], ")", collapse = ", ")
  if (length(b) > max_show) txt <- paste0(txt, ", ... ", length(b) - max_show, " more")
  txt
}

#' Strict whole-PSU folds (simulation only).
#'
#' @param strata,psu unit-level stratum and PSU ids of the realised sample.
#' @param V fixed fold count of the cell (: V=2 everywhere except the E3 V=5 limb).
#' @param fold_key list(cell, population_id, rep_id, fold_structure_id), the
#'   level-3 key. Passed to `draw`, recorded on the result.
#' @param draw function(fold_key, psus) returning one finite score per row of
#'   `psus` (data.frame(stratum, psu) in canonical order). This is the pluggable
#'   fold-stream point; the production convention is OPEN. An attribute
#'   "convention" on `draw` is recorded on the result.
#' @param m_h_design optional design m_h (scalar, or named by stratum) asserted
#'   against the realised PSU counts.
#' @param check_draw_purity reject a draw that depends on ambient RNG state.
#' @return integer fold labels in 1..V, one per unit, with attributes
#'   fold_scheme, V, m_h, m_vh, psu_table, fold_key, draw_convention.
make_strict_psu_folds <- function(strata, psu, V, fold_key, draw,
                                  m_h_design = NULL, check_draw_purity = TRUE) {
  n <- length(strata)
  if (!n || length(psu) != n)
    stop("strata and psu must be non-empty and of equal length", call. = FALSE)
  if (anyNA(strata) || anyNA(psu))
    stop("NA stratum or PSU id: every sampled PSU must carry exactly one fold ",
         "label", call. = FALSE)
  V <- .check_V(V)
  .check_fold_key(fold_key)
  if (!is.function(draw)) stop("draw must be a function(fold_key, psus)", call. = FALSE)

  idx <- .psu_index(strata, psu)
  tab <- idx$tab
  s_code <- match(tab$stratum, unique(tab$stratum))
  m_h <- tabulate(s_code)
  names(m_h) <- as.character(unique(tab$stratum))

  # U4 and : V stays fixed and an incompatible request stops.
  # This is the point where the shared package would silently shrink V instead.
  bad <- m_h %% V != 0L
  if (any(bad))
    stop("incompatible fold request: V = ", V, " does not divide m_h in ", sum(bad),
         " stratum/strata: ", .report_strata(m_h, bad),
         ". V | m_h is required per stratum and V is never reduced ",
         "(folds are exactly balanced)", call. = FALSE)
  if (!is.null(m_h_design)) .check_m_h_design(m_h, m_h_design)
  m_vh <- m_h %/% V

  score <- .draw_scores(draw, fold_key, tab[, c("stratum", "psu")], check_draw_purity)

  # Within each stratum, rank PSUs by score (ties broken by canonical PSU order,
  # so the result is a pure function of the scores) and give contiguous blocks of
  # exactly m_{v,h} = m_h/V whole PSUs to folds 1..V.
  o <- order(s_code, score, seq_along(score), method = "radix")
  r <- sequence(m_h)
  pfold <- integer(nrow(tab))
  pfold[o] <- (r - 1L) %/% m_vh[s_code[o]] + 1L
  tab$fold <- pfold

  fold <- pfold[idx$unit_psu]
  attr(fold, "fold_scheme") <- "strict_psu"
  attr(fold, "V") <- V
  attr(fold, "m_h") <- m_h
  attr(fold, "m_vh") <- m_vh
  attr(fold, "psu_table") <- tab
  attr(fold, "fold_key") <- fold_key
  attr(fold, "draw_convention") <- if (is.null(attr(draw, "convention")))
    NA_character_ else as.character(attr(draw, "convention"))
  # Recount on an independent path before returning (U4 is "exact").
  assert_strict_psu_folds(fold, strata, psu, V)
  fold
}

.check_m_h_design <- function(m_h, m_h_design) {
  if (!is.numeric(m_h_design) || anyNA(m_h_design))
    stop("m_h_design must be numeric without NA", call. = FALSE)
  if (length(m_h_design) == 1L && is.null(names(m_h_design))) {
    want <- rep(m_h_design, length(m_h))
    names(want) <- names(m_h)
  } else {
    if (is.null(names(m_h_design)) || !setequal(names(m_h_design), names(m_h)) ||
        anyDuplicated(names(m_h_design)))
      stop("m_h_design must be a scalar or be named by exactly the sampled strata",
           call. = FALSE)
    want <- m_h_design[names(m_h)]
  }
  bad <- m_h != want
  if (any(bad))
    stop("realised m_h differs from the design in ", sum(bad), " stratum/strata: ",
         .report_strata(m_h, bad), call. = FALSE)
  invisible(TRUE)
}

#' Independent U4 check on any unit-level fold vector (e.g. the folds an arm
#' actually used). Deliberately shares no code with the builder: grouping is by
#' sorted adjacent comparison, counting by table().
assert_strict_psu_folds <- function(fold, strata, psu, V) {
  V <- .check_V(V)
  n <- length(fold)
  if (!n || length(strata) != n || length(psu) != n)
    stop("fold, strata and psu must be non-empty and of equal length", call. = FALSE)
  f <- as.vector(fold)
  if (anyNA(f) || any(f != round(f)) || any(f < 1 | f > V))
    stop("every unit must carry exactly one fold label in 1..", V,
         " (whole-PSU folds)", call. = FALSE)
  o <- order(strata, psu, method = "radix")
  s <- strata[o]; p <- psu[o]; fo <- f[o]
  new_psu <- c(TRUE, s[-1L] != s[-n] | p[-1L] != p[-n])
  if (any(!new_psu & c(FALSE, fo[-1L] != fo[-n])))
    stop("a PSU is split across folds: folds must assign whole PSUs ",
         "(no PSU may be split)", call. = FALSE)
  counts <- table(as.character(s[new_psu]), factor(fo[new_psu], levels = seq_len(V)))
  m_h <- rowSums(counts)
  if (any(m_h %% V != 0L))
    stop("V = ", V, " does not divide m_h in stratum/strata: ",
         .report_strata(m_h, m_h %% V != 0L), " (exact balance is required)", call. = FALSE)
  unbalanced <- apply(counts, 1L, function(x) any(x != x[1L]))
  if (any(unbalanced))
    stop("fold counts are not exactly m_h/V per stratum in: ",
         .report_strata(m_h, unbalanced), " (exact balance is required)", call. = FALSE)
  invisible(TRUE)
}

.fold_V <- function(fold) {
  V <- attr(fold, "V")
  # "unit": individual-CF's unit-level folds (dev_unit_folds, R/estimators.R);
  # every PSU-level assertion still requires "strict_psu" (assert_strict_psu_folds).
  if (is.null(V) || !(identical(attr(fold, "fold_scheme"), "strict_psu") || identical(attr(fold, "fold_scheme"), "unit")))
    stop("expected a fold vector from make_strict_psu_folds() or dev_unit_folds()", call. = FALSE)
  V
}

#' Per-fold PSU sets: fold v scores `score` and its nuisance is trained on
#' `train` (the other folds' PSUs). Data frames of (stratum, psu).
strict_fold_psu_sets <- function(fold) {
  V <- .fold_V(fold)
  tab <- attr(fold, "psu_table")
  lapply(seq_len(V), function(v) list(
    score = tab[tab$fold == v, c("stratum", "psu")],
    train = tab[tab$fold != v, c("stratum", "psu")]))
}

#' PSUs, keyed by (stratum, psu), that occur in both row sets. Empty = disjoint.
psu_overlap <- function(strata, psu, fit_rows, score_rows) {
  idx <- .psu_index(strata, psu)
  both <- intersect(idx$unit_psu[fit_rows], idx$unit_psu[score_rows])
  idx$tab[sort(both), c("stratum", "psu"), drop = FALSE]
}

#' U5: for each fold, the PSUs whose units fitted the nuisance must be
#' disjoint from the PSUs it scores. PSU-CF only; for individual-CF the property
#' is false by design and for no-CF false by construction, so the check is
#' skipped there, never asserted.
#' @param fit_rows,score_rows lists (one element per fold) of unit row indices.
check_u5_disjoint <- function(strata, psu, fit_rows, score_rows,
                              cf_scheme = c("psu", "individual", "none")) {
  cf_scheme <- match.arg(cf_scheme)
  if (cf_scheme != "psu") return(invisible("skipped"))
  if (!is.list(fit_rows) || !is.list(score_rows) || length(fit_rows) != length(score_rows))
    stop("fit_rows and score_rows must be lists with one element per fold", call. = FALSE)
  for (v in seq_along(fit_rows)) {
    ov <- psu_overlap(strata, psu, fit_rows[[v]], score_rows[[v]])
    if (nrow(ov))
      stop("U5 violated in fold ", v, ": ", nrow(ov), " PSU(s) both fit and scored, e.g. ",
           paste0("(", utils::head(ov$stratum, 5L), ", ", utils::head(ov$psu, 5L), ")",
                  collapse = " "), " (fitting and scoring PSUs must be disjoint)", call. = FALSE)
  }
  invisible("pass")
}

#' Unit rows of each fold: the partition U3 needs for N_hat_v. U3 itself,
#' with the U1 design weights, is computed on its own code path elsewhere.
fold_rows <- function(fold) {
  V <- .fold_V(fold)
  unname(split(seq_along(fold), factor(as.vector(fold), levels = seq_len(V))))
}

#' Global relabel: new label = perm[old label]. A pure permutation of labels; the
#' PSU partition and every m_{v,h} are unchanged. Used for the U6 invariance limb,
#' which binds oracle and R-a only -- an arm-level test.
relabel_folds <- function(fold, perm) {
  V <- .fold_V(fold)
  if (!is.numeric(perm) || length(perm) != V || anyNA(perm) ||
      !setequal(perm, seq_len(V)) || anyDuplicated(perm))
    stop("perm must be a permutation of 1..", V, call. = FALSE)
  perm <- as.integer(perm)
  out <- perm[as.vector(fold)]
  attributes(out) <- attributes(fold)
  tab <- attr(fold, "psu_table")
  tab$fold <- perm[tab$fold]
  attr(out, "psu_table") <- tab
  prev <- attr(fold, "relabel")
  attr(out, "relabel") <- if (is.null(prev)) perm else perm[prev]
  out
}
