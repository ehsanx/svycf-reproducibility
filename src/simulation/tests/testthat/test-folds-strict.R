# Strict whole-PSU folds: U4/U5, G3(c), level 3.

source(test_path("..", "..", "R", "folds.R"), local = TRUE)

# TEST-ONLY fold-stream stand-ins. The production keying and traversal
# come from R/estimators.R and R/rng_streams.R; these neither propose nor
# fix it. One hashes each (key, stratum, psu) independently (order-free, touches
# no RNG); the other reseeds the global RNG, to exercise the builder's restore.
draw_hash_testonly <- structure(function(fold_key, psus) {
  base <- digest::digest(fold_key, algo = "sha256")
  vapply(seq_len(nrow(psus)), function(i) {
    h <- digest::digest(list(base, psus$stratum[i], psus$psu[i]), algo = "sha256")
    strtoi(substr(h, 1L, 6L), 16L) * 2^24 + strtoi(substr(h, 7L, 12L), 16L)
  }, numeric(1))
}, convention = "TEST-ONLY sha256 per PSU")

draw_reseed_testonly <- function(kind = "Mersenne-Twister") {
  function(fold_key, psus) {
    s <- strtoi(substr(digest::digest(fold_key, algo = "sha256"), 1L, 7L), 16L)
    set.seed(s, kind = kind)
    stats::runif(nrow(psus))
  }
}

key <- function(cell = "test-cell", pop = 1L, rep = 1L, fs = 1L)
  list(cell = cell, population_id = pop, rep_id = rep, fold_structure_id = fs)

# Synthetic sample: m_h PSUs per stratum, variable units per PSU, rows shuffled so
# nothing relies on PSU-contiguous rows. psu ids repeat across strata by default,
# as they may in a real design.
make_sample <- function(m_h, units = c(1L, 3L), repeat_ids = TRUE, seed = 1L) {
  set.seed(seed)
  s <- rep(seq_along(m_h), m_h)
  p <- if (repeat_ids) sequence(m_h) else seq_along(s)
  k <- units[1L] - 1L + sample.int(units[2L] - units[1L] + 1L, length(s), replace = TRUE)
  d <- data.frame(stratum = rep(s, k), psu = rep(p, k))
  d[sample(nrow(d)), , drop = FALSE]
}

# PSU-level fold counts by stratum, computed from the unit rows alone.
psu_counts <- function(d, fold) {
  u <- unique(data.frame(stratum = d$stratum, psu = d$psu, fold = as.vector(fold)))
  expect_equal(anyDuplicated(u[, c("stratum", "psu")]), 0L)  # whole PSUs
  table(u$stratum, u$fold)
}

# Every cell shape whose arms use PSU-CF, with its V:
# V=2 in every cell except the E3 sensitivity
# limb V=5. Fixed-strata ladder from including the
# cheap-arm (320) and R-d-only [640] rungs; many-strata from including (800).
spec_cells <- c(
  lapply(c(25L, 50L, 100L, 200L, 400L, 800L), function(H) list(
    id = paste0("E1 many-strata H=", H), m_h = rep(2L, H), V = 2L)),
  lapply(c(10L, 20L, 40L, 80L, 160L, 320L, 640L), function(mh) list(
    id = paste0("E1/R-d fixed-strata m_h=", mh), m_h = rep(mh, 5L), V = 2L)),
  list(
    list(id = "E2 strict fitted rung H=100 (:1142)", m_h = rep(2L, 100L), V = 2L),
    list(id = "E2 FPC sweep fitted f=0.5 H=50 (:1146)", m_h = rep(2L, 50L), V = 2L),
    list(id = "E3 R-c foil / R-a / power panel m_h=80 (:1144, :1148)",
         m_h = rep(80L, 5L), V = 2L),
    list(id = "E3 V-sensitivity m_h=40 V=2 (:1145)", m_h = rep(40L, 5L), V = 2L),
    list(id = "E3 V-sensitivity m_h=40 V=5 (:1145, :1152)", m_h = rep(40L, 5L), V = 5L)))

test_that("every PSU-CF spec cell gets exactly m_h/V whole PSUs per fold", {
  for (cell in spec_cells) {
    d <- make_sample(cell$m_h, units = c(1L, 2L))
    f <- make_strict_psu_folds(d$stratum, d$psu, cell$V, key(cell$id), draw_hash_testonly)
    expect_identical(attr(f, "V"), cell$V, info = cell$id)
    expect_true(all(attr(f, "m_h") == cell$m_h), info = cell$id)
    expect_true(all(attr(f, "m_vh") == cell$m_h / cell$V), info = cell$id)
    tc <- psu_counts(d, f)
    expect_equal(ncol(tc), cell$V, info = cell$id)
    expect_true(all(tc == cell$m_h[1L] / cell$V), info = cell$id)
    expect_silent(assert_strict_psu_folds(f, d$stratum, d$psu, cell$V))
  }
})

test_that("E3 V-sensitivity gives m_vh = 20 at V=2 and 8 at V=5 (:1152)", {
  d <- make_sample(rep(40L, 5L))
  f2 <- make_strict_psu_folds(d$stratum, d$psu, 2, key(), draw_hash_testonly)
  f5 <- make_strict_psu_folds(d$stratum, d$psu, 5, key(), draw_hash_testonly)
  expect_true(all(attr(f2, "m_vh") == 20L))
  expect_true(all(attr(f5, "m_vh") == 8L))
})

test_that("unequal-allocation pair 6/2 and its reversal are integral at V=2 (:1147, :1170)", {
  # Oracle-only cell; G3(c) integrality is still asserted for it at V=2.
  # K_B = 1 on odd h, 2 on even h; K_B = 2 is the designated high class.
  hi_on_hi <- ifelse(seq_len(50L) %% 2L == 0L, 6L, 2L)
  for (m_h in list(hi_on_hi, 8L - hi_on_hi)) {
    d <- make_sample(m_h)
    f <- make_strict_psu_folds(d$stratum, d$psu, 2, key(), draw_hash_testonly,
                               m_h_design = setNames(m_h, seq_along(m_h)))
    expect_equal(sum(m_h), 200L)
    expect_true(all(attr(f, "m_vh") == m_h / 2L))
    tc <- psu_counts(d, f)
    expect_true(all(tc[, 1L] == m_h / 2L & tc[, 2L] == m_h / 2L))
  }
})

test_that("generic arithmetic outside the spec grid: m_h=80, V=5 gives 16 per fold", {
  # Not a study cell: every m_h=80 cell runs at V=2.
  d <- make_sample(rep(80L, 5L))
  f <- make_strict_psu_folds(d$stratum, d$psu, 5, key(), draw_hash_testonly)
  expect_true(all(psu_counts(d, f) == 16L))
})

test_that("incompatible requests are hard errors naming the strata, never a reduced V", {
  d <- make_sample(rep(3L, 4L))
  expect_error(make_strict_psu_folds(d$stratum, d$psu, 2, key(), draw_hash_testonly),
               "V = 2 does not divide m_h in 4 stratum.*'1' \\(m_h = 3\\)")
  # Mixed m_h: only the failing strata are named.
  d <- make_sample(c(2L, 3L, 4L, 5L, 6L))
  err <- tryCatch(make_strict_psu_folds(d$stratum, d$psu, 2, key(), draw_hash_testonly),
                  error = conditionMessage)
  expect_match(err, "in 2 stratum/strata: '2' \\(m_h = 3\\), '4' \\(m_h = 5\\)\\.")
  expect_false(grepl("'1'|'3'|'5'", sub("\\. V \\|.*", "", err)))
  expect_match(err, "never reduced")
  # The legacy builder would shrink these to V_eff = 2; here they stop.
  d <- make_sample(rep(2L, 25L))
  expect_error(make_strict_psu_folds(d$stratum, d$psu, 5, key(), draw_hash_testonly),
               "V = 5 does not divide m_h in 25 stratum.*15 more")
  d <- make_sample(rep(40L, 5L))
  expect_error(make_strict_psu_folds(d$stratum, d$psu, 3, key(), draw_hash_testonly),
               "V = 3 does not divide")
  for (bad_V in list(1, 2.5, NA_real_, c(2, 5), "2", Inf))
    expect_error(make_strict_psu_folds(d$stratum, d$psu, bad_V, key(), draw_hash_testonly),
                 "V must be a single whole number >= 2")
})

test_that("realised m_h is checked against the design when supplied", {
  d <- make_sample(c(80L, 80L, 78L, 80L, 80L))  # V=2 divides 78: only the design check catches it
  expect_error(make_strict_psu_folds(d$stratum, d$psu, 2, key(), draw_hash_testonly,
                                     m_h_design = 80),
               "differs from the design in 1 stratum.*'3' \\(m_h = 78\\)")
  expect_error(make_strict_psu_folds(d$stratum, d$psu, 2, key(), draw_hash_testonly,
                                     m_h_design = c(a = 80)),
               "named by exactly the sampled strata")
})

test_that("malformed inputs stop", {
  d <- make_sample(rep(2L, 4L))
  expect_error(make_strict_psu_folds(d$stratum[-1], d$psu, 2, key(), draw_hash_testonly),
               "equal length")
  s <- d$stratum; s[3] <- NA
  expect_error(make_strict_psu_folds(s, d$psu, 2, key(), draw_hash_testonly), "NA stratum")
  expect_error(make_strict_psu_folds(d$stratum, d$psu, 2, key(), "not a function"),
               "draw must be a function")
  expect_error(make_strict_psu_folds(d$stratum, d$psu, 2, key(),
                                     function(k, p) stats::runif(nrow(p) - 1L)),
               "finite numeric scores")
  expect_error(make_strict_psu_folds(d$stratum, d$psu, 2, key(),
                                     function(k, p) rep(NA_real_, nrow(p))),
               "finite numeric scores")
})

test_that("fold_key carries exactly the level-3 factors (:560-567)", {
  d <- make_sample(rep(2L, 4L))
  mk <- function(k) make_strict_psu_folds(d$stratum, d$psu, 2, k, draw_hash_testonly)
  expect_error(mk(key()[-4]), "exactly the fields")
  expect_error(mk(c(key(), arm = "R-b")), "exactly the fields")
  expect_error(mk(key(cell = list(regime = "fixed", m_h = 40, V = 5))), "must not contain v")
  expect_error(mk(key(cell = list(regime = "fixed", alpha = 0.25))), "must not contain alpha")
  expect_error(mk(key(rep = 1:2)), "rep_id must be a single id")
  expect_error(mk(key(pop = NA)), "population_id is empty or contains NA")
  ok <- mk(key(cell = list(regime = "fixed", H = 5, m_h = 40, J = 400, m = 200, f = 0.1)))
  expect_identical(attr(ok, "fold_key")$cell$m_h, 40)
})

test_that("whole PSUs: every unit of a PSU carries its PSU's single fold", {
  d <- make_sample(rep(10L, 5L), units = c(1L, 7L))
  f <- make_strict_psu_folds(d$stratum, d$psu, 2, key(), draw_hash_testonly)
  per_psu <- tapply(as.vector(f), list(d$stratum, d$psu), function(x) length(unique(x)))
  expect_true(all(per_psu == 1L))
  expect_true(all(as.vector(f) %in% 1:2))
  expect_type(as.vector(f), "integer")
  tab <- attr(f, "psu_table")
  expect_identical(as.vector(f),
                   tab$fold[match(paste(d$stratum, d$psu), paste(tab$stratum, tab$psu))])
})

test_that("PSU ids repeating across strata are distinct PSUs", {
  d <- make_sample(rep(2L, 60L), repeat_ids = TRUE)  # psu ids 1,2 in every stratum
  f <- make_strict_psu_folds(d$stratum, d$psu, 2, key(), draw_hash_testonly)
  expect_true(all(attr(f, "m_h") == 2L))
  expect_equal(nrow(attr(f, "psu_table")), 120L)
  tab <- attr(f, "psu_table")
  # psu id 1 lands in fold 1 in some strata and fold 2 in others (P(all equal) = 2^-59).
  expect_length(unique(tab$fold[tab$psu == 1L]), 2L)
  # Same partition whether ids repeat or are globally unique (same canonical order).
  d2 <- d; d2$psu <- (d2$stratum - 1L) * 2L + d2$psu
  f2 <- make_strict_psu_folds(d2$stratum, d2$psu, 2,
                              key(), function(k, p) draw_hash_testonly(k, transform(p, psu = (p$psu - 1L) %% 2L + 1L)))
  expect_identical(as.vector(f2), as.vector(f))
  # Character ids, also repeating across strata.
  dc <- data.frame(stratum = paste0("h", d$stratum), psu = paste0("k", d$psu))
  fc <- make_strict_psu_folds(dc$stratum, dc$psu, 2, key(), draw_hash_testonly)
  expect_silent(assert_strict_psu_folds(fc, dc$stratum, dc$psu, 2))
})

test_that("determinism: same key -> identical folds; any key factor change -> new folds", {
  d <- make_sample(rep(2L, 100L))
  mk <- function(k, drw = draw_hash_testonly)
    make_strict_psu_folds(d$stratum, d$psu, 2, k, drw)
  f <- mk(key())
  expect_identical(mk(key()), f)
  expect_identical(attr(f, "draw_convention"), "TEST-ONLY sha256 per PSU")
  # Each level-3 factor moves the partition (P(no change) = 2^-100 for each).
  for (k in list(key(fs = 2L), key(rep = 2L), key(pop = 2L), key(cell = "other")))
    expect_false(identical(as.vector(mk(k)), as.vector(f)))
  # The reseeding stand-in is deterministic in its key as well.
  rs <- draw_reseed_testonly()
  expect_identical(mk(key(), rs), mk(key(), rs))
  expect_false(identical(as.vector(mk(key(fs = 2L), rs)), as.vector(mk(key(), rs))))
})

test_that("folds do not depend on ambient RNG state or on sample row order", {
  d <- make_sample(rep(10L, 5L))
  mk <- function(dd, drw = draw_hash_testonly)
    make_strict_psu_folds(dd$stratum, dd$psu, 2, key(), drw)
  for (drw in list(draw_hash_testonly, draw_reseed_testonly())) {
    set.seed(1); a <- mk(d, drw)
    set.seed(99); b <- mk(d, drw)
    rm(".Random.seed", envir = globalenv()); c0 <- mk(d, drw)
    expect_identical(a, b); expect_identical(a, c0)
    perm <- sample(nrow(d))
    p <- mk(d[perm, ], drw)
    expect_identical(attr(p, "psu_table"), attr(a, "psu_table"))
    expect_identical(as.vector(p), as.vector(a)[perm])
  }
})

test_that("a draw that reads the ambient RNG is rejected", {
  d <- make_sample(rep(2L, 10L))
  expect_error(make_strict_psu_folds(d$stratum, d$psu, 2, key(),
                                     function(k, p) stats::runif(nrow(p))),
               "depends on the ambient RNG state")
})

test_that("the caller's RNG state is untouched", {
  d <- make_sample(rep(2L, 20L))
  mk <- function(drw) make_strict_psu_folds(d$stratum, d$psu, 2, key(), drw)
  # Initialised state, including a draw that switches RNG kind.
  for (drw in list(draw_hash_testonly, draw_reseed_testonly(),
                   draw_reseed_testonly("L'Ecuyer-CMRG"))) {
    set.seed(2026); s0 <- .Random.seed
    mk(drw)
    expect_identical(.Random.seed, s0)
    expect_identical(RNGkind()[1L], "Mersenne-Twister")
    x <- stats::runif(3); set.seed(2026); expect_identical(stats::runif(3), x)
  }
  # Uninitialised state stays uninitialised, with the original kind.
  rm(".Random.seed", envir = globalenv())
  mk(draw_reseed_testonly("L'Ecuyer-CMRG"))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_identical(RNGkind()[1L], "Mersenne-Twister")
  # Restored on error too.
  set.seed(7); s0 <- .Random.seed
  expect_error(mk(function(k, p) { set.seed(1); stop("boom") }), "boom")
  expect_identical(.Random.seed, s0)
  # Also with the purity check off (the draw then sees the caller's state once).
  set.seed(8); s0 <- .Random.seed
  make_strict_psu_folds(d$stratum, d$psu, 2, key(), draw_reseed_testonly(),
                        check_draw_purity = FALSE)
  expect_identical(.Random.seed, s0)
})

test_that("U5: each fold's training PSUs are disjoint from the PSUs it scores", {
  for (V in c(2L, 5L)) {
    d <- make_sample(rep(40L, 5L))
    f <- make_strict_psu_folds(d$stratum, d$psu, V, key(), draw_hash_testonly)
    sets <- strict_fold_psu_sets(f)
    expect_length(sets, V)
    k <- function(x) paste(x$stratum, x$psu)
    all_psu <- k(attr(f, "psu_table"))
    for (v in seq_len(V)) {
      expect_length(intersect(k(sets[[v]]$train), k(sets[[v]]$score)), 0L)
      expect_setequal(c(k(sets[[v]]$train), k(sets[[v]]$score)), all_psu)
    }
    # Scoring sets partition the sampled PSUs: each PSU scored exactly once.
    scored <- unlist(lapply(sets, function(s) k(s$score)))
    expect_setequal(scored, all_psu); expect_equal(anyDuplicated(scored), 0L)
    # Unit-row form, as an arm would supply it.
    fit <- lapply(seq_len(V), function(v) which(f != v))
    sco <- lapply(seq_len(V), function(v) which(f == v))
    expect_identical(check_u5_disjoint(d$stratum, d$psu, fit, sco), "pass")
  }
})

test_that("U5 check can fail, keys on (stratum, psu), and is skipped off PSU-CF", {
  d <- make_sample(rep(2L, 30L), repeat_ids = TRUE)
  f <- make_strict_psu_folds(d$stratum, d$psu, 2, key(), draw_hash_testonly)
  fit <- list(which(f != 1L), which(f != 2L))
  sco <- list(which(f == 1L), which(f == 2L))
  # Leak one unit of a scored PSU into fold 1's training rows.
  leaky <- fit; leaky[[1]] <- c(leaky[[1]], sco[[1]][1])
  expect_error(check_u5_disjoint(d$stratum, d$psu, leaky, sco),
               "U5 violated in fold 1: 1 PSU")
  expect_equal(nrow(psu_overlap(d$stratum, d$psu, leaky[[1]], sco[[1]])), 1L)
  # Same psu id in a different stratum is a different PSU: not a leak.
  a <- which(d$stratum == 1L & d$psu == 1L)
  b <- which(d$stratum == 2L & d$psu == 1L)
  expect_identical(check_u5_disjoint(d$stratum, d$psu, list(b), list(a)), "pass")
  expect_equal(nrow(psu_overlap(d$stratum, d$psu, b, a)), 0L)
  # individual-CF: false by design; no-CF: false by construction. Skipped, never asserted.
  expect_identical(check_u5_disjoint(d$stratum, d$psu, leaky, sco, "individual"), "skipped")
  everything <- list(seq_len(nrow(d)))
  expect_identical(check_u5_disjoint(d$stratum, d$psu, everything, everything, "none"),
                   "skipped")
})

test_that("fold_rows gives the unit partition U3 needs, with no weights invented", {
  d <- make_sample(rep(40L, 5L), units = c(1L, 5L))
  f <- make_strict_psu_folds(d$stratum, d$psu, 5, key(), draw_hash_testonly)
  fr <- fold_rows(f)
  expect_length(fr, 5L)
  expect_identical(sort(unlist(fr)), seq_len(nrow(d)))
  for (v in 1:5) expect_identical(fr[[v]], which(f == v))
  expect_error(fold_rows(1:3), "expected a fold vector")
})

test_that("global relabel is a pure permutation of labels", {
  d <- make_sample(rep(40L, 5L))
  for (V in c(2L, 5L)) {
    f <- make_strict_psu_folds(d$stratum, d$psu, V, key(), draw_hash_testonly)
    perm <- if (V == 2L) c(2L, 1L) else c(3L, 5L, 1L, 2L, 4L)
    g <- relabel_folds(f, perm)
    expect_identical(as.vector(g), perm[as.vector(f)])
    expect_identical(attr(g, "m_vh"), attr(f, "m_vh"))
    expect_silent(assert_strict_psu_folds(g, d$stratum, d$psu, V))
    # Same PSU partition: new fold perm[v] scores exactly old fold v's PSUs.
    sf <- strict_fold_psu_sets(f); sg <- strict_fold_psu_sets(g)
    for (v in seq_len(V))
      expect_identical(sg[[perm[v]]]$score, sf[[v]]$score)
    # Inverse relabel restores the labels exactly.
    back <- relabel_folds(g, order(perm))
    expect_identical(as.vector(back), as.vector(f))
    expect_identical(attr(back, "psu_table"), attr(f, "psu_table"))
    expect_identical(attr(back, "relabel"), seq_len(V))
  }
  expect_error(relabel_folds(f, c(1L, 1L, 2L, 3L, 4L)), "permutation of 1..5")
  expect_error(relabel_folds(f, 1:4), "permutation of 1..5")
})

test_that("the independent U4 checker rejects approximate or split folds", {
  # Legacy-style approximate balance (rep_len over m_h = 3 at V = 2).
  d <- make_sample(rep(3L, 4L), units = c(1L, 1L))
  approx <- ave(d$psu, d$stratum, FUN = function(x) rep_len(1:2, length(x))[rank(x)])
  expect_error(assert_strict_psu_folds(approx, d$stratum, d$psu, 2), "does not divide")
  # V divides m_h but the split is 3-1.
  d <- make_sample(rep(4L, 3L), units = c(1L, 1L))
  lop <- ifelse(d$psu <= 3L, 1L, 2L)
  expect_error(assert_strict_psu_folds(lop, d$stratum, d$psu, 2), "not exactly m_h/V")
  # A PSU split across folds.
  d <- make_sample(rep(2L, 3L), units = c(2L, 2L))
  good <- make_strict_psu_folds(d$stratum, d$psu, 2, key(), draw_hash_testonly)
  split_psu <- as.vector(good); i <- which(d$stratum == 1L & d$psu == 1L)[1]
  split_psu[i] <- 3L - split_psu[i]
  expect_error(assert_strict_psu_folds(split_psu, d$stratum, d$psu, 2), "split across folds")
  bad_label <- as.vector(good); bad_label[1] <- 3L
  expect_error(assert_strict_psu_folds(bad_label, d$stratum, d$psu, 2), "fold label in 1..2")
})
