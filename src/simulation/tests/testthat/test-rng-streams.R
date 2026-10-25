# SYNTHETIC tests of the random-number key, hash and fill conventions: the
# prefix rule and seed levels
# fifth level and the size law.
# Every population below is a SYNTHETIC toy (placeholder variables,
# not the DGP); nothing here builds or inspects a production population.

source(test_path("..", "..", "R", "rng_streams.R"), local = TRUE)

ROOT <- "SYNTHETIC-TOY-rng-draft-1"

toy_cell <- function(...) {
  utils::modifyList(list(regime = "toy", H = 6, m_h = 2, J = 5, m = 12, f = "2/5",
                         e2c = "none"), list(...))
}

ctx_full <- function(...) {
  base <- list(root = ROOT, master_id = "TOYMASTER", population_id = 3, rep_id = 7,
               cell = toy_cell(), fold_structure_id = "psu-s1", arm = "toy-arm-1",
               fold = 1, statistic = "toy-stat", grading_pass = 1)
  utils::modifyList(base, list(...))
}

# Written out from the design, independently of the code's own table, so a
# drifted RNG_LEVEL_FIELDS cannot make the separation tests agree with itself.
SPEC_LEVEL_FIELDS <- list(
  population  = c("master_id", "population_id"),
  sample      = c("cell", "population_id", "rep_id"),
  fold        = c("cell", "population_id", "rep_id", "fold_structure_id"),
  learner     = c("cell", "population_id", "rep_id", "arm", "fold"),
  postprocess = c("cell", "statistic", "grading_pass"))

ENT <- list(population = list(type = "psu", h = 1, k = 1),
            sample = list(type = "stage1", h = 1),
            fold = list(type = "stratum", h = 1),
            learner = NULL, postprocess = NULL)

# TEST-ONLY negative control: one stream per population filled stratum-major --
# the shape of a direct-frame builder. Its fill offsets depend on J and on every
# earlier PSU's size, so it is expected to break the J-prefix (brief).
naive_population <- function(ctx, H, J) {
  st <- rng_state(rng_key("population", ctx))
  rows <- rng_scoped({
    rng_set_state(st)
    out <- list()
    for (h in seq_len(H)) for (k in seq_len(J)) out[[length(out) + 1L]] <- .toy_psu(h, k)
    out
  })
  pop <- as.data.frame(lapply(stats::setNames(nm = names(rows[[1]])), function(cl)
    unlist(lapply(rows, `[[`, cl), use.names = FALSE)))
  pop <- pop[order(pop$h, pop$k, pop$i), , drop = FALSE]
  rownames(pop) <- NULL
  toy_derive(pop)
}

test_that("SYNTHETIC (1): prefix views equal direct generation on every overlapping value", {
  ctx <- ctx_full()
  big <- toy_population(ctx, H = 9, J = 7)
  expect_true(length(unique(big$M)) > 3)  # variable sizes, so unit fill offsets vary
  for (s in list(c(4, 3), c(4, 7), c(9, 3), c(6, 5), c(1, 1), c(9, 7))) {
    direct <- toy_population(ctx, s[1], s[2])
    view <- toy_prefix_view(big, s[1], s[2])
    # IDs, covariates, exposure, outcome and the toy nuisances: every column.
    expect_identical(direct, view, info = paste(s, collapse = "x"))
    expect_identical(rng_store_hash(direct, TOY_IDENTITY_COLS, "SYNTHETIC-store"),
                     rng_store_hash(view, TOY_IDENTITY_COLS, "SYNTHETIC-store"))
  }
  # Nested non-maximal pair (H, J) = (4, 3) inside (H', J') = (6, 5).
  expect_identical(toy_prefix_view(toy_population(ctx, 6, 5), 4, 3),
                   toy_population(ctx, 4, 3))
})

test_that("SYNTHETIC (1): generation order cannot change any value", {
  ctx <- ctx_full()
  canon <- toy_population(ctx, 6, 5)
  n <- 6 * 5
  expect_identical(toy_population(ctx, 6, 5, traversal = rev(seq_len(n))), canon)
  expect_identical(toy_population(ctx, 6, 5, traversal = c(seq(2, n, 2), seq(1, n, 2))),
                   canon)
})

test_that("SYNTHETIC (1): negative control -- a single sequential fill breaks the J-prefix", {
  ctx <- ctx_full()
  nb <- naive_population(ctx, 6, 7)
  # Holding J, stratum-major fill happens to keep the H-prefix ...
  expect_identical(toy_prefix_view(nb, 4, 7), naive_population(ctx, 4, 7))
  # ... but changing J moves every later stratum's fill offset.
  expect_false(identical(toy_prefix_view(nb, 6, 5), naive_population(ctx, 6, 5)))
})

test_that("SYNTHETIC (2): each level ignores fields it must not depend on, reacts to those it must", {
  expect_identical(lapply(RNG_LEVEL_FIELDS, sort), lapply(SPEC_LEVEL_FIELDS, sort))
  expect_identical(RNG_CELL_FIELDS, c("regime", "H", "m_h", "J", "m", "f", "e2c"))
  base <- ctx_full()
  vary <- list(master_id = "TOYMASTER-B", population_id = 4, rep_id = 8,
               arm = "toy-arm-2", fold = 2, fold_structure_id = "psu-s2",
               statistic = "toy-stat-2", grading_pass = 2,
               cell = toy_cell(J = 10, f = "1/5"))
  for (lev in names(SPEC_LEVEL_FIELDS)) {
    u0 <- rng_uniforms(lev, base, 6, ENT[[lev]])
    for (fld in names(vary)) {
      ctx <- base
      ctx[[fld]] <- vary[[fld]]
      u1 <- rng_uniforms(lev, ctx, 6, ENT[[lev]])
      if (fld %in% SPEC_LEVEL_FIELDS[[lev]]) {
        expect_false(identical(u1, u0), info = paste(lev, fld))
      } else {
        expect_identical(u1, u0, info = paste(lev, fld))
      }
    }
  }
})

test_that("SYNTHETIC (2): cell never reaches the population stream; arm never the sample stream", {
  # Each cell field separately (f moves only with m_h/J).
  cells <- list(toy_cell(regime = "toy2"), toy_cell(H = 9, m = 18), toy_cell(m_h = 1, f = "1/5"),
                toy_cell(J = 10, f = "1/5"), toy_cell(m = 13), toy_cell(e2c = "toy-other"))
  p0 <- rng_uniforms("population", ctx_full(), 6, ENT$population)
  s0 <- rng_uniforms("sample", ctx_full(), 6, ENT$sample)
  for (cl in cells) {
    expect_identical(rng_uniforms("population", ctx_full(cell = cl), 6, ENT$population), p0)
    expect_false(identical(rng_uniforms("sample", ctx_full(cell = cl), 6, ENT$sample), s0))
  }
  # alpha lives inside arm: R-d rungs share sample and folds, differ in learner.
  ra <- ctx_full(arm = "Rd-AIPW-a030")
  rb <- ctx_full(arm = "Rd-AIPW-a035")
  expect_identical(rng_uniforms("sample", ra, 6, ENT$sample), rng_uniforms("sample", rb, 6, ENT$sample))
  expect_identical(rng_uniforms("fold", ra, 6, ENT$fold), rng_uniforms("fold", rb, 6, ENT$fold))
  expect_false(identical(rng_uniforms("learner", ra, 6), rng_uniforms("learner", rb, 6)))
  # Entity coordinates separate streams inside one level key.
  expect_false(identical(p0, rng_uniforms("population", ctx_full(), 6,
                                          list(type = "psu", h = 1, k = 2))))
})

test_that("SYNTHETIC (2): the level-2 toy sample is arm-free and prefix-consistent", {
  ctx <- ctx_full()
  big <- toy_population(ctx, 9, 8)
  s1 <- toy_sample(toy_population(ctx, 6, 5), ctx, n0 = 3)
  expect_identical(toy_sample(big, ctx, n0 = 3), s1)  # the master itself, frame = cell's prefix
  expect_identical(toy_sample(big, ctx_full(arm = "toy-arm-9", fold = 2,
                                            fold_structure_id = "unit-s3"), n0 = 3), s1)
  expect_false(identical(toy_sample(big, ctx_full(rep_id = 8), n0 = 3), s1))
  psu <- unique(s1[c("h", "k")])
  expect_identical(as.vector(table(psu$h)), rep(2L, 6))  # m_h per stratum
  expect_true(all(psu$k <= 5))                           # inside the first J PSUs
  expect_identical(anyDuplicated(s1), 0L)
})

test_that("SYNTHETIC (2): WR stage 2 is keyed by selection, so a repeated PSU gets its own subsample", {
  # (m_h = 2, J = 4: a stratum repeats a PSU with probability 1/J).
  # Toy coordinates only; docs.
  cell_wr <- toy_cell(H = 40, m_h = 2, J = 4, m = 80, f = "1/2", e2c = "toy-wr")
  ctx <- ctx_full(cell = cell_wr)
  pop <- toy_population(ctx, 40, 4)
  s <- toy_sample_wr(pop, ctx, n0 = 2)
  sel <- unique(s[c("h", "j", "k")])
  for (h in 1:40)  # stage 1: m_h with-replacement indices floor(J*U)+1 from stage1(h)
    expect_identical(sel$k[sel$h == h],
                     floor(4 * rng_uniforms("sample", ctx, 2, list(type = "stage1", h = h))) + 1)
  rep_h <- sel$h[duplicated(sel[c("h", "k")])]
  expect_true(length(rep_h) > 0)
  for (h in rep_h) {
    k <- sel$k[sel$h == h][1]
    Mk <- pop$M[pop$h == h & pop$k == k][1]
    u <- lapply(1:2, function(j) rng_uniforms("sample", ctx, Mk, list(type = "stage2wr", h = h, j = j)))
    expect_false(identical(u[[1]], u[[2]]), info = paste("stratum", h))
    for (j in 1:2)
      expect_identical(s$i[s$h == h & s$j == j],
                       sort(order(u[[j]], seq_len(Mk))[seq_len(min(2, Mk))]))
  }
  expect_identical(toy_sample_wr(pop, ctx_full(cell = cell_wr, arm = "toy-arm-9"), n0 = 2), s)
})

test_that("SYNTHETIC (2): canonical key bytes, and rejected encodings", {
  expect_identical(
    rng_key("population", ctx_full(), ENT$population),
    "SYNTHETIC-TOY-rng-draft-1;L=population;master_id=TOYMASTER;population_id=3;ent=psu(h=1,k=1)")
  expect_identical(
    rng_key("learner", ctx_full()),
    paste0("SYNTHETIC-TOY-rng-draft-1;L=learner;cell=(regime=toy,H=6,m_h=2,J=5,m=12,",
           "f=2/5,e2c=none);population_id=3;rep_id=7;arm=toy-arm-1;fold=1"))
  expect_identical(
    rng_encode_cell(list(regime = "toy", H = 4, m_h = c(K2 = 6, K1 = 2), J = 20, m = 16,
                         f = c(K1 = "1/10", K2 = "3/10"), e2c = "toy-unequal")),
    "(regime=toy,H=4,m_h=(K1=2,K2=6),J=20,m=16,f=(K1=1/10,K2=3/10),e2c=toy-unequal)")
  expect_identical(
    rng_key("sample", ctx_full(), list(type = "stage2wr", h = 3, j = 2)),
    paste0("SYNTHETIC-TOY-rng-draft-1;L=sample;cell=(regime=toy,H=6,m_h=2,J=5,m=12,",
           "f=2/5,e2c=none);population_id=3;rep_id=7;ent=stage2wr(h=3,j=2)"))
  expect_error(rng_key("sample", ctx_full(), list(type = "stage2wr", h = 1, k = 1)),
               "exactly coordinates")
  hk <- data.frame(h = c(1, 2, 3), k = c(4, 5, 6))
  expect_identical(rng_entity_keys("population", ctx_full(), "psu", hk),
                   vapply(1:3, function(r) rng_key("population", ctx_full(),
                     list(type = "psu", h = hk$h[r], k = hk$k[r])), ""))
  expect_error(rng_key("sample", ctx_full(cell = c(toy_cell(), V = 2))), "outside :565")
  expect_error(rng_key("sample", ctx_full(cell = c(toy_cell(), alpha = "0.30"))), "outside :565")
  expect_error(rng_key("sample", ctx_full(cell = toy_cell(f = 0.4))), "rational")
  expect_error(rng_key("sample", ctx_full(cell = toy_cell(f = "4/10"))), "reduced")
  expect_error(rng_key("sample", ctx_full(cell = toy_cell(f = "1/2"))), "not m_h/J")
  expect_error(rng_key("population", ctx_full(population_id = 3.5)), "whole number")
  expect_error(rng_key("learner", ctx_full(arm = "R-d(.30)")), "ASCII label")
  expect_error(rng_key("learner", ctx_full(arm = NULL)), "needs field")
  expect_error(rng_key("population", ctx_full(), list(type = "stage1", h = 1)), "not addressable")
  expect_error(rng_entity_keys("population", ctx_full(), "psu", data.frame(h = 1, k = 0.5)),
               "whole number")
  expect_error(toy_population(ctx_full(root = "PROD-ROOT"), 2, 2), "SYNTHETIC")
})

test_that("SYNTHETIC (3): the caller's RNG kind and state are untouched", {
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  saved_seed <- if (had) get(".Random.seed", envir = genv)
  saved_kind <- RNGkind()
  on.exit({
    suppressWarnings(RNGkind(saved_kind[1], saved_kind[2], saved_kind[3]))
    if (had) assign(".Random.seed", saved_seed, envir = genv)
  }, add = TRUE)

  kinds <- list(c("Mersenne-Twister", "Inversion", "Rejection"),
                c("Wichmann-Hill", "Box-Muller", "Rejection"),
                c("L'Ecuyer-CMRG", "Inversion", "Rejection"),
                c("Mersenne-Twister", "Inversion", "Rounding"))
  for (kd in kinds) {
    suppressWarnings(set.seed(20261004, kind = kd[1], normal.kind = kd[2], sample.kind = kd[3]))
    s0 <- get(".Random.seed", envir = genv)
    k0 <- RNGkind()
    invisible(toy_population(ctx_full(), 3, 3))
    invisible(toy_sample(toy_population(ctx_full(), 6, 5), ctx_full(), n0 = 2))
    invisible(toy_sample_wr(toy_population(ctx_full(), 6, 5), ctx_full(), n0 = 2))
    invisible(rng_uniforms("learner", ctx_full(), 4))
    expect_identical(get(".Random.seed", envir = genv), s0, info = kd[1])
    expect_identical(RNGkind(), k0, info = kd[1])
    next_draw <- stats::runif(2)
    suppressWarnings(set.seed(20261004, kind = kd[1], normal.kind = kd[2], sample.kind = kd[3]))
    expect_identical(stats::runif(2), next_draw, info = kd[1])
  }

  # Restored on error, too.
  set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  s0 <- get(".Random.seed", envir = genv)
  expect_error(rng_scoped({
    rng_set_state(rng_state("SYNTHETIC-x"))
    stats::runif(3)
    stop("boom")
  }), "boom")
  expect_identical(get(".Random.seed", envir = genv), s0)

  # A caller with no .Random.seed still has none afterwards, and keeps its kind.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  rm(".Random.seed", envir = genv)
  invisible(rng_uniforms("sample", ctx_full(), 4, ENT$sample))
  expect_false(exists(".Random.seed", envir = genv, inherits = FALSE))
  expect_identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection"))

  # A stream can be positioned only inside a scope.
  expect_error(rng_set_state(rng_state("SYNTHETIC-x")), "inside rng_scoped")
})

test_that("SYNTHETIC (4): rounding rule on the four exact ties of the §12.3 size law", {
  # -- M_k = 34 + 14*C_psu, rounded. C_psu = +-0.25, +-0.75 give exact .5
  # values in binary64 as well as in real arithmetic.
  cdrv <- c(0.25, -0.25, 0.75, -0.75)
  x <- 34 + 14 * cdrv
  expect_identical(x, c(37.5, 30.5, 44.5, 23.5))
  expect_identical(size_from_driver(cdrv, 34, 14), c(38, 30, 44, 24))
  # Half-to-even keeps the reflection M(-c) = 68 - M(c) at these ties.
  expect_identical(size_from_driver(cdrv, 34, 14) + size_from_driver(-cdrv, 34, 14), rep(68, 4))
  # Base round() (IEC 60559) agrees; half-up, the named alternative, differs at 30.5, 44.5.
  expect_identical(round_half_even(x), round(x))
  expect_identical(floor(x + 0.5), c(38, 31, 45, 24))
  expect_identical(round_half_even(c(37.49999999999999, 37.50000000000001, 20, 48, 0.5, 1.5, 2.5)),
                   c(37, 38, 20, 48, 0, 2, 2))
  expect_error(round_half_even(-0.5), "x >= 0")
})

test_that("SYNTHETIC (4): ties reachable on R's L'Ecuyer runif lattice are produced by runif itself", {
  # u = d * normc, C = -1 + 2u. Exactly these four d give a binary64 tie in
  # 34 + 14*C (enumerated arithmetically; docs). Only the first is a real-
  # arithmetic tie; for the other three the binary64 sum rounds onto the tie.
  d <- c(536870886, 1610612658, 2684354430, 3758096202)
  cl <- -1 + 2 * (d * .MRG_NORM)
  expect_identical(cl[1], -0.75)
  expect_true(all(cl[2:4] != c(-0.25, 0.25, 0.75)))
  expect_identical(34 + 14 * cl, c(23.5, 30.5, 37.5, 44.5))
  expect_identical(size_from_driver(cl, 34, 14), c(24, 30, 38, 44))

  # Build a state whose first output is d (component 1 emits p1 = d, component 2
  # emits p2 = 0) and draw it with runif() under the scoped stream.
  mulmod <- function(a, b, m) {
    bh <- b %/% 65536
    bl <- b %% 65536
    ((((a * bh) %% m) * 65536) %% m + (a * bl) %% m) %% m
  }
  powmod <- function(a, e, m) {
    r <- 1
    a <- a %% m
    while (e > 0) {
      if (e %% 2 == 1) r <- mulmod(r, a, m)
      a <- mulmod(a, a, m)
      e <- e %/% 2
    }
    r
  }
  inv_a12 <- powmod(1403580, .MRG_M1 - 2, .MRG_M1)  # m1 is prime
  expect_identical(mulmod(inv_a12, 1403580, .MRG_M1), 1)
  for (j in seq_along(d)) {
    st <- c(0, mulmod(d[j], inv_a12, .MRG_M1), 0, 0, 1, 0)
    u <- rng_scoped({
      rng_set_state(st)
      stats::runif(1)
    })
    expect_identical(u, rng_mrg32k3a_emulate(st, 1))
    expect_identical(size_from_driver(-1 + 2 * u, 34, 14), c(24, 30, 38, 44)[j])
  }
})

test_that("SYNTHETIC (5): repeated calls are bit-identical", {
  expect_identical(toy_population(ctx_full(), 5, 4), toy_population(ctx_full(), 5, 4))
  for (lev in names(RNG_LEVEL_FIELDS))
    expect_identical(rng_uniforms(lev, ctx_full(), 20, ENT[[lev]]),
                     rng_uniforms(lev, ctx_full(), 20, ENT[[lev]]))
  big <- toy_population(ctx_full(), 6, 5)
  expect_identical(toy_sample(big, ctx_full(), n0 = 3), toy_sample(big, ctx_full(), n0 = 3))
})

test_that("SYNTHETIC (5): known answers pin hash, state, stream and store bytes", {
  # FIPS 180-4 SHA-256("abc"): the hash is the standard function on raw bytes.
  expect_identical(digest::digest("abc", algo = "sha256", serialize = FALSE),
                   "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")
  key <- rng_key("population", ctx_full(), ENT$population)
  st <- rng_state(key)
  # Recorded on Windows R 4.5.1 and reproduced by an independent Python
  # implementation (hashlib + integer MRG32k3a). On another platform or R
  # version this block is the identity-layer equality probe (docs).
  expect_identical(st, c(1368465174, 3927175565, 1022830719, 3752302388, 3478437990, 775587058))
  u <- rng_scoped({
    rng_set_state(st)
    stats::runif(50)
  })
  expect_identical(sprintf("%a", u[1:3]),
                   c("0x1.f23c7ba4d1248p-3", "0x1.dd380f3dbd8c7p-1", "0x1.d15f59ee1d792p-1"))
  expect_identical(u, rng_mrg32k3a_emulate(st, 50))
  pop <- toy_population(ctx_full(), 5, 4)
  expect_identical(rng_store_hash(pop, TOY_IDENTITY_COLS, "SYNTHETIC-store"),
                   "8ff4849e0b31fa0c68be3932841b2b0b89fde695ca4cc4c4f9f690becec5a0b9")
})

test_that("SYNTHETIC: states R would reject are re-derived, never handed to R", {
  bad <- rbind(c(0, 0, 0, 1, 2, 3, 0, 0),             # all-zero first triple
               c(1, 2, 3, 0, 0, 0, 0, 0),             # all-zero second triple
               c(2^31, 1, 1, 1, 1, 1, 0, 0),          # no non-NA int form
               c(.MRG_M1, 2 * .MRG_M1, 0, 5, 5, 5, 0, 0))  # reduces to zeros
  expect_identical(rng_state_from_words(bad)$valid, rep(FALSE, 4))
  ok <- rng_state_from_words(c(.MRG_M1 + 5, 1, 1, .MRG_M2 + 7, 1, 1, 0, 0))
  expect_true(ok$valid)
  expect_identical(ok$state[c(1, 4)], c(5, 7))
})

test_that("SYNTHETIC: store identity is coefficient-free; derived layer and row order behave", {
  ctx <- ctx_full()
  p1 <- toy_population(ctx, 4, 4)
  alt <- list(a = c(-0.3, 0.1, 0.2, 0.9), y = c(0.4, -0.8, 0.1, 0.2))
  p2 <- toy_population(ctx, 4, 4, coef = alt)
  expect_identical(p1[TOY_IDENTITY_COLS], p2[TOY_IDENTITY_COLS])
  expect_identical(rng_store_hash(p1, TOY_IDENTITY_COLS, "SYNTHETIC-store"),
                   rng_store_hash(p2, TOY_IDENTITY_COLS, "SYNTHETIC-store"))
  expect_false(identical(rng_store_hash(p1, TOY_DERIVED_COLS, "SYNTHETIC-store"),
                         rng_store_hash(p2, TOY_DERIVED_COLS, "SYNTHETIC-store")))
  # Inversion coupling: the new exposure thresholds the same stored uniform.
  expect_identical(p2$a, as.integer(p1$va < p2$pa))
  # Canonical row order, not input order, defines the bytes.
  expect_identical(rng_store_hash(p1[rev(seq_len(nrow(p1))), ], TOY_IDENTITY_COLS, "SYNTHETIC-store"),
                   rng_store_hash(p1, TOY_IDENTITY_COLS, "SYNTHETIC-store"))
  expect_false(identical(
    rng_store_hash(toy_population(ctx_full(population_id = 4), 4, 4), TOY_IDENTITY_COLS, "SYNTHETIC-store"),
    rng_store_hash(p1, TOY_IDENTITY_COLS, "SYNTHETIC-store")))
  bad <- p1
  bad$x1[1] <- NA
  expect_error(rng_store_bytes(bad, TOY_IDENTITY_COLS, "SYNTHETIC-store"), "NA/NaN")
})
