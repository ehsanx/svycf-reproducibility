# DEV / PROVISIONAL tests of the two-stage sampler (R/sampling.R) and the UC
# variance (R/variance.R). SYNTHETIC fixtures only: the populations below are built
# here in the dgp lane's shape (pop$psu, pop$units) and are not the DGP
# (no A, Y or eta0; phi0 is a placeholder score with a PSU effect). Sample draws
# use the level-2 streams under a
# SYNTHETIC root. U1, U8 and the UC variance estimator of the paper.

for (f in c("rng_streams.R", "sampling.R", "variance.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

ROOT <- "SYNTHETIC-DEV-E2E"

# Sizes on the frozen law M_k = round(34 + 14 C_psu) in {20..48}; K_B by
# stratum parity; phi0 carries a PSU effect and a size-linked mean, so both
# between- and within-PSU variance are present. Local seed, restored on exit.
make_pop <- function(H, J, seed = 1L, size = c(34, 14)) {
  withr::with_seed(seed, {
    psu <- data.frame(h = rep(seq_len(H), each = J), k = rep(seq_len(J), times = H))
    psu$KB <- ifelse(psu$h %% 2L == 1L, 1L, 2L)
    psu$Cpsu <- stats::runif(nrow(psu), -1, 1)
    psu$M <- round_half_even(size[1] + size[2] * psu$Cpsu)
    u <- stats::rnorm(nrow(psu), 0, c(0.5, 1)[psu$KB])
    idx <- rep(seq_len(nrow(psu)), psu$M)
    units <- data.frame(h = psu$h[idx], k = psu$k[idx], i = sequence(psu$M))
    units$C1 <- stats::runif(nrow(units), -1, 1)
    units$phi0 <- 0.2 + 0.6 * units$C1 + u[idx] + 0.5 * psu$Cpsu[idx] + stats::rnorm(nrow(units))
    list(psu = psu, units = units)
  }, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
}

frac <- function(p, q) { g <- .gcd(p, q); paste0(p / g, "/", q / g) }

cell_for <- function(H, J, m_h) {
  if (length(m_h) == 1L) {
    return(list(regime = "dev", H = H, m_h = m_h, J = J, m = H * m_h, f = frac(m_h, J), e2c = "none"))
  }
  n_cls <- c(K1 = ceiling(H / 2), K2 = floor(H / 2))          # parity allocation,
  list(regime = "dev", H = H, m_h = m_h, J = J, m = sum(n_cls * m_h[names(n_cls)]),
       f = vapply(m_h, frac, "", q = J), e2c = "none")
}

key_for <- function(cell, population_id = 1, rep_id = 1, ...)
  list(cell = cell, population_id = population_id, rep_id = rep_id, ...)

pop48 <- make_pop(4, 8)
cell48 <- cell_for(4, 8, 3)

test_that("U1: weights are (M_pop,h/m_h)(M_k/n0) bit-for-bit, checked against the recorded design", {
  s <- draw_sample(pop48, 3, 20, key_for(cell48))
  expect_identical(names(s)[1:8], c("h", "k", "i", "w", "pi_psu", "pi_unit", "M_pop_h", "M_k"))
  expect_true(all(c("C1", "phi0", "KB", "Cpsu") %in% names(s)))
  u1 <- check_u1(s, pop48)
  expect_true(u1$pass, info = paste(u1$checks$check[!u1$checks$pass], collapse = "; "))
  expect_identical(s$w, (8 / 3) * (s$M_k / 20))
  expect_true(all(s$M_pop_h == 8) && all(s$pi_psu == 3 / 8))
  expect_true(all(s$pi_unit >= 20 / 48 & s$pi_unit <= 1))     # pi in [0.417, 1]
  expect_match(attr(s, "design")$label, "adopted RNG contract")
  # negative control: a weight off by 1e-12 is not the recorded design
  bad <- s
  bad$w[5] <- bad$w[5] + 1e-12
  expect_false(check_u1(bad, pop48)$pass)
  bad <- s
  bad$pi_unit[1] <- 1
  expect_false(check_u1(bad, pop48)$pass)
})

test_that("every stratum has exactly m_h PSUs and every PSU exactly n0 units", {
  pop <- make_pop(6, 8, seed = 2L)
  for (alloc in list(3, c(K1 = 4, K2 = 2))) {
    cell <- cell_for(6, 8, alloc)
    for (r in 1:5) {
      s <- draw_sample(pop, alloc, 20, key_for(cell, rep_id = r))
      kb <- ifelse(seq_len(6) %% 2L == 1L, "K1", "K2")
      want <- if (length(alloc) == 1L) rep(alloc, 6) else unname(alloc[kb])
      psus <- unique(s[, c("h", "k")])
      expect_identical(as.integer(table(factor(psus$h, levels = 1:6))), as.integer(want))
      expect_true(all(table(paste(s$h, s$k)) == 20))
      expect_true(all(psus$k >= 1 & psus$k <= 8))
      expect_false(anyDuplicated(s[, c("h", "k", "i")]) > 0)
      expect_identical(unname(attr(s, "design")$m_h), as.integer(want))
      expect_true(check_u1(s, pop)$pass)
    }
  }
  # the allocation drawn must be the one the sample key names
  expect_error(draw_sample(pop, 2, 20, key_for(cell_for(6, 8, 3))), "key\\$cell\\$m_h")
})

test_that("the frame is the cell's prefix of a larger master", {
  master <- make_pop(7, 12, seed = 3L)
  view <- list(psu = master$psu[master$psu$h <= 4 & master$psu$k <= 8, ],
               units = master$units[master$units$h <= 4 & master$units$k <= 8, ])
  view$units <- view$units[nrow(view$units):1, ]                 # row order must not matter
  for (r in 1:3) {
    a <- draw_sample(master, 3, 20, key_for(cell48, rep_id = r))
    b <- draw_sample(view, 3, 20, key_for(cell48, rep_id = r))
    expect_identical(a, b)
  }
  expect_error(draw_sample(view, 3, 20, key_for(cell_for(5, 8, 3))), "prefix frame")
})

test_that("M_k < n0 is a hard error; allow_take_all takes the whole PSU", {
  small <- make_pop(3, 6, size = c(6, 4))                       # M_k in 2..10 < n0
  cell <- cell_for(3, 6, 2)
  expect_error(draw_sample(small, 2, 20, key_for(cell)), "infeasible")
  s <- draw_sample(small, 2, 20, key_for(cell), allow_take_all = TRUE)
  expect_true(all(s$pi_unit == 1))
  expect_identical(as.numeric(table(paste(s$h, s$k))),
                   as.numeric(tapply(s$M_k, paste(s$h, s$k), `[`, 1)))
  expect_true(check_u1(s, small)$pass)
  expect_identical(attr(s, "design")$n_take_all, 6L)
})

test_that("deterministic given the key, different for different keys, and arm-free", {
  k1 <- key_for(cell48, population_id = 1, rep_id = 1)
  s1 <- draw_sample(pop48, 3, 20, k1)
  expect_identical(draw_sample(pop48, 3, 20, k1), s1)
  expect_identical(draw_sample(sampling_frame(pop48, 4, 8), 3, 20, k1), s1)
  # arm, fold, alpha never reach the sample stream
  expect_identical(draw_sample(pop48, 3, 20, key_for(cell48, arm = "R-a", fold = 2, alpha = 0.25,
                                                      fold_structure_id = "psu-s1")), s1)
  ids <- function(s) paste(s$h, s$k, s$i)
  expect_false(identical(ids(draw_sample(pop48, 3, 20, key_for(cell48, rep_id = 2))), ids(s1)))
  expect_false(identical(ids(draw_sample(pop48, 3, 20, key_for(cell48, population_id = 2))), ids(s1)))
  expect_false(identical(ids(draw_sample(pop48, 3, 20, k1, root = "SYNTHETIC-OTHER")), ids(s1)))
  # a different cell (here: another E2-construction label) is another stream
  c2 <- cell48; c2$e2c <- "other"
  expect_false(identical(ids(draw_sample(pop48, 3, 20, key_for(c2))), ids(s1)))
  # 20 reps give 20 distinct samples
  many <- vapply(1:20, function(r) paste(ids(draw_sample(pop48, 3, 20, key_for(cell48, rep_id = r))),
                                         collapse = ","), "")
  expect_equal(length(unique(many)), 20L)
  expect_error(draw_sample(pop48, 3, 20, k1, root = "PRODUCTION"), "SYNTHETIC")
  expect_error(draw_sample(pop48, 3, 20, c(k1, root = "SYNTHETIC-X")), "differs")
})

test_that("draws follow the PROPOSED level-2 rule exactly as the prototype toy_sample()", {
  ctx <- c(key_for(cell48, population_id = 4, rep_id = 9), root = ROOT)
  df <- pop48$units
  df$M <- pop48$psu$M[match(paste(df$h, df$k), paste(pop48$psu$h, pop48$psu$k))]
  toy <- toy_sample(df, ctx, n0 = 20)
  s <- draw_sample(pop48, 3, 20, key_for(cell48, population_id = 4, rep_id = 9))
  expect_identical(as.numeric(s$h), as.numeric(toy$h))
  expect_identical(as.numeric(s$k), as.numeric(toy$k))
  expect_identical(as.numeric(s$i), as.numeric(toy$i))
})

test_that("the caller's RNG is untouched", {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  old_kind <- RNGkind()
  withr::defer({
    suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
    if (is.null(old)) { if (exists(".Random.seed", globalenv())) rm(".Random.seed", envir = globalenv()) }
    else assign(".Random.seed", old, envir = globalenv())
  })
  frame <- sampling_frame(pop48, 4, 8)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
  set.seed(20261005)
  before <- get(".Random.seed", globalenv())
  kind <- RNGkind()
  s <- draw_sample(frame, 3, 20, key_for(cell48))
  expect_identical(get(".Random.seed", globalenv()), before)
  expect_identical(RNGkind(), kind)
  rm(".Random.seed", envir = globalenv())
  s2 <- draw_sample(frame, 3, 20, key_for(cell48))
  expect_false(exists(".Random.seed", globalenv()))
  expect_identical(s2, s)                    # and the ambient state never feeds the draw
})

test_that("U8: hajek_uc reproduces survey::svymean (fpc not declared) to <= 1e-12 relative", {
  pop <- make_pop(6, 8, seed = 5L)
  worst <- 0
  for (alloc in list(3, c(K1 = 4, K2 = 2))) {
    cell <- cell_for(6, 8, alloc)
    for (r in 1:4) {
      s <- draw_sample(pop, alloc, 20, key_for(cell, rep_id = r))
      for (phi in list(s$phi0, s$phi0^2, as.numeric(s$phi0 > 0.3), 50 + s$C1)) {
        u8 <- check_u8(phi, s$w, s$h, s$k)
        expect_true(u8$pass, info = sprintf("rel %.3g psi %.3g", u8$rel_diff, u8$psi_rel_diff))
        worst <- max(worst, u8$rel_diff)
        # PSU ids unique across strata or repeating (k): same PSUs, same answer (TC-V1)
        alt <- hajek_uc(phi, s$w, s$h, 1000 * s$h + s$k)
        ref <- hajek_uc(phi, s$w, s$h, s$k)
        expect_equal(alt$se2, ref$se2, tolerance = 1e-14)
        expect_identical(alt$psi, ref$psi)
      }
    }
  }
  expect_lt(worst, 1e-12)
})

test_that("hajek_uc: eq:ucvar by explicit loops, N_hat is sum(w), shift invariance, m_h >= 2", {
  d <- data.frame(h = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2), k = c(1, 1, 2, 2, 3, 3, 1, 1, 2, 2),
                  w = c(2, 2, 3, 3, 1.5, 1.5, 4, 4, 2.5, 2.5),
                  phi = c(0.3, 1.1, -0.4, 0.8, 2.0, 0.1, 0.5, 0.7, -1.2, 0.9))
  r <- hajek_uc(d$phi, d$w, d$h, d$k)
  Nh <- sum(d$w); ps <- sum(d$w * d$phi) / Nh
  se2 <- 0
  for (h in unique(d$h)) {
    ks <- unique(d$k[d$h == h]); mh <- length(ks)
    t <- vapply(ks, function(k) { j <- d$h == h & d$k == k; sum(d$w[j] * (d$phi[j] - ps)) / Nh }, 0)
    se2 <- se2 + mh / (mh - 1) * sum((t - mean(t))^2)
  }
  expect_equal(r$se2, se2, tolerance = 1e-14)
  expect_identical(r$N_hat, Nh)
  expect_equal(r$psi, ps, tolerance = 1e-15)
  expect_identical(unname(r$m_h), c(3L, 2L))
  expect_equal(hajek_uc(d$phi + 7, d$w, d$h, d$k)$se2, r$se2, tolerance = 1e-12)
  expect_error(hajek_uc(d$phi[1:8], d$w[1:8], d$h[1:8], d$k[1:8]), ">= 2 sampled PSUs")
  expect_error(hajek_uc(d$phi, -d$w, d$h, d$k), "w must be")
})

test_that("U8 two-stage limb helper runs (dry-run object; not gated here)", {
  s <- draw_sample(pop48, 3, 20, key_for(cell48, rep_id = 3))
  r <- se2_two_stage_wor(s$phi0, s$w, s$h, s$k, s$i, s$M_pop_h, s$M_k)
  expect_true(is.finite(r$se2_wor) && r$se2_wor >= 0)
  expect_equal(r$psi, hajek_uc(s$phi0, s$w, s$h, s$k)$psi, tolerance = 1e-12)
  expect_match(r$label, "dry-run object")
})

test_that("HT total, N_hat and the Hajek mean are design-unbiased; inclusion rates are pi_k pi_{i|k}", {
  pop <- make_pop(6, 10, seed = 7L)
  cell <- cell_for(6, 10, 4)
  frame <- sampling_frame(pop, 6, 10)
  R <- 1000L
  T <- sum(frame$units$phi0); N <- nrow(frame$units); Psi_N <- T / N
  est <- matrix(NA_real_, R, 3)
  hits <- numeric(N)
  for (r in seq_len(R)) {
    s <- draw_sample(frame, 4, 20, key_for(cell, population_id = 1, rep_id = r))
    est[r, ] <- c(sum(s$w * s$phi0), sum(s$w), hajek_uc(s$phi0, s$w, s$h, s$k)$psi)
    row <- frame$offset[(s$h - 1) * 10 + s$k] + s$i
    hits[row] <- hits[row] + 1
  }
  z <- (colMeans(est) - c(T, N, Psi_N)) / (apply(est, 2, stats::sd) / sqrt(R))
  expect_true(all(abs(z) < 4), info = paste(format(z, digits = 3), collapse = " "))
  p <- (4 / 10) * (20 / frame$units$M_k)
  zi <- (hits - R * p) / sqrt(R * p * (1 - p))
  expect_lt(max(abs(zi)), 5)
  psu_hits <- tapply(hits, (frame$units$h - 1) * 10 + frame$units$k, sum) / 20
  zp <- (psu_hits - R * 0.4) / sqrt(R * 0.4 * 0.6)
  expect_lt(max(abs(zp)), 4.5)
})
