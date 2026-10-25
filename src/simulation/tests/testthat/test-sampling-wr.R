# Tests of the E2 WR code-check pieces:
# draw_sample_wr() (selection-keyed stage 2), check_u1_wr(), census_wr_target() and
# lin_var_wr(). The exact-enumeration test is the independent verification: under
# the WR design the census-centred linearised score with KNOWN N has design expectations
# E[total] = 0 and E[se2] = sigma2_m,WR = sum_h [(1 - 1/M_pop,h) B_h + V_h] EXACTLY, whereas the
# Hajek WR estimator does not (O(1/m) apart).

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R", "dgp.R",
            "census.R", "sampling.R", "variance.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

ROOT <- "SYNTHETIC-DEV-E2E"

# ---- exact enumeration on a tiny two-stage population (H = 2, J = 3, m_h = 2 WR, n0 = 2) ----
PSU <- data.frame(h = rep(1:2, each = 3), k = rep(1:3, 2), M = c(2L, 3L, 2L, 3L, 2L, 3L))
PHI <- c(0.3, -0.7, 1.9, 0.2, -1.1, 2.4, 0.8, 1.5, -0.4, 0.05, 3.1, -2.2, 0.6, 1.7, -0.9)
idx <- rep(seq_len(nrow(PSU)), PSU$M)
UNITS <- data.frame(h = PSU$h[idx], k = PSU$k[idx], i = sequence(PSU$M), phi0 = PHI, Q1 = PHI, Q0 = 0)
POP <- list(psu = PSU, units = UNITS, meta = list(label = "WR enumeration fixture"))
J <- 3L; M_H <- 2L; N0 <- 2L
N <- nrow(UNITS); PSI_N <- mean(UNITS$phi0)

# One stratum's WR outcomes: ordered selections (k1, k2) with probability (1/J)^2, each selection
# independently taking one of its choose(M_k, n0) subsets with equal probability.
stratum_outcomes <- function(h) {
  out <- list()
  for (k1 in seq_len(J)) for (k2 in seq_len(J)) {
    ks <- c(k1, k2)
    Ms <- PSU$M[PSU$h == h][ks]
    subs <- lapply(Ms, function(M) utils::combn(M, N0, simplify = FALSE))
    for (a in seq_along(subs[[1]])) for (b in seq_along(subs[[2]])) {
      rows <- do.call(rbind, lapply(1:2, function(j) {
        u <- UNITS[UNITS$h == h & UNITS$k == ks[j], ][subs[[j]][[c(a, b)[j]]], ]
        u$j <- j; u$M_k <- Ms[j]; u$w <- (J / M_H) * (Ms[j] / N0); u
      }))
      out[[length(out) + 1L]] <- list(p = (1 / J)^2 / (length(subs[[1]]) * length(subs[[2]])), rows = rows)
    }
  }
  out
}
O1 <- stratum_outcomes(1L); O2 <- stratum_outcomes(2L)
E <- local({
  acc <- c(p = 0, total = 0, se2 = 0, hajek = 0, Nhat = 0, Nhat2 = 0)
  for (a in O1) for (b in O2) {
    s <- rbind(a$rows, b$rows); p <- a$p * b$p
    u <- (s$phi0 - PSI_N) / N
    se2 <- .uc_var_of_totals(s$w * u, s$h, s$j)
    nh <- sum(s$w)
    acc <- acc + p * c(1, sum(s$w * u), se2, hajek_uc(s$phi0, s$w, s$h, s$j)$se2, nh, nh^2)
  }
  acc
})
CEN <- census_targets(POP, m_h = M_H, n0 = N0)
WR <- census_wr_target(CEN)

test_that("the WR enumeration is complete", {
  expect_equal(unname(E["p"]), 1, tolerance = 1e-14)
})

test_that("linearised known-N WR estimator: E[total] = 0 and E[se2] = sigma2_m,WR exactly (R14)", {
  expect_lt(abs(E[["total"]]), 1e-15)
  expect_equal(E[["se2"]], WR$sigma2_m_wr, tolerance = 1e-12)
  ps <- CEN$per_stratum
  expect_equal(WR$sigma2_m_wr, sum((1 - 1 / ps$M_pop_h) * ps$B_h + ps$V_h), tolerance = 1e-15)
  expect_equal(CEN$sigma2_m - WR$sigma2_m_wr, WR$wr_gap, tolerance = 1e-14)   # the WR gap
  expect_gt(WR$wr_gap, 0)
})

test_that("the Hajek WR estimator does not reproduce sigma2_m,WR exactly", {
  expect_gt(abs(E[["hajek"]] / WR$sigma2_m_wr - 1), 1e-6)
})

test_that("WR N_hat: E[N_hat] = N and Var(N_hat) = sum_h (M_h^2/m_h)(1 - 1/M_h) S2_M,h (U2, WR form)", {
  expect_equal(E[["Nhat"]], N, tolerance = 1e-12)
  S2M <- as.vector(tapply(PSU$M, PSU$h, stats::var)); Mh <- as.vector(table(PSU$h))
  expect_equal(E[["Nhat2"]] - E[["Nhat"]]^2, sum(Mh^2 / M_H * (1 - 1 / Mh) * S2M), tolerance = 1e-10)
})

# ---- the sampler on a synthetic DEV population ---------------------------------------------
CO <- dgp_coefficients(1.25)
BIG <- build_population(4, 4, CO, "WRTEST", 1, ROOT)
CELL <- list(regime = "many", H = 4, m_h = 2, J = 4, m = 8, f = "1/2", e2c = "wr")
key <- function(r) list(cell = CELL, population_id = 1, rep_id = r)

test_that("draw_sample_wr: U1 (WR) holds, keyed, scoped, and repeated selections draw independent subsamples", {
  s <- draw_sample_wr(BIG, 2, 20, key(1), ROOT)
  u1 <- check_u1_wr(s, BIG)
  expect_true(u1$pass)
  expect_identical(draw_sample_wr(BIG, 2, 20, key(1), ROOT), s)
  expect_false(identical(draw_sample_wr(BIG, 2, 20, key(2), ROOT)$i, s$i))
  withr::with_seed(9, {
    s0 <- get(".Random.seed", envir = globalenv())
    draw_sample_wr(BIG, 2, 20, key(3), ROOT)
    expect_identical(get(".Random.seed", envir = globalenv()), s0)
  })
  dups <- 0L; differ <- 0L
  for (r in 1:40) {
    x <- draw_sample_wr(BIG, 2, 20, key(r), ROOT)
    for (h in 1:4) {
      xs <- x[x$h == h, ]
      # Same PSU selected twice. Two independent n0-subsets of M_k units coincide with probability
      # 1/choose(M_k, n0), which is 1 at M_k = n0: compare only PSUs with M_k >= 25 (< 2e-5).
      if (length(unique(xs$k)) == 1L && xs$M_k[1] >= 25) {
        dups <- dups + 1L
        differ <- differ + !identical(sort(xs$i[xs$j == 1]), sort(xs$i[xs$j == 2]))
      }
    }
  }
  expect_gt(dups, 0L)                          # P(same PSU twice) = 1/4 per stratum
  expect_identical(differ, dups)               # selection keying: independent subsamples
  expect_identical(attr(s, "design")$wr, TRUE)
  expect_true(all(s$w == (4 / 2) * (s$M_k / 20)))
})

test_that("lin_var_wr: survey's WR svytotal equals the direct cluster-total formula", {
  s <- draw_sample_wr(BIG, 2, 20, key(5), ROOT)
  cen <- census_targets(population_prefix(BIG, 4, 4), m_h = 2, n0 = 20)
  lv <- lin_var_wr(s$phi0, cen$Psi_N, cen$N, s$w, s$h, s$j)
  expect_equal(lv$se2, lv$se2_direct, tolerance = 1e-12)
  expect_equal(lv$total, lv$total_direct, tolerance = 1e-12)
  expect_true(check_u8(s$phi0, s$w, s$h, s$j)$pass)    # Hajek UC on selections = survey, fpc not declared
})
