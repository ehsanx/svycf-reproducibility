# Tests of the E3 power-control panel pieces: the
# wrong-clustering se2 variants and their census predictions, U10, the engineered aligned-error
# construction (C-AE-1..4), the C1 shift Z and G6 limb (a)'s boundary distance. The census
# predictions are verified by EXACT ENUMERATION on a tiny two-stage population.

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R", "folds.R",
            "dgp.R", "census.R", "sampling.R", "variance.R", "estimators.R", "controls.R", "learner_rd.R",
            "gates.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

ROOT <- "SYNTHETIC-DEV-E2E"

# ---- exact enumeration: H = 2, J = 3, m_h = 2 (SRS-WOR), n0 = 2 -----------------------------
PSU <- data.frame(h = rep(1:2, each = 3), k = rep(1:3, 2), M = c(2L, 3L, 3L, 3L, 2L, 3L))
PHI <- c(0.3, -0.7, 1.9, 0.2, -1.1, 2.4, 0.8, 1.5, -0.4, 0.05, 3.1, -2.2, 0.6, 1.7, -0.9, 0.35)
idx <- rep(seq_len(nrow(PSU)), PSU$M)
UNITS <- data.frame(h = PSU$h[idx], k = PSU$k[idx], i = sequence(PSU$M), phi0 = PHI)
POP <- list(psu = PSU, units = UNITS, meta = list(label = "controls enumeration fixture"))
M_H <- 2L; N0 <- 2L; J <- 3L
N <- nrow(UNITS); PSI_N <- mean(UNITS$phi0)

E <- local({
  pairs <- utils::combn(J, M_H, simplify = FALSE)
  acc <- c(p = 0, uc = 0, nomh = 0, unit = 0, nostrat = 0, u10_psi = 0, u10_se2 = 0)
  for (a in pairs) for (b in pairs) {
    sel <- rbind(PSU[PSU$h == 1 & PSU$k %in% a, ], PSU[PSU$h == 2 & PSU$k %in% b, ])
    subs <- lapply(sel$M, function(Mk) utils::combn(Mk, N0, simplify = FALSE))
    grid <- expand.grid(lapply(subs, seq_along))
    for (g in seq_len(nrow(grid))) {
      parts <- lapply(seq_len(nrow(sel)), function(j) {
        u <- UNITS[UNITS$h == sel$h[j] & UNITS$k == sel$k[j], ][subs[[j]][[grid[g, j]]], ]
        u$M_k <- sel$M[j]; u$M_pop_h <- J; u$w <- (J / M_H) * (sel$M[j] / N0); u
      })
      s <- do.call(rbind, parts)
      attr(s, "design") <- list(m_h = c(M_H, M_H))
      p <- (1 / choose(J, M_H))^2 * prod(1 / choose(sel$M, N0))
      wc <- wrong_clustering_from_z(s$w * (s$phi0 - PSI_N) / N, s$h, s$k)
      ue <- u10_estimate(s)
      acc <- acc + p * c(1, wc[["uc"]], wc[["nomh"]], wc[["unit"]], wc[["nostrat"]], ue$psi, ue$se2)
    }
  }
  acc
})
WCC <- wrong_clustering_census(POP, m_h = M_H, n0 = N0)
U10 <- u10_census(POP, m_h = M_H, n0 = N0)

test_that("enumeration complete", expect_equal(unname(E["p"]), 1, tolerance = 1e-14))

test_that("wrong-clustering census predictions are EXACT for the linearised known-N score", {
  expect_equal(E[["uc"]], WCC$expect[["uc"]], tolerance = 1e-12)
  expect_equal(E[["nomh"]], WCC$expect[["nomh"]], tolerance = 1e-12)
  expect_equal(E[["unit"]], WCC$expect[["unit"]], tolerance = 1e-12)
  expect_equal(E[["nostrat"]], WCC$expect[["nostrat"]], tolerance = 1e-12)
  expect_equal(WCC$R_var[["nomh"]], (M_H - 1) / M_H, tolerance = 1e-12)   # 1/2 here; 79/80 in the panel
})

test_that("U10: E[psi_U10] = Psi_U10 (mean of PSU means) and E[se2_U10] = sigma2_m,U10 exactly", {
  expect_equal(E[["u10_psi"]], U10$Psi_U10, tolerance = 1e-12)
  expect_equal(E[["u10_se2"]], U10$sigma2_m_u10, tolerance = 1e-12)
  phibar <- tapply(UNITS$phi0, paste(UNITS$h, UNITS$k), mean)
  expect_equal(U10$Psi_U10, mean(phibar), tolerance = 1e-14)
  expect_equal(U10$N_hat1, N0 * nrow(PSU))
  expect_gt(abs(U10$bias), 1e-3)                                        # sizes informative here
})

# ---- a synthetic DEV population for the sampled pieces --------------------------------------
CO <- dgp_coefficients(1.25)
BIG <- build_population(5, 40, CO, "CTLTEST", 1, ROOT)
CELL <- list(regime = "fixed", H = 5, m_h = 4, J = 40, m = 20, f = "1/10", e2c = "none")
KEY <- list(cell = CELL, population_id = 1, rep_id = 1)
S <- draw_sample(BIG, 4, 20, KEY, ROOT)

test_that("wrong-clustering: survey's wrong designs reproduce the direct variants; uc = hajek_uc", {
  wc <- wrong_clustering_se2(S$phi0, S$w, S$h, S$k)
  sv <- wrong_clustering_survey(S$phi0, S$w, S$h, S$k)
  expect_equal(wc[["uc"]], hajek_uc(S$phi0, S$w, S$h, S$k)$se2, tolerance = 1e-12)
  for (v in c("uc", "nomh", "unit", "nostrat")) expect_equal(wc[[v]], sv[[v]], tolerance = 1e-10, info = v)
  expect_equal(wc[["nomh"]], wc[["uc"]] * 3 / 4, tolerance = 1e-12)
})

test_that("U10 estimator: stage-1 weights only, N_hat1 = n0 M_pop-sample constant, keyed to the oracle sample", {
  ue <- u10_estimate(S)
  expect_equal(ue$N_hat, 20 * nrow(BIG$psu), tolerance = 1e-9)        # sum w1 = n0 * M_pop exactly in algebra
  expect_true(all(ue$w1 == 40 / 4))
  expect_true(is.finite(ue$psi) && ue$se2 > 0)
})

# ---- aligned-error construction (C-AE-1..4) --------------------------------------------------
CTX <- p0_context(5L, CO, primary = list(rule = "gauss-legendre", n_c = 6L, n_psu = 2L),
                  second = list(rule = "clenshaw-curtis", n_c = 7L, n_psu = 3L))

test_that("AE: the closed-form map equals the generic T2 code, is positive, increasing and quadratic near 0", {
  for (k in c(0.05, 0.2)) {
    fit <- ae_fit(k, CO)
    expect_equal(ae_t2_map(k, CTX$primary_truth, CTX$primary$w), t2_closed_form(CTX, list(fit), 1)$T2, tolerance = 1e-12)
  }
  t2 <- vapply(c(0.01, 0.02, 0.1, 0.2), ae_t2_map, numeric(1), truth = CTX$primary_truth, w = CTX$primary$w)
  expect_true(all(t2 > 0) && all(diff(t2) > 0))
  expect_equal(t2[2] / t2[1], 4, tolerance = 0.02)                      # ~kappa^2 near 0
  expect_lt(abs(ae_t2_map(0, CTX$primary_truth, CTX$primary$w)), 1e-15)
})

test_that("AE: U9 brute force agrees with the closed form; the arm passes U3/U7 with no clip hits", {
  fit <- ae_fit(0.15, CO)
  xc <- u9_cross_check(CTX, list(fit))
  expect_true(all(xc$pass))
  fold <- dev_psu_folds(S, KEY, 2, ROOT)
  cf <- psu_cf_aipw(S, fold, function(train) ae_fit(0.15, CO), arm = ARM_AE_AIPW)
  or <- oracle_aipw(S)
  expect_true(all(check_u3(S$w, cf$fold, cf$omega, cf$N_hat_v)$pass))
  expect_true(check_u7(S, cf$fold, cf$fits, cf$psi, or$psi)$pass)
  expect_identical(cf$clip_hits, 0L)
  expect_equal(cf$g, ae_shift(S$g, 0.15), tolerance = 1e-15)
})

test_that("AE: clip bound and kappa solve (target reached; not buildable when out of reach)", {
  cb <- ae_clip_bound(CO)
  expect_true(cb$inside); expect_gt(cb$kappa_max, 0)
  sw <- c(0.012, 0.014)
  ks <- ae_kappa_solve(CTX$primary_truth, CTX$primary$w, sw, cb$kappa_max)
  expect_true(ks$buildable && ks$monotone)
  expect_equal(ks$stat, PANEL_KAPPA, tolerance = 1e-9)
  expect_equal(ks$stat, unname(stats::quantile(ks$T2 / sw, 0.9, type = 7)), tolerance = 1e-12)
  expect_lt(ks$kappa, cb$kappa_max)
  nb <- ae_kappa_solve(CTX$primary_truth, CTX$primary$w, c(10, 20), cb$kappa_max)   # sigma so large 0.30 is out of reach
  expect_false(nb$buildable); expect_match(nb$reason, "within the clip bound")
})

test_that("C1 shift and g-formula Z; G6 limb (a) boundary reading (R5)", {
  z <- panel_shift_z(psi = 0.2, se = 0.01, Psi_N = 0.2, Psi_N_gf = 0.19, sigma_wor = 0.01)
  expect_equal(z[["Z_shift"]], -0.30, tolerance = 1e-12)
  expect_equal(z[["Z_gf"]], 1, tolerance = 1e-12)
  g <- g6_boundary(c(0.9875, 1.4, -0.29, 0.30), c(0.85, 0.85, -0.10, 0), c(1.15, 1.15, 0.10, 0.10), c(NA, NA, 0.05, 0.1))
  expect_identical(g$inside, c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(g$distance, c(-0.1375, 0.25, 0.19, 0.20), tolerance = 1e-12)
  expect_match(g$status[1], "inside the compliant region")
  expect_match(g$status[2], "no MCSE given")
  expect_match(g$status[3], "graded power control")
  expect_match(g$status[4], "at most 3 MCSE")
})
