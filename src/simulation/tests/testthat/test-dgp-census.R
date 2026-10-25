# Tests of the population generator (R/dgp.R) and census targets (R/census.R).
# Every population here is a test population: the random-number conventions
# under a SYNTHETIC root, and a_cpsu = 1.25 is the
# study's value, passed explicitly. Checks:
# — prefix rule; law;
# targets and normings, U8; Delta_m identity;
# Supplement S1. Nothing here is a production population.

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R",
            "dgp.R", "census.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

ROOT <- "SYNTHETIC-DEV-E2E"
CO <- dgp_coefficients(1.25)   # explicit, never defaulted
MASTER <- "MSDEV"

# Clustered (by PSU) ratio mean and its MC standard error: units within a PSU share
# C_psu and the PSU effects, PSUs are independent under the generating law.
cl_mean_se <- function(x, psu_index) {
  t <- as.vector(rowsum(x, psu_index)); Mk <- as.vector(rowsum(rep(1, length(x)), psu_index))
  K <- length(t); xbar <- sum(t) / sum(Mk)
  c(mean = xbar, se = sqrt(K / (K - 1) * sum((t - xbar * Mk)^2)) / sum(Mk))
}
psu_index <- function(pop) rep(seq_len(nrow(pop$psu)), pop$psu$M)

# Large population shared by the law and identity tests (16,000 PSUs, ~544k units).
BIG <- build_population(800, 20, CO, MASTER, 2, ROOT)

test_that("guards: SYNTHETIC or the adopted production root only, a_cpsu explicit, outputs labelled", {
  # The production root is accepted by the guard; no production population is
  # built here (the guard is tested on its own), and any other root is refused.
  expect_error(build_population(2, 2, CO, MASTER, 1, root = "svytmle-eif-clt-sim-v9.7-rng2"), "SYNTHETIC")
  expect_error(build_population(2, 2, CO, MASTER, 1, root = "PRODUCTION"), "production root")
  expect_identical(assert_engine_root(SIM_PROD_ROOT), SIM_PROD_ROOT)
  expect_true(is_prod_root("svytmle-eif-clt-sim-v9.7-rng1"))
  expect_error(build_population(2, 2, list(g = CO$g), MASTER, 1, ROOT), "a_cpsu")
  expect_error(build_population(2, 2, root = ROOT, master_id = MASTER, pop_id = 1), "a_cpsu")
  p <- build_population(2, 3, CO, MASTER, 1, ROOT)
  expect_identical(p$meta[c("H", "J", "a_cpsu", "master_id", "pop_id", "root", "rng")],
                   list(H = 2L, J = 3L, a_cpsu = 1.25, master_id = MASTER, pop_id = 1, root = ROOT,
                        rng = "adopted RNG contract (docs/RNG-CONVENTIONS-DRAFT.md, R2)"))
  expect_match(p$meta$label, "spec v9.8 DGP")
  expect_match(census_targets(p, 2)$label, "^census targets")
  expect_named(p$psu, c("h", "k", "KB", "Cpsu", "M", "uA", "uY"))
  expect_named(p$units, c("h", "k", "i", "KB", "C1", "C2", "C3", "C4", "Cpsu", "A", "Y", "Q1", "Q0", "g", "phi0"))
  expect_error(census_targets(p, 1), "m_h")          # sigma2_m needs m_h >= 2
  expect_error(census_targets(p, 4), "m_h")          # m_h > M_pop,h
  expect_error(census_targets(p, 2, n0 = 49), "n0")  # pi_{i|k} > 1
})

test_that("caller's RNG state is untouched by a build", {
  set.seed(7); before <- .Random.seed; kind <- RNGkind()
  build_population(2, 2, CO, MASTER, 1, ROOT)
  expect_identical(.Random.seed, before); expect_identical(RNGkind(), kind)
})

test_that("prefix stability: (H, J) equals the prefix view of any larger build (:514-519, :532-533)", {
  master <- build_population(6, 7, CO, MASTER, 1, ROOT)
  for (hj in list(c(3, 4), c(6, 4), c(3, 7), c(1, 1), c(6, 7))) {
    direct <- build_population(hj[1], hj[2], CO, MASTER, 1, ROOT)
    view <- population_prefix(master, hj[1], hj[2])
    expect_identical(view$psu, direct$psu)
    expect_identical(view$units, direct$units)        # identity AND derived columns
    expect_identical(view$meta, direct$meta)
  }
  # A nested pair where neither is the 6 x 7 master, and the census on the view.
  a <- build_population(4, 5, CO, MASTER, 1, ROOT)
  expect_identical(population_prefix(a, 2, 3)$units, build_population(2, 3, CO, MASTER, 1, ROOT)$units)
  expect_identical(census_targets(population_prefix(master, 4, 5), 2), census_targets(a, 2))
  # The ids are load-bearing: another population or master is a different draw.
  expect_false(identical(build_population(3, 4, CO, MASTER, 2, ROOT)$psu$Cpsu, a$psu$Cpsu[a$psu$h <= 3 & a$psu$k <= 4]))
  expect_false(identical(build_population(3, 4, CO, "FSDEV", 1, ROOT)$psu$Cpsu, a$psu$Cpsu[a$psu$h <= 3 & a$psu$k <= 4]))
})

test_that("identity tier is coefficient-free; a_cpsu can flip only Y (RNG draft §6)", {
  p125 <- build_population(4, 5, CO, MASTER, 3, ROOT)
  p0 <- build_population(4, 5, dgp_coefficients(0), MASTER, 3, ROOT)
  id_cols <- c("h", "k", "i", "KB", "C1", "C2", "C3", "C4", "Cpsu", "A", "g")
  expect_identical(p0$units[id_cols], p125$units[id_cols])     # a_cpsu is not in g0
  expect_identical(p0$psu, p125$psu)
  expect_false(identical(p0$units$Q1, p125$units$Q1))
  expect_lt(mean(p0$units$Y != p125$units$Y), 0.25)            # thresholds: only near-threshold Y flip
})

test_that("law: K_B by parity, PSU sizes, PSU effects, covariates (§12.3 :1221-1243)", {
  ps <- BIG$psu; H <- BIG$meta$H
  expect_identical(ps$KB, ifelse(ps$h %% 2L == 1L, 1L, 2L))
  expect_identical(mean(ps$KB[ps$k == 1] == 1L), ceiling(H / 2) / H)
  expect_identical(kb_of_stratum(1:5), c(1L, 2L, 1L, 2L, 1L))  # odd H: class 1 holds ceil(H/2)
  # size law: M = rhe(34 + 14 C_psu) on {20..48}; exact discrete E M = 34, SD 8.093
  expect_identical(as.numeric(ps$M), round_half_even(34 + 14 * ps$Cpsu))
  expect_true(all(ps$M >= 20 & ps$M <= 48))
  cl <- size_law_cells(); p <- (cl$hi - cl$lo) / 2
  sd_law <- sqrt(sum(p * (cl$M - 34)^2))
  expect_equal(sum(p * cl$M), 34, tolerance = 1e-14)
  expect_equal(sd_law, 8.093, tolerance = 5e-4)
  n <- nrow(ps)
  expect_lt(abs(mean(ps$M) - 34), 4 * sd(ps$M) / sqrt(n))
  d2 <- (ps$M - mean(ps$M))^2
  expect_lt(abs(var(ps$M) - sd_law^2), 4 * sd(d2) / sqrt(n))
  expect_identical(BIG$meta$n_size_ties, 0L)
  # PSU effects by class: mean 0, SD sigma_{.,K_B}
  for (kb in 1:2) {
    for (v in c("uA", "uY")) {
      x <- ps[[v]][ps$KB == kb]; s <- if (v == "uA") CO$sigma_A[kb] else CO$sigma_Y[kb]
      expect_lt(abs(mean(x)), 4 * s / sqrt(length(x)))
      expect_lt(abs(var(x) - s^2), 4 * sd((x - mean(x))^2) / sqrt(length(x)))
    }
  }
  expect_lt(abs(cor(ps$uA, ps$uY)), 4 / sqrt(n))               # independent PSU effects
  expect_lt(abs(cor(ps$uA, ps$M)), 4 / sqrt(n))                # M_k free of u^A, u^Y
  u <- BIG$units
  for (v in c("C1", "C2", "C3", "C4")) {
    expect_true(all(u[[v]] > -1 & u[[v]] < 1))
    expect_lt(abs(mean(u[[v]])), 4 * sqrt(1 / 3 / nrow(u)))
    expect_lt(abs(var(u[[v]]) - 1 / 3), 4 * sd((u[[v]] - mean(u[[v]]))^2) / sqrt(nrow(u)))
  }
})

test_that("law: census means of A, Y, g, Q and Psi_N match the P0 grid within MC error", {
  ctx <- p0_context(BIG$meta$H, CO, second = list(rule = "clenshaw-curtis", n_c = 3L, n_psu = 2L))
  gp <- ctx$primary; tp <- ctx$primary_truth
  P0 <- c(g = sum(gp$w * tp$g), A = sum(gp$w * tp$g), Q1 = sum(gp$w * tp$Q1), Q0 = sum(gp$w * tp$Q0),
          # u^A and u^Y are independent, so P(Y=1|c) = g Q1 + (1-g) Q0 with marginal eta0
          Y = sum(gp$w * (tp$g * tp$Q1 + (1 - tp$g) * tp$Q0)),
          phi0 = sum(gp$w * (tp$Q1 - tp$Q0)),                   # Psi = P0 phi(eta0) (U9)
          KB1 = p0_class_shares(BIG$meta$H)[1])                 # unit share of class 1
  u <- BIG$units; ix <- psu_index(BIG)
  u$KB1 <- as.numeric(u$KB == 1L)
  for (v in names(P0)) {
    ms <- cl_mean_se(as.numeric(u[[v]]), ix)
    expect_lt(abs(ms[["mean"]] - P0[[v]]), 4 * ms[["se"]], label = paste("census mean of", v, "vs P0"))
  }
  ct <- census_targets(BIG, 2)
  expect_identical(ct$Psi_N, mean(u$phi0))
  ms <- cl_mean_se(u$Q1 - u$Q0, ix)
  expect_identical(ct$Psi_N_gf, mean(u$Q1 - u$Q0))
  expect_lt(abs(ct$Psi_N_gf - P0[["phi0"]]), 4 * ms[["se"]])
})

# Delta_m by the alternative identity N^-2 sum_h M_pop,h S2_t,h, computed
# here on a separate code path (tapply/var on per-PSU totals from aggregate()).
delta_alt <- function(pop, Psi_N) {
  z <- pop$units$phi0 - Psi_N
  tk <- aggregate(z, by = list(h = pop$units$h, k = pop$units$k), FUN = sum)
  S2 <- tapply(tk$x, tk$h, var); Mp <- tapply(tk$x, tk$h, length)
  sum(Mp * S2) / nrow(pop$units)^2
}

test_that("Delta_m = sum_h f_h B_h equals N^-2 sum_h M_pop,h S2_t,h to 1e-12, and is allocation-invariant", {
  cells <- list(list(pop = BIG, m_h = 2),                                   # f = 0.10
                list(pop = population_prefix(BIG, 50, 4), m_h = 2),         # strict host, f = 0.50
                list(pop = population_prefix(BIG, 50, 20), m_h = c(2, 6)[kb_of_stratum(1:50)]))
  for (cl in cells) {
    ct <- census_targets(cl$pop, cl$m_h)
    expect_lt(abs(ct$Delta_m - delta_alt(cl$pop, ct$Psi_N)) / ct$Delta_m, 1e-12)
    expect_lt(abs((ct$sigma2_m - ct$sigma2_m_wor) - ct$Delta_m) / ct$Delta_m, 1e-12)
    expect_true(ct$sigma2_m_wor < ct$sigma2_m)                              # dominance, Supplement S1
    expect_true(all(ct$per_stratum$beta_h >= 0 & ct$per_stratum$beta_h <= 1))
  }
  # the unequal-allocation pair: same population, allocations reversed
  p <- population_prefix(BIG, 50, 20)
  hi <- census_targets(p, c(2, 6)[kb_of_stratum(1:50)]); lo <- census_targets(p, c(6, 2)[kb_of_stratum(1:50)])
  expect_identical(hi$m, lo$m)
  expect_lt(abs(hi$Delta_m - lo$Delta_m) / hi$Delta_m, 1e-12)
  expect_false(isTRUE(all.equal(hi$sigma2_m, lo$sigma2_m)))                 # the contrast lives in sigma2_m
})

# Exact design moments by ENUMERATING every stage-1 SRS-WOR subset, with the stage-2
# variance from first/second-order inclusion probabilities (no S2 formula), and
# E[se^2 | s] = m/(m-1) {sum_j (y_j - ybar)^2 + (1 - 1/m) sum_j v_j}.
exact_moments <- function(pop, m_h, n0, Psi_N) {
  ps <- pop$psu; N <- nrow(pop$units)
  zs <- split(pop$units$phi0 - Psi_N, psu_index(pop))
  t_k <- vapply(zs, sum, 0)
  vt <- vapply(zs, function(x) {
    M <- length(x); p1 <- n0 / M; p2 <- n0 * (n0 - 1) / (M * (M - 1))
    (M / n0)^2 * (p1 * (1 - p1) * sum(x^2) + (p2 - p1^2) * (sum(x)^2 - sum(x^2)))
  }, 0)
  out <- NULL
  for (h in sort(unique(ps$h))) {
    r <- which(ps$h == h); Mh <- length(r); m <- m_h[h]; S <- combn(r, m)
    cm <- cv <- es <- numeric(ncol(S))
    for (j in seq_len(ncol(S))) {
      s <- S[, j]; y <- (Mh / m) * t_k[s] / N; v <- (Mh / m)^2 * vt[s] / N^2
      cm[j] <- sum(y); cv[j] <- sum(v)
      es[j] <- m / (m - 1) * (sum((y - mean(y))^2) + (1 - 1 / m) * sum(v))
    }
    out <- rbind(out, data.frame(h = h, var_wor = mean((cm - mean(cm))^2) + mean(cv), e_se2 = mean(es),
                                 mean_W = mean(cm)))
  }
  out
}

test_that("normings equal exact enumerated design moments (proofs.tex:95-107; U8 :346)", {
  p <- build_population(3, 6, CO, MASTER, 4, ROOT)
  for (m_h in list(c(2, 3, 4), c(2, 2, 2), c(5, 6, 3))) {
    ct <- census_targets(p, m_h)
    ex <- exact_moments(p, m_h, 20, ct$Psi_N)
    ps <- ct$per_stratum
    expect_lt(max(abs(ex$var_wor - ((1 - ps$f_h) * ps$B_h + ps$V_h)) / ex$var_wor), 1e-12)
    expect_lt(max(abs(ex$e_se2 - (ps$B_h + ps$V_h)) / ex$e_se2), 1e-12)
    expect_lt(abs(sum(ex$mean_W)), 1e-12 * sqrt(ct$sigma2_m_wor))          # L is design-unbiased for 0
  }
})

# Minimal two-stage SRS-WOR sampler, local to this test (not the sampling lane):
# stage 1 m_h of each stratum's PSUs, stage 2 n0 of each sampled PSU's units, by
# sample.int under a fixed local seed. Returns the linearised Hajek term
# L = N^-1 sum w z and the UC variance estimate se^2.
mc_two_stage <- function(pop, m_h, n0, Psi_N, R, seed) {
  ps <- pop$psu; N <- nrow(pop$units)
  zs <- split(pop$units$phi0 - Psi_N, psu_index(pop)); rows_h <- split(seq_len(nrow(ps)), ps$h)
  one <- function() {
    L <- 0; se2 <- 0
    for (h in seq_along(rows_h)) {
      r <- rows_h[[h]]; Mh <- length(r); m <- m_h[h]
      ks <- r[sample.int(Mh, m)]
      y <- vapply(ks, function(k) {
        Mk <- ps$M[k]; (Mh / m) * (Mk / n0) * sum(zs[[k]][sample.int(Mk, n0)]) / N
      }, 0)
      L <- L + sum(y); se2 <- se2 + m / (m - 1) * sum((y - mean(y))^2)
    }
    c(L = L, se2 = se2)
  }
  withr::with_seed(seed, t(replicate(R, one())), .rng_kind = "Mersenne-Twister",
                   .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
}

test_that("MC: Var(L) over 4000 two-stage draws matches sigma2_m_wor; E[se^2] matches sigma2_m", {
  p <- build_population(10, 8, CO, MASTER, 3, ROOT)
  m_h <- rep(4, 10)                                            # f_h = 0.5
  ct <- census_targets(p, m_h)
  R <- 4000
  X <- mc_two_stage(p, m_h, 20, ct$Psi_N, R, seed = 20261005)
  L <- X[, "L"]; s2 <- X[, "se2"]
  se_var <- sd((L - mean(L))^2) / sqrt(R)
  z_wor <- (var(L) - ct$sigma2_m_wor) / se_var
  z_m <- (mean(s2) - ct$sigma2_m) / (sd(s2) / sqrt(R))
  message(sprintf("MC R=%d: Var(L)/sigma2_m_wor = %.4f (z = %.2f); E[se2]/sigma2_m = %.4f (z = %.2f); Var(L)/sigma2_m = %.4f; Delta_m/sigma2_m = %.4f",
                  R, var(L) / ct$sigma2_m_wor, z_wor, mean(s2) / ct$sigma2_m, z_m,
                  var(L) / ct$sigma2_m, ct$Delta_m / ct$sigma2_m))
  expect_lt(abs(mean(L)), 4 * sd(L) / sqrt(R))
  expect_lt(abs(z_wor), 4)
  expect_lt(abs(z_m), 4)
  # Power: this census has Delta_m/sigma2_m well above the MC resolution, so the same
  # draws must reject sigma2_m as Var(L) (the FPC acts, U8's negative-limb logic).
  expect_gt(ct$Delta_m / ct$sigma2_m, 0.2)
  expect_gt(abs(var(L) - ct$sigma2_m) / se_var, 4)
})
