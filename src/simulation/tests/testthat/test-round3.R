# Unit tests of the production machinery:
# run modes and identities (R/run_mode.R), the grading primitives (R/grading.R), the frozen kappa_AE record
# (R/controls.R), the production-population gate statistics (R/census.R), the store hashes (R/dgp.R), the
# refit-then-FAILED rule (R/learner_ra.R) and R-b's nested rate gates (R/gates.R). Code checks on SYNTHETIC
# fixtures only: no production population is built here (the production root is exercised on the guard alone).

for (f in c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R", "dgp.R", "census.R",
            "variance.R", "controls.R", "learner_rd.R", "gates.R", "grading.R", "run_mode.R", "learner_ra.R"))
  source(test_path("..", "..", "R", f), local = TRUE)
CO <- dgp_coefficients(1.25)
ROOT <- "SYNTHETIC-ROUND3-TEST"

# ---- run modes ---------------------------------------------------------------------------------------
test_that("parse_ids reads ranges and lists; rejects junk", {
  expect_identical(parse_ids("1:3"), 1:3)
  expect_identical(parse_ids("1:400,100001:100012"), c(1:400, 100001:100012))
  expect_identical(parse_ids("5,1,3"), c(1L, 3L, 5L))
  expect_error(parse_ids("3:2"), "empty range"); expect_error(parse_ids("1:0"), "cannot read")
  expect_error(parse_ids("a"), "cannot read")
  expect_error(parse_ids("1,1"), "duplicated")
})

fake_repo <- function(spec_text, sha_ok = TRUE) {   # a frozen-design record over one frozen input
  d <- tempfile("repo-"); dir.create(file.path(d, "config", "simulation"), recursive = TRUE)
  writeLines(spec_text, file.path(d, "frozen-input.txt"))
  h <- digest::digest(file = file.path(d, "frozen-input.txt"), algo = "sha256")
  jsonlite::write_json(list(specification_version = "v9.8", specification_sha256 = strrep("a", 64),
                            frozen_files = list("frozen-input.txt" = if (sha_ok) h else strrep("0", 64))),
                       file.path(d, SIM_FROZEN_DESIGN), auto_unbox = TRUE)
  d
}
args_fun <- function(...) { a <- list(...); function(name, default) if (!is.null(a[[name]])) a[[name]] else default }

test_that("run_mode_setup: dev defaults; dry/prod need the frozen spec, explicit populations and disjoint reps", {
  ok <- fake_repo("spec v9.8"); bad <- fake_repo("spec v9.8", sha_ok = FALSE)
  d <- run_mode_setup(args_fun(), ok, B_default = 2L, R_default = 3L)
  expect_identical(d$mode, "dev"); expect_identical(d$pops, 1:2); expect_identical(d$reps, 1:3); expect_true(startsWith(d$root, "SYNTHETIC"))
  expect_error(run_mode_setup(args_fun("--mode" = "dry", "--pops" = "1:10"), bad), "frozen spec")
  expect_error(run_mode_setup(args_fun("--mode" = "dry"), ok), "explicit --pops")
  y <- run_mode_setup(args_fun("--mode" = "dry", "--pops" = "1:10", "--R" = "2"), ok, master_id = "FS")
  expect_identical(y$root, SIM_PROD_ROOT); expect_identical(y$reps, DRY_ID_BASE + 1:2)
  expect_error(run_mode_setup(args_fun("--mode" = "dry", "--pops" = "1:10", "--reps" = "1:2"), ok), "disjoint rep-ID range")
  expect_error(run_mode_setup(args_fun("--mode" = "prod", "--pops" = "1:10", "--reps" = "1000001:1000002"), ok), "rep IDs <=")
  expect_error(run_mode_setup(args_fun("--mode" = "prod", "--pops" = "1:500"), ok, master_id = "FS"), "production range")
  expect_error(run_mode_setup(args_fun("--mode" = "prod", "--pops" = "1:161"), ok, master_id = "FPC"), "production range")
  # census-only cells: reserved dry-run populations, no rep axis
  z <- run_mode_setup(args_fun("--mode" = "dry", "--pops" = "1000001:1000005"), ok, census_only = TRUE, master_id = "FPC")
  expect_identical(z$pops, 1000000L + 1:5); expect_length(z$reps, 0L)
  expect_error(run_mode_setup(args_fun("--mode" = "dry", "--pops" = "1:5"), ok, census_only = TRUE), "reserved population IDs")
  # array blocks
  t2 <- run_mode_setup(args_fun("--mode" = "prod", "--pops" = "1:25", "--task" = "3", "--block" = "10"), ok, master_id = "MS")
  expect_identical(t2$pops, 21:25)
  expect_error(run_mode_setup(args_fun("--mode" = "prod", "--pops" = "1:25", "--task" = "4", "--block" = "10"), ok), "beyond")
})

# ---- grading primitives ------------------------------------------------------------------------------
CELL <- list(regime = "many", H = 25L, m_h = 2L, J = 20L, m = 50L, f = "1/10", e2c = "none")
test_that("boot_index: deterministic, keyed, in range", {
  a <- boot_index(ROOT, CELL, "C1-meanZ", "pass-1", 7L, n_boot = 50L)
  expect_identical(a, boot_index(ROOT, CELL, "C1-meanZ", "pass-1", 7L, n_boot = 50L))
  expect_false(identical(a, boot_index(ROOT, CELL, "C1-meanZ", "pass-2", 7L, n_boot = 50L)))
  expect_true(all(a >= 1L & a <= 7L)); expect_identical(dim(a), c(50L, 7L))
})

test_that("boot_joint resamples whole population blocks with all their draws, populations equally weighted (R6)", {
  pop <- rep(1:4, times = c(5, 5, 4, 5))                    # population 3 lost a draw (a FAILED arm)
  x <- seq_along(pop) / 10
  idx <- boot_index(ROOT, CELL, "test-joint", "p", 4L, n_boot = 200L)
  seen <- list()
  stat <- function(xx, w) { seen[[length(seen) + 1L]] <<- list(n = nrow(xx), w = w); sum(w * xx[, 1]) }
  b <- boot_joint(x, pop, stat, idx)
  w0 <- pop_equal_weights(pop)
  expect_equal(as.vector(tapply(w0, pop, sum)), rep(0.25, 4))
  expect_equal(b$estimate, sum(w0 * x))
  sizes <- c(5, 5, 4, 5)
  for (k in 1:5) expect_identical(seen[[k + 1L]]$n, as.integer(sum(sizes[idx[k, ]])))   # whole blocks, never single draws
  expect_length(b$boot, 200L)
})

test_that("decide() is §6.1's operator; boot_census bootstraps populations only", {
  expect_identical(decide(c(0.01, 0.05), 0, 0.10), "PASS")
  expect_identical(decide(c(0.11, 0.2), 0, 0.10), "RED-FLAG")
  expect_identical(decide(c(0.05, 0.12), 0, 0.10), "UNGRADABLE (MC precision)")
  expect_identical(decide(c(-0.03, -0.021), lo = -0.02), "RED-FLAG")      # one-sided
  v <- c(0.02, 0.04, 0.06, 0.08, 0.05)
  b <- boot_census(v, function(z) unname(stats::quantile(z, 0.9, type = 7)), boot_index(ROOT, CELL, "bridge", "p", 5L, 300L))
  expect_equal(b$estimate, unname(stats::quantile(v, 0.9, type = 7)))
  expect_true(all(b$boot <= max(v) & b$boot >= min(v)))
  g <- grade(b, 0, 0.10); expect_identical(g$disposition, "PASS")
})

# ---- C-AE-3: the frozen record -----------------------------------------------------------------------
test_that("the frozen kappa_AE record round-trips, verifies its provenance and is never overwritten (R23)", {
  f <- file.path(tempfile("frozen-"), "kappa_ae.json")
  rec <- list(rule = AE_RULE, root = SIM_PROD_ROOT, master_id = "FS", population_ids = as.integer(AE_PANEL_POPS),
              cell = AE_PANEL_CELL, target = 0.30, kappa = 0.1337, T2 = 3.87e-3, stat = 0.3, buildable = TRUE,
              monotone = TRUE, kappa_max = 1.18, reason = "", clip = list(inside = TRUE, min_g = 0.1, max_g = 0.9, min_Q = 0.02, max_Q = 0.98),
              sigma_wor = seq(0.010, 0.012, length.out = 160), p0_primary_hash = "abc", spec_sha256 = "s1")
  ae_frozen_write(f, rec)
  got <- ae_frozen_load(f, spec_sha = "s1", grid_hash = "abc")
  expect_identical(got$kappa, 0.1337); expect_identical(got$sigma_wor, rec$sigma_wor)
  expect_error(ae_frozen_load(f, spec_sha = "s2"), "spec hash")
  expect_error(ae_frozen_load(f, grid_hash = "zzz"), "P0 grid hash")
  expect_error(ae_frozen_write(f, rec), "derived once")
  txt <- readLines(f); writeLines(sub("0.13370000000000001", "0.135", txt, fixed = TRUE), f)
  expect_error(ae_frozen_load(f), "content hash")
  expect_error(ae_frozen_load(file.path(tempdir(), "nope.json")), "derives it once")
})

# ---- production-population gate statistics, on a SYNTHETIC population --------------------------------
POP <- build_population(6, 30, CO, "MS", 1, ROOT, keep_identity = TRUE)
test_that("informativeness: the stratum-demeaned pooled correlation", {
  idx <- rep(seq_len(nrow(POP$psu)), POP$psu$M)
  pb <- as.vector(rowsum(POP$units$phi0, idx)) / POP$psu$M
  x <- POP$psu$M - ave(POP$psu$M, POP$psu$h); y <- pb - ave(pb, POP$psu$h)
  expect_equal(census_informativeness(POP), sum(x * y) / sqrt(sum(x^2) * sum(y^2)), tolerance = 1e-14)
})
test_that("DEFF_aug: m_h cancels with equal allocation; per-class values; clip check", {
  d2 <- census_deff_aug(POP, m_h = 2L); d4 <- census_deff_aug(POP, m_h = 4L)
  expect_equal(unname(d2), unname(d4), tolerance = 1e-10)
  expect_named(d2, c("pooled", "K1", "K2")); expect_true(all(is.finite(d2) & d2 > 0))
  cl <- census_clip(POP$units); expect_identical(cl[["pass"]], 1)
})
test_that("store hashes stream exactly rng_store_bytes()'s bytes", {
  tag <- "SYNTHETIC-store-test"
  hs <- population_store_hashes(POP, tag)
  idx <- rep(seq_len(nrow(POP$psu)), POP$psu$M)
  u <- POP$units
  df <- data.frame(h = as.integer(u$h), k = as.integer(u$k), i = as.integer(u$i), M = POP$psu$M[idx], KB = as.integer(POP$psu$KB[idx]),
                   Cpsu = POP$psu$Cpsu[idx], UA = POP$psu$UA[idx], UY = POP$psu$UY[idx], C1 = u$C1, C2 = u$C2, C3 = u$C3, C4 = u$C4,
                   VA = u$VA, VY = u$VY)
  expect_identical(hs[["identity"]], rng_store_hash(df, POP_IDENTITY_COLS, paste0(tag, ";tier=identity")))
  expect_identical(hs[["identity"]], population_store_hashes(build_population(6, 30, CO, "MS", 1, ROOT, keep_identity = TRUE), tag)[["identity"]])
  expect_error(population_store_hashes(build_population(2, 2, CO, "MS", 1, ROOT), tag), "keep_identity")
  # identity tier is coefficient-free (contract): a different a_cpsu leaves it unchanged, the derived tier moves
  h0 <- population_store_hashes(build_population(6, 30, dgp_coefficients(0), "MS", 1, ROOT, keep_identity = TRUE), tag)
  expect_identical(h0[["identity"]], hs[["identity"]]); expect_false(identical(h0[["derived"]], hs[["derived"]]))
})

# ---- refit once with Nelder-Mead, then ra_fit_failed ------------------------------------------------
test_that("R4: a hard fitting error is refitted once and then raises ra_fit_failed; a clean fit records no refit", {
  d <- data.frame(Y = rep(0:1, 20), A = rep(0:1, each = 20), x = seq(-1, 1, length.out = 40), psu = factor(rep(1:8, 5)))
  e <- tryCatch(.ra_glmer(Y ~ A + not_a_column + (1 | psu), d, 1L, "record", "Q[test]"), error = function(e) e)
  expect_s3_class(e, "ra_fit_failed")
  expect_match(conditionMessage(e), "FAILED after one Nelder-Mead refit")
  tr <- ra_synthetic_sample(CO, 12L, seed = 3L)
  fit <- fit_ra(tr, a_cpsu_term = TRUE)
  expect_true(all(!fit$diagnostics$refit)); expect_true(all(fit$diagnostics$optimizer == "bobyqa,bobyqa"))
})

# ---- R-b gates and G6 helpers -----------------------------------------------------------------
test_that("rb_rate_gates: the one-sided upper bound and the two nested verdicts", {
  m <- rep(c(50, 100, 200, 400, 800), each = 10); pop <- rep(1:10, 5)
  r <- exp(-0.5 * log(m) + rep(seq(-0.02, 0.02, length.out = 10), 5) + 0.001 * sin(seq_along(m)))
  g <- rb_rate_gates(r, m, pop)
  expect_true(g$pass_L2); expect_true(g$pass_quarter); expect_equal(g$slope, -0.5, tolerance = 0.01)
  g2 <- rb_rate_gates(exp(-0.1 * log(m) + 0.001 * sin(seq_along(m))), m, pop)
  expect_true(g2$pass_L2); expect_false(g2$pass_quarter)
  b <- g6_boundary(c(-0.29, 0.99, 0.30), c(-0.10, 0.85, 0), c(0.10, 1.15, 0.10), c(0.035, 0.01, 0.2))
  expect_identical(b$graded, c(TRUE, FALSE, FALSE)); expect_identical(b$inside, c(FALSE, TRUE, FALSE))
  expect_equal(g6_mcse_joint(c(1, 1, 1, 1), 1, 5), sqrt(1 / 20))
})
