# Unit tests of the grading pass:
# known-result fixtures for every statistic and the decision operator's three outcomes, the population-block bootstrap's invariants
# (blocks intact, populations equally weighted, keyed and reproducible), and the handling of incomplete,
# superseded, FAILED and wrongly-provenanced output on a SYNTHETIC production tree. No production output, no
# population, no sample: everything here is constructed numbers.

for (f in c("rng_streams.R", "run_mode.R", "learner_rd.R", "estimators.R", "tmle.R", "learner_rb.R", "controls.R",
            "grading.R", "grade_plan.R", "grade_io.R", "grade_run.R"))
  source(test_path("..", "..", "R", f), local = TRUE)
ROOT <- "SYNTHETIC-GRADING-TEST"
CELL <- list(regime = "many", H = 25L, m_h = 2L, J = 20L, m = 50L, f = "1/10", e2c = "none")
KEY <- function(s) list(root = ROOT, cell = CELL, statistic = s, pass = "test-pass")
NB <- 300L
# A deterministic standard-normal sample spread over populations (no RNG): exact normal quantiles, permuted.
zfix <- function(n) qnorm(ppoints(n))[order((seq_len(n) * 7919L) %% n)]

# ---- frequency-weighted statistics -------------------------------------------------------------------------
test_that("frequency statistics are the plain statistics at f = 1 and the expanded-sample statistics otherwise", {
  x <- c(3.1, -0.2, 1.5, 0.7, 2.2, -1.4, 0.05); f <- c(2, 1, 1, 3, 1, 2, 4); xe <- rep(x, f)
  for (p in c(0.1, 0.5, 0.9, 0.975)) {
    expect_equal(f_q7(x, rep(1, 7), p), unname(quantile(x, p, type = 7)))
    expect_equal(f_q7(x, f, p), unname(quantile(xe, p, type = 7)))
  }
  expect_equal(f_sd(x, rep(1, 7)), sd(x)); expect_equal(f_sd(x, f), sd(xe))
  expect_equal(f_mean(x, f), mean(xe)); expect_equal(f_ks(x, f), f_ks(xe, rep(1, length(xe))))
  ks_direct <- function(z) { z <- sort(z); n <- length(z); max(pmax(abs(seq_len(n) / n - pnorm(z)), abs((seq_len(n) - 1) / n - pnorm(z)))) }
  expect_equal(f_ks(x, rep(1, 7)), ks_direct(x))
  # Every population gets equal total frequency, whatever its number of retained draws.
  pop <- rep(c(11, 12, 13, 14), c(5, 5, 3, 1)); fr <- pop_freq(pop)
  expect_equal(as.vector(tapply(fr, pop, sum)), rep(15, 4))
  expect_true(all(pop_freq(rep(1:4, each = 5)) == 1))
})

# ---- the keyed bootstrap -------------------------------------------------------------------------------------
test_that("boot_stream reproduces boot_index exactly, in any chunking; the key changes the resamples", {
  idx <- boot_index(ROOT, CELL, "s1", "p1", 7L, n_boot = 40L)
  for (ch in c(1L, 3L, 40L)) {
    got <- matrix(NA_integer_, 40, 7); k <- 0L
    boot_stream(ROOT, CELL, "s1", "p1", 7L, function(i) { k <<- k + 1L; got[k, ] <<- i; 0 }, n_boot = 40L, chunk = ch)
    expect_identical(got, idx)
  }
  a <- boot_stream(ROOT, CELL, "s1", "p1", 9L, function(i) sum(i * seq_along(i)), n_boot = 25L)
  expect_identical(a, boot_stream(ROOT, CELL, "s1", "p1", 9L, function(i) sum(i * seq_along(i)), n_boot = 25L))
  expect_false(identical(a, boot_stream(ROOT, CELL, "s1", "p2", 9L, function(i) sum(i * seq_along(i)), n_boot = 25L)))
  expect_false(identical(a, boot_stream(ROOT, CELL, "s2", "p1", 9L, function(i) sum(i * seq_along(i)), n_boot = 25L)))
})

test_that("population-block bootstrap: blocks intact, never resampled within, populations equally weighted", {
  rb <- c(5L, 5L, 5L, 2L, 5L, 1L); pop <- rep(1:6, rb)
  d <- data.frame(pop_id = pop, v = as.numeric(pop), row = seq_along(pop))
  b <- boot_joint_stat(d, d$pop_id, function(d, f) f_mean(d$v, f), KEY("blk"), n_boot = 200L)
  expect_equal(b$estimate, mean(1:6))                       # equal weights, not row weights (which give 3.30)
  want <- numeric(200); k <- 0L
  boot_stream(ROOT, CELL, "blk", "test-pass", 6L, function(i) { k <<- k + 1L; want[k] <<- mean(i); 0 }, n_boot = 200L)
  expect_equal(b$boot, want)                               # each resample = mean of the selected population ids
  # every selected block brings ALL of its own rows, unchanged, as many times as it is selected
  intact <- function(d, f) {
    ok <- vapply(split(d$row, d$pop_id), function(r) { p <- d$pop_id[match(r[1], d$row)]; own <- which(pop == p)
      length(r) %% length(own) == 0 && all(table(r) == length(r) / length(own)) && setequal(r, own) }, TRUE)
    as.numeric(all(ok))
  }
  expect_true(all(boot_joint_stat(d, d$pop_id, intact, KEY("intact"), n_boot = 150L)$boot == 1))
})

# ---- known-result fixtures for 's three outcomes -----------------------------------------------------
joint_fixture <- function(shift = 0, scale = 1, B = 400L, R = 5L)
  data.frame(pop_id = rep(seq_len(B), each = R), Z_wor = shift + scale * zfix(B * R), R_var = 1 + 0.03 * zfix(B * R))
c1_grades <- function(d) vapply(c(mean_z = "mean_z", sd = "sd_z_minus1", ks = "ks_z"), function(s) {
  reg <- c(mean_z = "z_location", sd_z_minus1 = "z_scale", ks_z = "z_ks")[[s]]
  grade_boot(boot_joint_stat(d, d$pop_id, STAT[[s]]("Z_wor"), KEY(paste0("c1-", s)), NB), GRADE_MARGINS[[reg]])$disposition }, "")

test_that("C1 margins: a standard-normal Z PASSes, a shifted Z RED-FLAGs, a boundary shift is UNGRADABLE (MC precision)", {
  ok <- c1_grades(joint_fixture()); expect_identical(unname(ok), rep("PASS", 3)); expect_identical(compose(ok), "PASS")
  bad <- c1_grades(joint_fixture(shift = 0.5))
  expect_identical(unname(bad[c("mean_z", "ks")]), c("RED-FLAG", "RED-FLAG")); expect_identical(compose(bad), "RED-FLAG")
  edge <- c1_grades(joint_fixture(shift = 0.10))
  expect_identical(unname(edge[["mean_z"]]), "UNGRADABLE (MC precision)")
  expect_identical(compose(edge), "UNGRADABLE (MC precision)")
  wide <- c1_grades(joint_fixture(scale = 1.3)); expect_identical(unname(wide[["sd"]]), "RED-FLAG")
})

test_that("C2: fitted mean/SD/transfer known results; the transfer is the mean of per-replicate absolute differences", {
  d <- joint_fixture(); d$R_var_or <- 1
  g <- grade_boot(boot_joint_stat(d, d$pop_id, STAT$mean_minus1("R_var"), KEY("c2m"), NB), GRADE_MARGINS$rvar_fit_mean)
  expect_identical(g$disposition, "PASS"); expect_equal(g$estimate, mean(d$R_var) - 1)
  d$t <- d$R_var - d$R_var_or
  tr <- grade_boot(boot_joint_stat(d, d$pop_id, STAT$mean_abs("t"), KEY("c2t"), NB), GRADE_MARGINS$transfer)
  expect_equal(tr$estimate, mean(abs(d$t))); expect_identical(tr$disposition, "PASS")
  d$R_var <- 1.4 + 0.03 * zfix(nrow(d))
  expect_identical(grade_boot(boot_joint_stat(d, d$pop_id, STAT$mean_minus1("R_var"), KEY("c2m2"), NB), GRADE_MARGINS$rvar_fit_mean)$disposition, "RED-FLAG")
})

test_that("C3b coverage, C3a variance agreement and Q_0.9 negligibility: known results", {
  n <- 2000L; sm <- 1.2; sw <- 1
  q <- zfix(n); q <- (q - mean(q)) / sd(q)
  dc <- data.frame(psi_hat = 0.3 + sm * q, Psi_N = 0.3, sigma_m = sm, sigma_m_wor = sw)
  cv <- grade_boot(boot_draws_stat(dc, STAT$cov_det_sigma_m(0.95), KEY("c3b"), NB), GRADE_MARGINS$coverage)
  expect_equal(cv$estimate, mean(abs(sm * q) <= qnorm(0.975) * sm) - 0.95); expect_identical(cv$disposition, "PASS")
  dc$psi_hat <- 0.3 + 1.3 * sm * q                             # coverage P(|Z| <= 1.96/1.3) = 0.868 at the 0.95 level
  expect_identical(grade_boot(boot_draws_stat(dc, STAT$cov_det_sigma_m(0.95), KEY("c3b2"), NB), GRADE_MARGINS$coverage)$disposition, "RED-FLAG")
  d3 <- data.frame(psi_hat = 0.3 + sw * q, Psi_N = 0.3, sigma_m = sm, sigma_m_wor = sw)
  a <- grade_boot(boot_draws_stat(d3, STAT$c3a(), KEY("c3a"), NB), GRADE_MARGINS$c3a)
  expect_equal(a$estimate, 0, tolerance = 1e-12); expect_identical(a$disposition, "PASS")
  d3$psi_hat <- 0.3 + 1.3 * sw * q                             # (1.69 - 1) / 1.44 = 0.479
  b <- grade_boot(boot_draws_stat(d3, STAT$c3a(), KEY("c3a2"), NB), GRADE_MARGINS$c3a)
  expect_equal(b$estimate, 0.69 / 1.44, tolerance = 1e-12); expect_identical(b$disposition, "RED-FLAG")
  dq <- data.frame(pop_id = rep(1:400, each = 5), t = 0.05 * ppoints(2000)[order((1:2000 * 13L) %% 2000)])
  qa <- grade_boot(boot_joint_stat(dq, dq$pop_id, STAT$q90_abs("t"), KEY("q90"), NB), GRADE_MARGINS$negligible)
  expect_equal(qa$estimate, unname(quantile(dq$t, 0.9, type = 7))); expect_identical(qa$disposition, "PASS")
  dq$t <- dq$t + 0.25
  expect_identical(grade_boot(boot_joint_stat(dq, dq$pop_id, STAT$q90_abs("t"), KEY("q90b"), NB), GRADE_MARGINS$negligible)$disposition, "RED-FLAG")
  # census-only: the bridge's Q_0.9 over populations (type 7)
  v <- seq(0.01, 0.06, length.out = 120)
  bc <- grade_boot(boot_census_stat(v, function(z) unname(quantile(z, 0.9, type = 7)), KEY("bridge"), NB), GRADE_MARGINS$negligible)
  expect_equal(bc$estimate, unname(quantile(v, 0.9, type = 7))); expect_identical(bc$disposition, "PASS")
  # a statistic with no MC spread gives a degenerate interval and an exact decision
  const <- data.frame(pop_id = rep(1:10, each = 2), t = 0.2)
  expect_identical(grade_boot(boot_joint_stat(const, const$pop_id, STAT$q90_abs("t"), KEY("const"), 50L), GRADE_MARGINS$negligible)$disposition, "RED-FLAG")
})

test_that("composition (§6.1a) and power-control dispositions", {
  expect_identical(compose(c("PASS", "PASS")), "PASS")
  expect_identical(compose(c("PASS", "RED-FLAG", "UNGRADABLE (MC precision)")), "RED-FLAG")
  expect_identical(compose(c("PASS", "UNGRADABLE (MC precision)")), "UNGRADABLE (MC precision)")
  expect_identical(compose(c("UNGRADABLE (premise: x)", "PASS")), "UNGRADABLE (premise: x)")
  expect_identical(compose(c("NOT GRADED (incomplete output: y)", "RED-FLAG")), "NOT GRADED (incomplete output: y)")
  expect_identical(control_disposition("RED-FLAG"), "FAILS AS PREDICTED")
  expect_match(control_disposition("PASS"), "^NO POWER")
  expect_identical(control_disposition("UNGRADABLE (MC precision)"), "UNGRADABLE (MC precision)")
})

# ---- a synthetic production tree: completeness, superseded reruns, FAILED arms, provenance -------------------
FROZEN <- strrep("a", 64); COMMIT <- "c0ffee0"
# Required-check rows a production run writes for its results: U1 per replicate (arm "all") and U8 per replicate
# and arm marked OK; census: the FPC algebra per population and rung; E4: U1, U8 and finite outputs per rung.
slice_checks <- function(res) {
  reps <- unique(res[, c("pop_id", "rep_id")])
  okr <- res[res$arm_status %in% "OK", c("pop_id", "rep_id", "arm")]
  rbind(data.frame(pop_id = reps$pop_id, rep_id = reps$rep_id, arm = "all", check = "U1 weights bit-exact", value = 1, reference = 1,
                   abs_diff = NA, tolerance = NA, pass = TRUE, required = TRUE, stringsAsFactors = FALSE),
        if (nrow(okr)) data.frame(pop_id = okr$pop_id, rep_id = okr$rep_id, arm = okr$arm, check = "U8 se2 = survey::svymean", value = 1, reference = 1,
                                  abs_diff = 0, tolerance = 1e-12, pass = TRUE, required = TRUE, stringsAsFactors = FALSE))
}
census_checks <- function(tr) data.frame(pop_id = tr$pop_id, H = tr$H, alloc = NA, check = "census FPC algebra: sigma2_m - sigma2_m_wor = Delta_m",
                                         value = 0, reference = 0, pass = TRUE, required = TRUE, stringsAsFactors = FALSE)
write_run <- function(dir, results = NULL, traj = NULL, cell = NULL, arms = NULL, spec = FROZEN, census = FALSE, checks = NULL, declared = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(checks)) checks <- if (!is.null(traj)) census_checks(traj) else if (!is.null(results)) slice_checks(results) else data.frame()
  n_req <- if (nrow(checks)) sum(checks$required) else 0L
  if (!is.null(declared)) n_req <- declared
  man <- if (census) list(mode = "prod", spec_freeze = list(actual = spec, ok = TRUE), package_commit = COMMIT, config = list(root = ROOT),
                          checks = list(n_required = n_req))
         else list(mode = "prod", spec_freeze = list(sha256 = spec), git = list(package_commit = COMMIT), config = list(root = ROOT, cell = cell, arms = arms),
                   checks = list(n_required = n_req))
  jsonlite::write_json(man, file.path(dir, "manifest.json"), auto_unbox = TRUE, digits = I(17))
  if (!is.null(results)) utils::write.csv(results, file.path(dir, "results.csv"), row.names = FALSE)
  if (!is.null(traj)) utils::write.csv(traj, file.path(dir, "census_trajectory.csv"), row.names = FALSE)
  if (nrow(checks)) utils::write.csv(checks, file.path(dir, "checks.csv"), row.names = FALSE)
}
fitted_rows <- function(pops, R, arms, shift = 0, reps = seq_len(R)) {
  g <- expand.grid(rep_id = reps, pop_id = pops, arm = arms, stringsAsFactors = FALSE)
  n <- nrow(g); z <- zfix(n)
  data.frame(pop_id = g$pop_id, rep_id = g$rep_id, arm = g$arm, psi_hat = 0.3 + 0.01 * (z + shift), se_hat = 0.01, Psi_N = 0.3,
             sigma_m = 0.011, sigma_m_wor = 0.01, Z = z + shift, cover95 = abs(z) < 1.96, Z_wor = z + shift,
             R_var = 1 + 0.03 * zfix(n), arm_status = "OK", paired_tmle_minus_aipw = 1e-5 * z, se2_diff_over_sigma2_m = 1e-3 * abs(z),
             gram_premise_met = TRUE, lambda_min = 0.2, T1_over_sigma_wor = 0.01 * z, T2_over_sigma_wor = 0.002 * z, T1 = 1e-4 * z, T2 = 2e-5 * z,
             stringsAsFactors = FALSE)
}
fail_rows <- function(res, arm, pop, rep) {
  i <- which(res$arm == arm & res$pop_id == pop & res$rep_id == rep)
  res$arm_status[i] <- "FAILED"; res$psi_hat[i] <- NA; res$Z_wor[i] <- NA; res$R_var[i] <- NA; res
}
CK <- function(H, m) list(regime = "many", H = H, m_h = 2L, J = 20L, m = m, f = "1/10", e2c = "none")
MS5 <- c("oracle", "Ra-AIPW", "Ra-TMLE", "Rb-AIPW", "Rb-TMLE")
# The test-sized frozen design for the tree (the production table is frozen_cells(); tested against cells.tsv below).
TEST_FROZEN <- list(
  "P-ms-H800" = list(pops = 1:400, R = 5L, arms = c("oracle", "Ra-AIPW"), rungs = NULL),
  "P-ms-H400" = list(pops = 1:20, R = 5L, arms = MS5, rungs = NULL),
  "P-fs-m1600" = list(pops = 1:4, R = 5L, arms = c("oracle", "Ra-AIPW"), rungs = NULL),
  "P-fpc-or-J4" = list(pops = 1:2, R = 40L, arms = "oracle", rungs = NULL),
  "P-vanish-census" = list(pops = 1:20, R = 0L, arms = character(0), rungs = data.frame(H = c(25L, 200L), J = c(10L, 800L))))
TEST_CELLS <- data.frame(id = c("P-ms-H800", "P-ms-H400", "P-fs-m1600", "P-fpc-or-J4", "P-vanish-census"),
                         minutes = "30", tests = "0", script = c("slice", "slice", "slice", "slice", "census"),
                         args = c("--mode prod --regime many --H 800 --pops 1:400 --R 5", "--mode prod --regime many --H 400 --pops 1:20 --R 5",
                                  "--mode prod --regime fixed --m 1600 --pops 1:4 --R 5", "--mode prod --regime many --e2c fpc --J 4 --oracle-only 1 --pops 1:2 --R 40",
                                  "--traj vanish --mode prod --pops 1:20"),
                         mem = "8G", tasks = c("2", "2", "1", "2", "2"), block = c("200", "10", "1", "1", "10"), after = "-", dry_ref = "-", reps = "5",
                         stringsAsFactors = FALSE)
build_r3 <- function(r3) {
  ad <- file.path(r3, "assemble", "assemble-20261012T120000Z"); dir.create(ad, recursive = TRUE)
  jsonlite::write_json(list(verdict = "PASS", package_commit = COMMIT), file.path(ad, "execution_gates.json"), auto_unbox = TRUE)
  utils::write.csv(data.frame(alpha = RD_ALPHAS, pass = TRUE), file.path(ad, "g5.csv"), row.names = FALSE)
  utils::write.csv(data.frame(regime = "many", rungs = 5, populations = 10, slope = -0.5, upper_95 = -0.40, L2 = "x", quarter = "y"), file.path(ad, "rb_gates.csv"), row.names = FALSE)
  pd <- file.path(r3, "prep-panel", "prep-panel-20261012T100000Z"); dir.create(pd, recursive = TRUE)
  utils::write.csv(data.frame(control = "C1 shift: mean Z", predicted = -0.27, lo = -0.1, hi = 0.1, inside = FALSE, graded = TRUE, status = "graded"), file.path(pd, "g6.csv"), row.names = FALSE)
  jsonlite::write_json(list(content_sha256 = "k", buildable = TRUE), file.path(pd, "kappa_ae.json"), auto_unbox = TRUE)
}
# H800: complete (400 x 5, two tasks); ONE FAILED Ra-AIPW replicate; an older superseded folder per task.
write_h800 <- function(prod, fail = TRUE, fail_in_superseded = FALSE, mutate = NULL) {
  for (t in 1:2) {
    pops <- (t - 1L) * 200L + 1:200
    res <- fitted_rows(pops, 5L, c("oracle", "Ra-AIPW"))
    if (fail && t == 1L) res <- fail_rows(res, "Ra-AIPW", 3, 2)
    if (!is.null(mutate) && t == 1L) res <- mutate(res)
    write_run(file.path(prod, "P-ms-H800", sprintf("prod-e2e-x-t%03d-20261013T000000Z", t)), res, cell = CK(800L, 1600L), arms = c("oracle", "Ra-AIPW"))
    old <- res; old$Z_wor <- 99
    if (fail_in_superseded && t == 2L) old <- fail_rows(old, "Ra-AIPW", 250, 4)
    write_run(file.path(prod, "P-ms-H800", sprintf("prod-e2e-x-t%03d-20261012T000000Z", t)), old, cell = CK(800L, 1600L), arms = c("oracle", "Ra-AIPW"))
  }
}
build_tree <- function() {
  top <- tempfile("grade-tree-"); prod <- file.path(top, "production"); r3 <- file.path(top, "round3")
  write_h800(prod)
  # P-ms-H400: task 2 missing -> INCOMPLETE, NOT GRADED; C6 therefore NOT GRADED
  write_run(file.path(prod, "P-ms-H400", "prod-e2e-x-t001-20261013T000000Z"), fitted_rows(1:10, 5L, MS5), cell = CK(400L, 800L), arms = MS5)
  # P-fs-m1600: wrong specification identifier -> REFUSED
  write_run(file.path(prod, "P-fs-m1600", "prod-e2e-x-t001-20261013T000000Z"), fitted_rows(1:4, 5L, c("oracle", "Ra-AIPW")),
            cell = list(regime = "fixed", H = 5L, m_h = 320L, J = 3200L, m = 1600L, f = "1/10", e2c = "none"), arms = c("oracle", "Ra-AIPW"), spec = strrep("b", 64))
  # P-fpc-or-J4: 2 populations x 40 draws, one task each
  for (p in 1:2) {
    res <- fitted_rows(p, 40L, "oracle"); res$psi_hat <- 0.3 + 0.01 * res$Z_wor
    write_run(file.path(prod, "P-fpc-or-J4", sprintf("prod-e2e-x-t%03d-20261013T000000Z", p)), res,
              cell = list(regime = "many", H = 50L, m_h = 2L, J = 4L, m = 100L, f = "1/2", e2c = "fpc"), arms = "oracle")
  }
  # P-vanish-census: 20 populations x 2 frozen rungs, 2 tasks
  for (t in 1:2) {
    pops <- (t - 1L) * 10L + 1:10
    tr <- rbind(data.frame(pop_id = pops, H = 25L, J = 10L, bridge = 0.2, terminal = FALSE),
                data.frame(pop_id = pops, H = 200L, J = 800L, bridge = seq(0.01, 0.05, length.out = 10), terminal = TRUE))
    tr$Delta_m <- 1e-5; tr$Psi_N_minus_Psi <- 1e-3 * zfix(nrow(tr)); tr$Delta_over_sigma2_m <- 0.1; tr$sigma_ratio <- 1.05
    write_run(file.path(prod, "P-vanish-census", sprintf("prod-e2e-vanish-t%03d-20261013T000000Z", t)), traj = tr, census = TRUE)
  }
  build_r3(r3)
  list(prod = prod, r3 = r3, cells = TEST_CELLS, out = file.path(top, "out"))
}
quiet <- function(...) invisible(NULL)
run_tree <- function(tr, out = tr$out, cells = tr$cells) grade_run(tr$prod, tr$r3, cells, out, FROZEN, n_boot = 200L, root = ROOT, frozen = TEST_FROZEN, say = quiet)
rd_cell <- function(tr, i, frozen = TEST_FROZEN) read_prod_cell(tr$prod, tr$cells[i, ], FROZEN, ROOT, frozen)

test_that("the frozen design table is the production cells.tsv, cell for cell (populations and R)", {
  pc <- utils::read.delim(test_path("..", "..", "..", "..", "config", "simulation", "production-cells.tsv"), stringsAsFactors = FALSE, colClasses = "character")
  fz <- frozen_cells(ae_buildable = TRUE)
  expect_setequal(names(fz), pc$id)
  for (i in seq_len(nrow(pc))) expect_identical(cell_expectation(pc[i, ], fz)$mismatch, "", info = pc$id[i])
  expect_true(ARM_AE_AIPW %in% fz[["P-e3-ctrl-m400"]]$arms)
  expect_false(ARM_AE_AIPW %in% frozen_cells(ae_buildable = FALSE)[["P-e3-ctrl-m400"]]$arms)
  # The R-d rungs alpha = 0.20 and 0.35 are STOPPED after G5 and not run in production
  expect_identical(fz[["P-fs-m3200"]]$arms, c(ARM_ORACLE, vapply(c(0.25, 0.30), rd_arm, "")))
  expect_identical(fz[["P-fs-m200"]]$arms, c(ARM_ORACLE, ARM_RA_AIPW, ARM_RA_TMLE, vapply(c(0.25, 0.30), rd_arm, "")))
  expect_identical(grade_rd_run(), c(0.25, 0.30))
  expect_identical(fz[["P-ms-H800"]]$arms, c(ARM_ORACLE, ARM_RA_AIPW))
  bad <- pc[pc$id == "P-ms-H25", ]; bad$args <- sub("--pops 1:400", "--pops 1:399", bad$args)
  expect_match(cell_expectation(bad, fz)$mismatch, "populations")
})

test_that("completeness: superseded reruns ignored, a missing task found, a wrong spec refused, FAILED counted", {
  tr <- build_tree()
  a <- rd_cell(tr, 1); expect_identical(cell_completeness(a)$status, "COMPLETE")
  expect_length(a$folders$used, 2L); expect_length(a$folders$superseded, 2L)
  expect_false(any(a$results$Z_wor == 99, na.rm = TRUE))            # the newest folder of each task wins
  ar <- arm_rows(a$results, "Ra-AIPW", "Z_wor"); expect_identical(ar$failed, 1L); expect_identical(nrow(ar$rows), 1999L)
  expect_identical(nrow(a$failed), 1L)
  b <- cell_completeness(rd_cell(tr, 2)); expect_identical(b$status, "INCOMPLETE"); expect_identical(b$missing_tasks, 2L)
  expect_identical(b$missing_pops, 11:20)
  expect_identical(cell_completeness(rd_cell(tr, 3))$status, "REFUSED")
  expect_identical(cell_completeness(rd_cell(tr, 5))$status, "COMPLETE")
  d <- rd_cell(tr, 1); d$results <- rbind(d$results, d$results[1, ]); expect_identical(cell_completeness(d)$status, "INCOMPLETE")
})

test_that("completeness is judged against the FROZEN arms, rungs, populations and replicate IDs", {
  tr <- build_tree()
  base <- rd_cell(tr, 1)
  x <- base; x$results <- x$results[x$results$arm != "Ra-AIPW", ]               # an expected arm absent
  cx <- cell_completeness(x); expect_identical(cx$status, "INCOMPLETE"); expect_match(cx$detail, "missing arms \\[Ra-AIPW\\]")
  x <- base; extra <- x$results[x$results$arm == "oracle", ]; extra$arm <- "Rb-AIPW"; x$results <- rbind(x$results, extra)
  expect_match(cell_completeness(x)$detail, "unexpected arms \\[Rb-AIPW\\]")
  x <- base; x$results$rep_id <- x$results$rep_id + 1L                          # reps 2..6, not the frozen 1..5
  cx <- cell_completeness(x); expect_identical(cx$status, "INCOMPLETE"); expect_match(cx$detail, "without exactly reps 1..5")
  x <- base; x$results <- x$results[!(x$results$pop_id == 7 & x$results$rep_id == 3 & x$results$arm == "oracle"), ]
  expect_identical(cell_completeness(x)$missing_pops, 7L)
  x <- base; x$results <- x$results[x$results$pop_id != 400, ]
  cx <- cell_completeness(x); expect_identical(cx$missing_pops, 400L); expect_identical(cx$missing_tasks, 2L)
  v <- rd_cell(tr, 5); v$trajectory <- v$trajectory[!(v$trajectory$pop_id == 4 & v$trajectory$H == 200L), ]   # a census rung missing
  cv <- cell_completeness(v); expect_identical(cv$status, "INCOMPLETE"); expect_identical(cv$missing_pops, 4L)
  e <- list(expect = list(script = "e4", rungs = data.frame(H = c(25L, 50L)), R = 3L, mismatch = "", pops = 1L, tasks = 1L, block = 1L),
            folders = list(used = "x", superseded = character(0)), provenance = data.frame(),
            results = data.frame(H = c(25L, 25L, 25L), rep_id = 1:3))
  ce <- cell_completeness(e); expect_identical(ce$status, "INCOMPLETE"); expect_match(ce$detail, "missing rungs 50")
  e$results <- data.frame(H = rep(c(25L, 50L), each = 3), rep_id = c(1:3, 2:4)); expect_match(cell_completeness(e)$detail, "reps 1..3: 50")
  e$results <- data.frame(H = rep(c(25L, 50L), each = 3), rep_id = rep(1:3, 2)); expect_identical(cell_completeness(e)$status, "COMPLETE")
  # the production cells.tsv must agree with the frozen design, or the cell is refused
  wrong <- tr$cells; wrong$args[2] <- sub("--pops 1:20", "--pops 1:30", wrong$args[2])
  m <- read_prod_cell(tr$prod, wrong[2, ], FROZEN, ROOT, TEST_FROZEN); cm <- cell_completeness(m)
  expect_identical(cm$status, "REFUSED"); expect_match(cm$detail, "design mismatch")
})

test_that("required-check audit: failed, missing and truncated required checks block the grades they affect", {
  tr <- build_tree(); base <- rd_cell(tr, 1)
  expect_length(required_check_audit(base)$all, 0L); expect_length(required_check_audit(base)$arms, 0L)
  x <- base; i <- which(x$checks$arm == "Ra-AIPW")[5]; x$checks$pass[i] <- FALSE
  au <- required_check_audit(x); expect_length(au$all, 0L); expect_match(au$arms[["Ra-AIPW"]], "failed required check: U8")
  x <- base; i <- which(x$checks$arm == "all")[1]; x$checks$pass[i] <- NA          # NA is not a pass
  expect_match(required_check_audit(x)$all, "failed required check: U1")
  x <- base; x$checks <- x$checks[!(x$checks$arm == "Ra-AIPW" & x$checks$pop_id == 9 & x$checks$rep_id == 1), ]
  au <- required_check_audit(x); expect_match(au$arms[["Ra-AIPW"]], "no U8 rows for 1 replicate"); expect_match(au$all, "holds")
  x <- base; x$checks <- x$checks[!(x$checks$arm == "all" & x$checks$pop_id == 9 & x$checks$rep_id == 1), ]
  expect_true(any(grepl("no U1 rows for 1 replicate", required_check_audit(x)$all)))
  x <- base; x$required_declared[1] <- x$required_declared[1] + 5                 # checks.csv truncated against its manifest
  expect_match(required_check_audit(x)$all, "missing required checks: .* holds")
  v <- rd_cell(tr, 5); v$checks <- v$checks[!(v$checks$pop_id == 3 & v$checks$H == 25L), ]; v$required_declared <- rep(NA_real_, 2)
  expect_match(required_check_audit(v)$all, "without the census FPC algebra")
  e <- list(expect = list(script = "e4", rungs = data.frame(H = c(25L, 50L))), folders = list(used = "f"), required_declared = NA_real_,
            results = data.frame(), checks = data.frame(H = c(25L, 25L, 25L, 50L), check = c("U1 a", "U8 b", "finite outputs: c", "U1 a"),
                                                       required = TRUE, pass = TRUE, run_folder = "f"))
  ae <- required_check_audit(e); expect_length(ae$all, 0L); expect_match(ae$rungs[["50"]], "U8, finite outputs")
})

test_that("R25: FAILED fits make the arm's grades UNGRADABLE (fitting failures), survivors descriptive; failures preserved", {
  tr <- build_tree(); res <- run_tree(tr); g <- res$grades
  c1 <- g[g$claim == "C1" & g$cell == "P-ms-H800", ]
  expect_identical(nrow(c1), 3L)
  expect_true(all(startsWith(c1$disposition, "UNGRADABLE (fitting failures: Ra-AIPW 1 FAILED")))
  expect_true(all(c1$basis == "descriptive (survivors only, R25)")); expect_true(all(is.finite(c1$estimate) & is.finite(c1$ci_lo)))
  c2 <- g[g$claim == "C2" & g$cell == "P-ms-H800", ]
  expect_true(all(startsWith(c2$disposition, "UNGRADABLE (fitting failures")))
  cl <- res$claims; expect_match(cl$disposition[cl$claim == "C1" & cl$cell == "P-ms-H800"], "^UNGRADABLE \\(fitting failures")
  expect_identical(res$failed_arms$failed[res$failed_arms$cell == "P-ms-H800" & res$failed_arms$arm == "Ra-AIPW"], 1L)
  # a failure recorded only in a superseded folder is still a failure (never erased by a rerun)
  top <- tempfile("grade-sup-"); prod <- file.path(top, "production"); r3 <- file.path(top, "round3")
  write_h800(prod, fail = FALSE, fail_in_superseded = TRUE); build_r3(r3)
  t2 <- list(prod = prod, r3 = r3, cells = TEST_CELLS[1, ], out = file.path(top, "out"))
  x <- rd_cell(t2, 1); expect_identical(nrow(x$failed), 1L); expect_identical(x$failed_superseded_only, 1L)
  r2 <- run_tree(t2); expect_true(all(startsWith(r2$grades$disposition[r2$grades$claim == "C1"], "UNGRADABLE (fitting failures")))
  # with no failure anywhere, the same cell is graded
  top <- tempfile("grade-ok-"); prod <- file.path(top, "production"); r3 <- file.path(top, "round3")
  write_h800(prod, fail = FALSE); build_r3(r3)
  r3g <- run_tree(list(prod = prod, r3 = r3, cells = TEST_CELLS[1, ], out = file.path(top, "out")))
  expect_true(all(r3g$grades$disposition[r3g$grades$claim == "C1"] == "PASS")); expect_true(all(r3g$grades$basis[r3g$grades$claim == "C1"] == "graded"))
})

test_that("blocking end to end: a failed required check and a non-finite value marked OK block the affected grades", {
  top <- tempfile("grade-blk-"); prod <- file.path(top, "production"); r3 <- file.path(top, "round3")
  write_h800(prod, fail = FALSE, mutate = function(res) { i <- which(res$arm == "Ra-AIPW" & res$pop_id == 5 & res$rep_id == 1); res$Z_wor[i] <- NA; res })
  build_r3(r3)
  r <- run_tree(list(prod = prod, r3 = r3, cells = TEST_CELLS[1, ], out = file.path(top, "out")))
  c1 <- r$grades[r$grades$claim == "C1", ]
  expect_true(all(startsWith(c1$disposition, "BLOCKED (non-finite required value in 1 replicate(s) marked OK")))
  expect_true(all(startsWith(r$grades$disposition[r$grades$claim == "C2"], "BLOCKED")))
  expect_match(r$claims$disposition[r$claims$claim == "C1"], "^BLOCKED")
  # a failed required U8 row of Ra-AIPW: every Ra-AIPW grade of the cell is BLOCKED
  top <- tempfile("grade-blk2-"); prod <- file.path(top, "production"); r3 <- file.path(top, "round3")
  write_h800(prod, fail = FALSE); build_r3(r3)
  f <- list.files(file.path(prod, "P-ms-H800"), pattern = "^checks\\.csv$", recursive = TRUE, full.names = TRUE)
  f <- f[grepl("20261013T", f)][1]; ck <- utils::read.csv(f, stringsAsFactors = FALSE)
  ck$pass[which(ck$arm == "Ra-AIPW")[1]] <- FALSE; utils::write.csv(ck, f, row.names = FALSE)
  r <- run_tree(list(prod = prod, r3 = r3, cells = TEST_CELLS[1, ], out = file.path(top, "out")))
  gg <- r$grades[r$grades$cell == "P-ms-H800", ]
  expect_true(all(startsWith(gg$disposition, "BLOCKED (failed required check: U8")))
  expect_true(nrow(r$blocked) >= 1L)
})

test_that("grade_run end to end on the synthetic tree: NOT GRADED cells, C6 gating, the bridge and C3a; keyed and reproducible", {
  tr <- build_tree()
  res <- run_tree(tr)
  cmp <- setNames(res$completeness$status, res$completeness$cell)
  expect_identical(unname(cmp[c("P-ms-H800", "P-ms-H400", "P-fs-m1600", "P-fpc-or-J4", "P-vanish-census")]),
                   c("COMPLETE", "INCOMPLETE", "REFUSED", "COMPLETE", "COMPLETE"))
  g <- res$grades
  expect_true(all(startsWith(g$disposition[g$cell == "P-ms-H400" & g$claim %in% c("C1", "C2")], "NOT GRADED")))
  expect_true(all(startsWith(g$disposition[g$claim == "C6"], "NOT GRADED")))
  expect_true(all(startsWith(g$disposition[g$cell == "P-fs-m1600"], "NOT GRADED")))
  expect_identical(nrow(g[g$claim == "C3a", ]), 2L)
  expect_identical(nrow(g[g$claim == "C3b", ]), 6L)                 # 2 populations x 3 levels
  br <- g[g$claim == "sec. Psi", ]; expect_identical(br$disposition, "PASS")
  expect_equal(br$estimate, unname(quantile(rep(seq(0.01, 0.05, length.out = 10), 2), 0.9, type = 7)))
  pc <- res$provenance_checks; expect_true(pc$ok[1]); expect_false(pc$ok[2]); expect_true(pc$ok[4])
  expect_true(file.exists(file.path(tr$out, "grades.csv"))); expect_true(file.exists(file.path(tr$out, "displays", "e2_vanishing_trajectory.csv")))
  res2 <- run_tree(tr, out = file.path(tr$out, "again"))
  expect_identical(res2$grades[, c("estimate", "ci_lo", "ci_hi", "disposition")], g[, c("estimate", "ci_lo", "ci_hi", "disposition")])
})

# ---- stopped R-d rungs, C5, and the UNRESOLVED-T2 rule ---------------------
test_that("an UNRESOLVED T2 is counted apart (never non-finite-OK) and only for the grades that use T1 or T2", {
  r <- fitted_rows(1:4, 5L, "Rd-AIPW-a030"); r$T2_status <- "primary"
  i <- which(r$pop_id == 2 & r$rep_id == 3)
  r$T2_status[i] <- "unresolved"; r$T2[i] <- NA; r$T1[i] <- NA; r$T2_over_sigma_wor[i] <- NA; r$T1_over_sigma_wor[i] <- NA
  a <- arm_rows(r, "Rd-AIPW-a030", c("T1_over_sigma_wor", "T2_over_sigma_wor"))
  expect_identical(a$unresolved, 1L); expect_identical(a$nonfinite_ok, 0L); expect_identical(nrow(a$rows), 19L)
  z <- arm_rows(r, "Rd-AIPW-a030", "Z_wor")                        # an estimation grade: the replicate is used
  expect_identical(z$unresolved, 0L); expect_identical(nrow(z$rows), 20L)
  r$T2[i + 1L] <- NA; r$T2_over_sigma_wor[i + 1L] <- NA              # a missing T2 NOT marked unresolved still blocks
  expect_identical(arm_rows(r, "Rd-AIPW-a030", c("T1_over_sigma_wor", "T2_over_sigma_wor"))$nonfinite_ok, 1L)
  expect_identical(compose(c("UNGRADABLE (unresolved T2: 1 replicate(s) of x)", "UNGRADABLE (MC precision)")), "UNGRADABLE (unresolved T2: 1 replicate(s) of x)")
})

test_that("grade_run under R27: stopped rungs UNGRADABLE (premise), C5 UNGRADABLE, an unresolved T2 makes its remainder grade UNGRADABLE", {
  top <- tempfile("grade-r3d-"); prod <- file.path(top, "production"); r3 <- file.path(top, "round3")
  arms <- c("oracle", "Rd-AIPW-a025", "Rd-AIPW-a030")
  res <- fitted_rows(1:4, 5L, arms); res$T2_status <- ifelse(res$arm == "oracle", NA, "primary")
  i <- which(res$arm == "Rd-AIPW-a030" & res$pop_id == 2 & res$rep_id == 3)
  res$T2_status[i] <- "unresolved"; res$T2[i] <- NA; res$T1[i] <- NA; res$T2_over_sigma_wor[i] <- NA; res$T1_over_sigma_wor[i] <- NA
  write_run(file.path(prod, "P-fs-m3200", "prod-e2e-x-t001-20261020T000000Z"), res,
            cell = list(regime = "fixed", H = 5L, m_h = 640L, J = 6400L, m = 3200L, f = "1/10", e2c = "none"), arms = arms)
  build_r3(r3)
  fz <- list("P-fs-m3200" = list(pops = 1:4, R = 5L, arms = arms, rungs = NULL))
  cells <- data.frame(id = "P-fs-m3200", minutes = "30", tests = "0", script = "slice", args = "--mode prod --regime fixed --m 3200 --pops 1:4 --R 5 --rd-alphas 0.25,0.30",
                      mem = "8G", tasks = "1", block = "4", after = "-", dry_ref = "-", reps = "5", stringsAsFactors = FALSE)
  out <- grade_run(prod, r3, cells, file.path(top, "out"), FROZEN, n_boot = 200L, root = ROOT, frozen = fz, say = quiet)
  g <- out$grades
  for (al in c("Rd-AIPW-a020", "Rd-AIPW-a035")) {
    x <- g[g$arm == al & g$claim %in% c("C4a", "C4b"), ]
    expect_identical(nrow(x), 2L)
    expect_true(all(grepl("^UNGRADABLE [(]premise: G5 alpha 0[.](20|35) STOPPED after G5 .*ruling R27", x$disposition)))
  }
  c4a30 <- g[g$arm == "Rd-AIPW-a030" & g$claim == "C4a", ]
  expect_true(startsWith(c4a30$disposition, "UNGRADABLE (unresolved T2: 1 replicate(s) of Rd-AIPW-a030"))
  expect_match(c4a30$basis, "descriptive"); expect_true(is.finite(c4a30$estimate)); expect_identical(c4a30$rows, 19L)
  expect_match(g$disposition[g$arm == "Rd-AIPW-a025" & g$claim == "C4b"], "^BOUNDARY")
  c4a25 <- g[g$arm == "Rd-AIPW-a025" & g$claim == "C4a", ]
  expect_false(grepl("unresolved", c4a25$disposition)); expect_identical(c4a25$rows, 20L)
  c5 <- g[g$claim == "C5", ]
  expect_match(c5$disposition, "only premise-violating rung, alpha < 1/4 (spec :128), is 0.20, STOPPED after G5 (ruling R27); alpha = 0.25 is BOUNDARY", fixed = TRUE)
  expect_identical(c5$arm, "R-d (alpha < 1/4)")
  # the display counts the unresolved replicate
  d <- out$displays$e1_fixed_strata_by_arm
  expect_identical(as.integer(d$t2_unresolved[d$arm == "Rd-AIPW-a030"]), 1L)
})

test_that("ruling R28: Ra-TMLE's T1 and T2 enter no display and no grade; Ra-AIPW's stay", {
  top <- tempfile("grade-r28-"); prod <- file.path(top, "production"); r3 <- file.path(top, "round3")
  arms <- c("oracle", "Ra-AIPW", "Ra-TMLE")
  res <- fitted_rows(1:4, 5L, arms); res$T2_status <- ifelse(res$arm == "oracle", NA, "primary")
  write_run(file.path(prod, "P-fs-m200", "prod-e2e-x-t001-20261020T000000Z"), res,
            cell = list(regime = "fixed", H = 5L, m_h = 40L, J = 400L, m = 200L, f = "1/10", e2c = "none"), arms = arms)
  build_r3(r3)
  fz <- list("P-fs-m200" = list(pops = 1:4, R = 5L, arms = arms, rungs = NULL))
  cells <- data.frame(id = "P-fs-m200", minutes = "30", tests = "0", script = "slice", args = "--mode prod --regime fixed --m 200 --pops 1:4 --R 5",
                      mem = "8G", tasks = "1", block = "4", after = "-", dry_ref = "-", reps = "5", stringsAsFactors = FALSE)
  out_dir <- file.path(top, "out")
  out <- grade_run(prod, r3, cells, out_dir, FROZEN, n_boot = 200L, root = ROOT, frozen = fz, say = quiet)
  tcols <- c("q50_abs_T1_std", "q90_abs_T1_std", "q50_abs_T2_std", "q90_abs_T2_std", "q50_sqrt_m_abs_T1", "q50_sqrt_m_abs_T2")
  for (d in list(out$displays$e1_fixed_strata_by_arm,
                 utils::read.csv(file.path(out_dir, "displays", "e1_fixed_strata_by_arm.csv"), stringsAsFactors = FALSE))) {
    expect_true(all(tcols %in% names(d)))
    tm <- d[d$arm == "Ra-TMLE", , drop = FALSE]; ra <- d[d$arm == "Ra-AIPW", , drop = FALSE]
    expect_identical(nrow(tm), 1L); expect_identical(nrow(ra), 1L)
    expect_true(all(is.na(unlist(tm[, tcols]))))
    expect_true(is.na(tm$t2_unresolved))
    expect_identical(tm$remainder_note, GRADE_R28_NOTE)
    expect_true(all(is.finite(unlist(ra[, tcols]))))                     # Ra-AIPW's remainders are still displayed
    expect_true(is.na(ra$remainder_note))
    expect_true(is.finite(tm$mean_Z))                                    # Ra-TMLE's estimation summaries stay
  }
  g <- out$grades
  expect_false(any(g$arm == "Ra-TMLE" & grepl("T1|T2", g$statistic)))
  expect_true(grade_remainder_excluded("Ra-TMLE")); expect_false(grade_remainder_excluded("Ra-AIPW"))
})
