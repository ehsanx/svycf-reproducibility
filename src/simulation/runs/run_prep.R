# PREPARATION: census algebra on the PRODUCTION populations, before any fitted cell runs.
#
#   Rscript src/simulation/env/run.R src/simulation/runs/run_prep.R --job census --master MS|FS|FPC --pops a:b [--task i --block k]
#   Rscript src/simulation/env/run.R src/simulation/runs/run_prep.R --job panel
#   Rscript src/simulation/env/run.R src/simulation/runs/run_prep.R --job g1g4
#   Rscript src/simulation/env/run.R src/simulation/runs/run_prep.R --job summary --in <prep output root>
#
# --job census   every listed production population built ONCE on its whole master frame (MS 800 x 20, FS 5 x 6400,
#                FPC 200 x 800): the clip check and the informativeness requirement (STOPPING gates),
#                DEFF_aug pooled and per class (REPORTED against its unmet 1.50 target), and the
#                store record (identity- and derived-tier SHA-256; exact size-rounding ties).
# --job panel    C-AE-3: kappa_AE derived ONCE by census algebra from the complete production panel
#                population set -- FS populations 1..160 at the panel cell (5, 80, 800, 400) -- and saved with its
#                provenance in runs/frozen/kappa_ae.json (an existing record is verified and reused, never
#                re-derived). G6 limb (a) for every stochastic control at the planned count: the panel (160 x 5)
#                and the U10 cell (12 x 15, dry run), with their Monte Carlo standard errors.
# --job g1g4     G1: C3a's host (FPC f = 0.50, populations 1-12, R = 10,000; binding) and the reported
#                FPC contrasts; G4 (reported): the single-linear-index R^2 of logit Q0(a, .) and logit g0 on the P0
#                grid at H = 5, 25, 50.
# --job summary  reads every prep output: completeness (every production population exactly once), the stopping
#                gates, the store manifest, DEFF_aug, G1, G4, G6 and the kappa_AE record; writes PREP-GATES.md and
#                prep_gates.json; exits nonzero iff a stopping gate fails or an output is missing, which keeps the
#                fitted cells from starting (a Slurm afterok dependency).
# No sample is drawn and no estimator runs. Only census quantities are computed; no confirmatory output (the
# census-only cells' production statistics, e.g. the Psi-bridge) is computed here.

via_wrapper <- endsWith(getOption("sim.engine.script", ""), "/runs/run_prep.R")
if (!via_wrapper) stop("run through the locked wrapper: Rscript simulation/engine/env/run.R simulation/engine/runs/run_prep.R --job ...")
now <- function() proc.time()[["elapsed"]]
t_start <- now()
utc_start <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
engine_dir <- getOption("sim.engine.dir")
repo_root <- normalizePath(file.path(engine_dir, "..", ".."), winslash = "/", mustWork = TRUE)
args <- getOption("sim.engine.args")
arg_of <- function(name, default) { i <- which(args == name); if (length(i) && length(args) > i[1]) args[i[1] + 1L] else default }
ENGINE_FILES <- c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R", "dgp.R", "census.R",
                  "quad_rules_io.R", "variance.R", "controls.R", "gates.R", "grading.R", "run_mode.R")
for (f in ENGINE_FILES) source(file.path(engine_dir, "R", f), local = TRUE)

JOB <- arg_of("--job", "")
if (!JOB %in% c("census", "panel", "g1g4", "summary")) stop("--job must be census, panel, g1g4 or summary", call. = FALSE)
spec <- spec_freeze_check(repo_root)
if (!isTRUE(spec$ok)) stop("the preparation runs only against the recorded design specification: SHA-256 ", spec$actual, " vs frozen ", spec$expected, call. = FALSE)
ROOT <- SIM_PROD_ROOT
co <- dgp_coefficients(1.25)                                       # a_cpsu = 1.25
N0 <- 20L; Q0C <- 0.005; DELTA <- 0.01
INFO_FLOOR <- 0.30; DEFF_TARGET <- 1.50                            # stopping floor; reported target (not met)
MASTERS <- list(MS = c(H = 800L, J = 20L), FS = c(H = 5L, J = 6400L), FPC = c(H = 200L, J = 800L))
out_root <- arg_of("--out-root", file.path(engine_dir, "runs", "out"))
pkg_commit <- if (file.exists(file.path(repo_root, "PACKAGE-COMMIT.txt"))) readLines(file.path(repo_root, "PACKAGE-COMMIT.txt"), warn = FALSE)[1] else NA_character_
sha <- function(p) if (file.exists(p)) digest::digest(file = p, algo = "sha256") else NA_character_
engine_hashes <- stats::setNames(lapply(file.path(engine_dir, "R", ENGINE_FILES), sha), paste0("R/", ENGINE_FILES))
peak_mb <- function() { f <- "/proc/self/status"; l <- if (file.exists(f)) grep("^VmHWM:", readLines(f, warn = FALSE), value = TRUE) else character(0)
  if (length(l)) as.numeric(gsub("[^0-9]", "", l[1])) / 1024 else NA_real_ }
w17 <- function(df, path) {
  for (cl in names(df)) if (is.double(df[[cl]])) df[[cl]] <- ifelse(is.na(df[[cl]]), NA_character_, sprintf("%.17g", df[[cl]]))
  utils::write.csv(df, path, row.names = FALSE, na = "NA")
}
finish <- function(out_dir, md, status, extra = list()) {
  jsonlite::write_json(c(list(job = JOB, created_utc = utc_start, wall_seconds = now() - t_start, peak_rss_mb = peak_mb(),
                              spec_freeze = spec, root = ROOT, package_commit = pkg_commit, engine_files_sha256 = engine_hashes,
                              runtime = list(R = R.version.string, platform = R.version$platform), exit_status = status), extra),
                       file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = I(17), null = "null", na = "null")
  writeLines(md, file.path(out_dir, "SUMMARY.md"), useBytes = TRUE)
  if (JOB == "summary") writeLines(md, file.path(out_dir, "PREP-GATES.md"), useBytes = TRUE)
  cat(sprintf("[prep] %s done in %.1f s; exit %d; wrote %s\n", JOB, now() - t_start, status, out_dir))
  if (status != 0L) quit(save = "no", status = status)
}
tag_of <- function(master, b) sprintf("svytmle-sim-store-v1;root=%s;master=%s;pop=%d", ROOT, master, b)

# ==== census: the production-population gates =========================================================
if (JOB == "census") {
  MASTER <- arg_of("--master", "")
  if (!MASTER %in% names(MASTERS)) stop("--master must be MS, FS or FPC", call. = FALSE)
  pops <- parse_ids(arg_of("--pops", ""))
  if (!all(pops %in% PROD_POPS[[MASTER]])) stop("--pops must be production population IDs of master ", MASTER, call. = FALSE)
  task <- as.integer(arg_of("--task", "0")); block <- as.integer(arg_of("--block", "0"))
  if (task > 0L) {
    first <- (task - 1L) * block + 1L
    if (first > length(pops)) stop("--task beyond the population list", call. = FALSE)
    pops <- pops[first:min(task * block, length(pops))]
  }
  out_dir <- file.path(out_root, sprintf("prep-census-%s%s-%s", MASTER, if (task > 0L) sprintf("-t%03d", task) else "", stamp))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  MH <- MASTERS[[MASTER]][["H"]]; MJ <- MASTERS[[MASTER]][["J"]]
  cat(sprintf("[prep] census: master %s (%d x %d), %d populations (%s)\n", MASTER, MH, MJ, length(pops), paste(range(pops), collapse = "..")))
  rows <- list()
  for (b in pops) {
    t0 <- now()
    pop <- build_population(MH, MJ, co, MASTER, b, ROOT, keep_identity = TRUE)
    t_build <- now() - t0
    hs <- population_store_hashes(pop, tag_of(MASTER, b))
    kb_h <- tapply(pop$psu$KB, pop$psu$h, `[`, 1L)
    info <- c(pooled = census_informativeness(pop),
              K1 = census_informativeness(pop_subset_strata(pop, as.integer(names(kb_h)[kb_h == 1L]))),
              K2 = census_informativeness(pop_subset_strata(pop, as.integer(names(kb_h)[kb_h == 2L]))))
    deff <- census_deff_aug(pop, m_h = 2L, n0 = N0)
    cl <- census_clip(pop$units, Q0C, DELTA)
    rows[[length(rows) + 1L]] <- data.frame(
      master = MASTER, pop_id = b, H = MH, J = MJ, N = nrow(pop$units), M_pop = nrow(pop$psu),
      clip_min_g = cl[["min_g"]], clip_max_g = cl[["max_g"]], clip_min_Q = cl[["min_Q"]], clip_max_Q = cl[["max_Q"]],
      clip_pass = cl[["pass"]] == 1,
      info_pooled = info[["pooled"]], info_K1 = info[["K1"]], info_K2 = info[["K2"]],
      info_pass = is.finite(info[["pooled"]]) && info[["pooled"]] >= INFO_FLOOR,
      deff_aug_pooled = deff[["pooled"]], deff_aug_K1 = deff[["K1"]], deff_aug_K2 = deff[["K2"]],
      deff_target_met = is.finite(deff[["pooled"]]) && deff[["pooled"]] >= DEFF_TARGET,
      size_ties = pop$meta$n_size_ties, identity_sha256 = hs[["identity"]], derived_sha256 = hs[["derived"]],
      build_seconds = t_build, seconds = now() - t0, stringsAsFactors = FALSE)
    cat(sprintf("[prep]   %s pop %d: N = %d, info %.4f (K1 %.4f, K2 %.4f), DEFF_aug %.4f (K1 %.4f, K2 %.4f), clip %s, ties %d, %.1f s, peak %.0f MB\n",
                MASTER, b, nrow(pop$units), info[["pooled"]], info[["K1"]], info[["K2"]], deff[["pooled"]], deff[["K1"]], deff[["K2"]],
                if (cl[["pass"]] == 1) "ok" else "FAIL", pop$meta$n_size_ties, now() - t0, peak_mb()))
    flush(stdout())
    rm(pop); invisible(gc())
  }
  df <- do.call(rbind, rows)
  w17(df, file.path(out_dir, "prep_census.csv"))
  md <- c(sprintf("# Preparation: production-population census gates, master %s (run %s)", MASTER, stamp), "",
    "> Spec v9.8 (frozen SHA-256 verified). Census algebra on each production population's whole master frame (A3 §A3.2). Stopping gates: clip, informativeness (>= 0.30). DEFF_aug is REPORTED against its original, unmet 1.50 target (R24), never PASS/FAIL.", "",
    sprintf("- %d populations (%s); clip failures %d; informativeness below 0.30: %d; min informativeness %.4f; DEFF_aug pooled %.4f-%.4f (target 1.50 met by %d); size ties %d.",
            nrow(df), paste(range(df$pop_id), collapse = ".."), sum(!df$clip_pass), sum(!df$info_pass), min(df$info_pooled),
            min(df$deff_aug_pooled), max(df$deff_aug_pooled), sum(df$deff_target_met), sum(df$size_ties)),
    sprintf("- Wall %.1f s; peak RSS %.0f MB.", now() - t_start, peak_mb()))
  finish(out_dir, md, 0L, list(master = MASTER, population_ids = pops, task = task, block = block))
}

# ==== panel: kappa_AE (C-AE-3, derived once) and G6 limb (a) ===========================================
if (JOB == "panel") {
  out_dir <- file.path(out_root, sprintf("prep-panel-%s", stamp))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  QR <- quad_rules_dev(engine_dir); quad_rules_activate(QR, required = TRUE)
  ctx <- p0_context(5L, co)
  gh <- QR$manifest$expected_grid_hashes; gi <- which(gh$H == 5L)
  grid_ok <- length(gi) == 1L && identical(ctx$primary_hash, gh$primary_sha256[gi])
  if (!grid_ok) stop("the H = 5 P0 grid does not reproduce grids/MANIFEST.json", call. = FALSE)
  C <- AE_PANEL_CELL
  cens <- list(); u10 <- list()
  for (b in AE_PANEL_POPS) {
    pop <- build_population(C$H, C$J, co, "FS", b, ROOT)
    cen <- census_targets(pop, m_h = C$m_h, n0 = N0)
    wc <- wrong_clustering_census(pop, C$m_h, N0)
    cens[[length(cens) + 1L]] <- data.frame(pop_id = b, sigma2_m = cen$sigma2_m, sigma2_m_wor = cen$sigma2_m_wor,
      Psi_N = cen$Psi_N, Psi_N_gf = cen$Psi_N_gf, R_var_nomh = wc$R_var[["nomh"]], R_var_unit = wc$R_var[["unit"]],
      R_var_nostrat = wc$R_var[["nostrat"]], stringsAsFactors = FALSE)
    if (b <= 12L) {                                   # the U10 cell: FS populations 1-12 at the panel coordinates
      uc <- u10_census(pop, C$m_h, N0)
      u10[[length(u10) + 1L]] <- data.frame(pop_id = b, Psi_U10 = uc$Psi_U10, sigma2_m_u10 = uc$sigma2_m_u10, Psi_N = uc$Psi_N,
                                            pred_mean_Z = uc$bias / sqrt(uc$sigma2_m_u10), stringsAsFactors = FALSE)
    }
    rm(pop)
  }
  cdf <- do.call(rbind, cens); udf <- do.call(rbind, u10)
  sw <- sqrt(cdf$sigma2_m_wor)
  fz <- file.path(engine_dir, AE_FROZEN_FILE)
  cb <- ae_clip_bound(co, Q0C, DELTA)
  if (file.exists(fz)) {
    rec <- ae_frozen_load(fz, spec_sha = spec$actual, grid_hash = ctx$primary_hash)
    same_sigma <- isTRUE(all.equal(as.numeric(rec$sigma_wor), sw, tolerance = 0))
    reused <- TRUE
    cat(sprintf("[prep] the frozen kappa_AE record exists and verifies: kappa_AE = %.10f (reused, never re-derived; census sigma_wor identical: %s)\n",
                rec$kappa, same_sigma))
  } else {
    ae <- ae_kappa_solve(ctx$primary_truth, ctx$primary$w, sw, cb$kappa_max)
    rec <- list(rule = AE_RULE, root = ROOT, master_id = "FS", population_ids = as.integer(AE_PANEL_POPS), cell = C,
                target = PANEL_KAPPA, kappa = if (ae$buildable) ae$kappa else NA_real_, T2 = if (ae$buildable) ae$T2 else NA_real_,
                stat = if (ae$buildable) ae$stat else NA_real_, buildable = ae$buildable, monotone = ae$monotone, at_max = ae$at_max,
                scale = ae$scale, kappa_max = ae$kappa_max, reason = ae$reason,
                clip = list(kappa_max = cb$kappa_max, min_g = cb$min_g, max_g = cb$max_g, min_Q = cb$min_Q, max_Q = cb$max_Q,
                            inside = cb$inside, q0 = cb$q0, delta = cb$delta),
                sigma_wor = sw, p0_primary_hash = ctx$primary_hash, quad_rules_sha256 = QR$sha256, spec_sha256 = spec$actual,
                package_commit = pkg_commit, engine_files_sha256 = engine_hashes[c("R/controls.R", "R/census.R", "R/dgp.R", "R/nuisance_truth.R", "R/remainders.R")],
                created_utc = utc_start, R = R.version.string, platform = R.version$platform)
    rec <- ae_frozen_write(fz, rec)
    rec <- ae_frozen_load(fz, spec_sha = spec$actual, grid_hash = ctx$primary_hash)      # read back: the saved value is the value
    reused <- FALSE; same_sigma <- TRUE
    cat(sprintf("[prep] kappa_AE derived ONCE from FS populations 1..160 and frozen: %s\n",
                if (isTRUE(rec$buildable)) sprintf("kappa_AE = %.10f, T2 = %.6e, Q_0.9 = %.6f", rec$kappa, rec$T2, rec$stat) else paste("NOT BUILDABLE:", rec$reason)))
  }
  file.copy(fz, file.path(out_dir, "kappa_ae.json"), overwrite = TRUE)
  # ---- G6 limb (a), boundary reading, MCSE at the planned count ----
  B <- length(AE_PANEL_POPS); Rp <- 5L
  shift_b <- -PANEL_KAPPA * sw / sqrt(cdf$sigma2_m)
  m_psu <- C$m; n_units <- C$m * N0; Hc <- C$H
  ae_vals <- if (isTRUE(rec$buildable)) rec$T2 / sw else rep(NA_real_, B)
  ae_mcse <- if (isTRUE(rec$buildable)) {
    idx <- boot_index(ROOT, C, "G6-AE-Q090", "prep-v9.8", B)
    stats::sd(boot_census(ae_vals, function(v) unname(stats::quantile(v, 0.9, type = 7)), idx)$boot)
  } else NA_real_
  g6 <- rbind(
    cbind(control = "C1 shift: mean Z keyed to Psi_N + 0.30 sigma_m,wor (R8)", statistic = "mean Z",
          g6_boundary(mean(shift_b), -0.10, 0.10, g6_mcse_joint(shift_b, 1, Rp))),
    cbind(control = "C2 wrong clustering: no m_h/(m_h-1), mean R_var", statistic = "mean R_var",
          g6_boundary(mean(cdf$R_var_nomh), 0.85, 1.15, g6_mcse_joint(cdf$R_var_nomh, mean(cdf$R_var_nomh)^2 * 2 / (m_psu - Hc), Rp))),
    cbind(control = "C2 wrong clustering: unit-level clusters, mean R_var", statistic = "mean R_var",
          g6_boundary(mean(cdf$R_var_unit), 0.85, 1.15, g6_mcse_joint(cdf$R_var_unit, mean(cdf$R_var_unit)^2 * 2 / (n_units - Hc), Rp))),
    cbind(control = "C2 wrong clustering: no stratification, mean R_var", statistic = "mean R_var",
          g6_boundary(mean(cdf$R_var_nostrat), 0.85, 1.15, g6_mcse_joint(cdf$R_var_nostrat, mean(cdf$R_var_nostrat)^2 * 2 / (m_psu - 1), Rp))),
    cbind(control = "C4 aligned error: Q_0.9(|T2| / sigma_m,wor) (kappa_AE frozen, C-AE-3)", statistic = "Q_0.9",
          g6_boundary(if (isTRUE(rec$buildable)) unname(stats::quantile(ae_vals, 0.9, type = 7)) else NA_real_, 0, 0.10, ae_mcse)))
  g6$cell <- "E3 power-control panel (160 x 5, R7)"; g6$planned <- "160 x 5"
  g6u <- cbind(control = sprintf("U10 (design-conditional, population %d): mean Z_U10", udf$pop_id), statistic = "mean Z",
               g6_boundary(udf$pred_mean_Z, -0.10, 0.10, g6_mcse_design(1, 15L)))
  g6u$cell <- "G8 U10 cell (dry run 12 x 15)"; g6u$planned <- "12 x 15"
  g6all <- rbind(g6, g6u)
  w17(cdf, file.path(out_dir, "panel_census.csv")); w17(udf, file.path(out_dir, "u10_census.csv")); w17(g6all, file.path(out_dir, "g6.csv"))
  role <- function(d) ifelse(d$inside, "diagnostic (inside its region, R5)", ifelse(d$graded, "graded power control", "implementation check (<= 3 MCSE)"))
  md <- c(sprintf("# Preparation: kappa_AE (C-AE-3) and G6 limb (a) (run %s)", stamp), "",
    "> Spec v9.8 (frozen SHA-256 verified). Census algebra on the production panel populations, FS 1..160 at (H, m_h, J, m) = (5, 80, 800, 400). kappa_AE is derived once from the complete set (R23 as amended) and reused by every later job; G6 limb (a) is read from the region boundary (R5) with the MCSE at the planned count (A3 §A3.8).", "",
    sprintf("- kappa_AE record: %s; %s.", if (reused) "EXISTING record verified and reused (never re-derived)" else "derived once now and frozen",
            if (isTRUE(rec$buildable)) sprintf("kappa_AE = %.10f in [0, kappa_max = %.6f]; T2 = %.6e; predicted Q_0.9(|T2|/sigma_m,wor) = %.6f (target 0.30); content sha256 %s",
                                               rec$kappa, rec$kappa_max, rec$T2, rec$stat, substr(rec$content_sha256, 1, 16)) else paste("NOT BUILDABLE --", rec$reason)),
    sprintf("- sigma_m,wor over the 160 populations: %.6f-%.6f (median %.6f).", min(sw), max(sw), stats::median(sw)), "",
    "| cell | control | predicted | region | distance from boundary | MCSE (planned count) | role |", "|---|---|---|---|---|---|---|",
    sprintf("| %s | %s | %.4f | [%.2f, %.2f] | %.4f | %.4f | %s |", g6all$cell, g6all$control, g6all$predicted, g6all$lo, g6all$hi,
            g6all$distance, g6all$mcse, role(g6all)))
  finish(out_dir, md, 0L,
         list(kappa_ae = rec[c("kappa", "T2", "stat", "buildable", "content_sha256")], reused = reused, sigma_identical_if_reused = same_sigma))
}

# ==== g1g4: G1 (C3a's host binding; FPC contrasts reported) and G4 (reported) ==========================
if (JOB == "g1g4") {
  out_dir <- file.path(out_root, sprintf("prep-g1g4-%s", stamp))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  g1 <- list()
  for (J in c(100L, 20L, 8L, 4L)) for (b in 1:12) {
    pop <- build_population(50L, J, co, "FPC", b, ROOT)
    cen <- census_targets(pop, m_h = 2L, n0 = N0); fs <- census_fpc_summary(cen)
    R_pl <- if (J == 4L) 10000L else 2000L
    mcse <- (cen$sigma2_m_wor / cen$sigma2_m) * sqrt(2 / (R_pl - 1))
    g1[[length(g1) + 1L]] <- data.frame(cell = sprintf("E2 FPC oracle J = %d (f = %s)", J, c("100" = "0.02", "20" = "0.10", "8" = "0.25", "4" = "0.50")[[as.character(J)]]),
      J = J, pop_id = b, R_planned = R_pl, contrast = fs$Delta_over_sigma2_m, c3a_target = 1 - fs$Delta_over_sigma2_m,
      mcse = mcse, clears_3mcse = fs$Delta_over_sigma2_m > 3 * mcse, binding = J == 4L, beta_bar = fs$beta_pooled,
      sigma_ratio = fs$sigma_ratio, overlay_080 = fs$overlay[["0.80"]], overlay_090 = fs$overlay[["0.90"]], overlay_095 = fs$overlay[["0.95"]],
      stringsAsFactors = FALSE)
  }
  g1 <- do.call(rbind, g1)
  QR <- quad_rules_dev(engine_dir); quad_rules_activate(QR, required = TRUE)
  g4 <- list()
  for (H in c(5L, 25L, 50L)) {
    ctx <- p0_context(H, co)
    gp <- ctx$primary; tr <- ctx$primary_truth
    X <- cbind(1, gp$C1, gp$C2, gp$C3, gp$C4, gp$Cpsu, as.numeric(gp$KB == 2L))
    for (nm in c("Q1", "Q0", "g")) {
      y <- stats::qlogis(tr[[nm]]); w <- gp$w / sum(gp$w)
      fit <- stats::lm.wfit(X, y, w)
      r2 <- 1 - sum(w * fit$residuals^2) / sum(w * (y - sum(w * y))^2)
      g4[[length(g4) + 1L]] <- data.frame(H = H, nuisance = c(Q1 = "logit Q0(1, .)", Q0 = "logit Q0(0, .)", g = "logit g0")[[nm]],
                                          R2_single_linear_index = r2, grid_points = nrow(gp), stringsAsFactors = FALSE)
    }
  }
  g4 <- do.call(rbind, g4)
  w17(g1, file.path(out_dir, "g1.csv")); w17(g4, file.path(out_dir, "g4.csv"))
  bind <- g1[g1$binding, ]
  md <- c(sprintf("# Preparation: G1 and G4 (run %s)", stamp), "",
    "> Spec v9.8 (frozen SHA-256 verified). G1 binds C3a's graded host (FPC f = 0.50, populations 1-12, R = 10,000; R9): its census FPC contrast Delta_m/sigma2_m must exceed 3 MCSE of Var(S_des)/sigma2_m at R (A3 §A3.12). The other fractions' contrasts are reported (displayed fractions; C3b's counts are budget-time choices). G4 is reported, no cutoff.", "",
    sprintf("- **G1 (binding, C3a host):** %d of %d populations clear 3 MCSE (contrast %.4f-%.4f against 3 MCSE %.4f-%.4f).",
            sum(bind$clears_3mcse), nrow(bind), min(bind$contrast), max(bind$contrast), 3 * min(bind$mcse), 3 * max(bind$mcse)), "",
    "| cell | binding | populations clearing 3 MCSE | contrast Delta_m/sigma2_m (min-max) | 3 MCSE (max) | C3a target 1 - Delta/sigma2 (min-max) | beta_bar (min-max) |",
    "|---|---|---|---|---|---|---|",
    vapply(split(g1, g1$J), function(d) sprintf("| %s | %s | %d/%d | %.4f-%.4f | %.4f | %.4f-%.4f | %.4f-%.4f |", d$cell[1], d$binding[1],
      sum(d$clears_3mcse), nrow(d), min(d$contrast), max(d$contrast), 3 * max(d$mcse), min(d$c3a_target), max(d$c3a_target),
      min(d$beta_bar), max(d$beta_bar)), ""), "",
    "**G4 (reported):** R^2 of the best single linear index (P0-weighted least squares on C1..C4, C_psu, 1{K_B = 2}); 1 - R^2 is the share the interactions carry.", "",
    "| H | nuisance | R^2 |", "|---|---|---|", sprintf("| %d | %s | %.4f |", g4$H, g4$nuisance, g4$R2_single_linear_index))
  finish(out_dir, md, 0L, list(g1_binding_clear = sum(bind$clears_3mcse), g1_binding_n = nrow(bind)))
}

# ==== summary: completeness, stopping gates, the store manifest ========================================
if (JOB == "summary") {
  in_root <- arg_of("--in", dirname(out_root))          # the round folder: runs/out/<round>
  out_dir <- file.path(out_root, sprintf("prep-summary-%s", stamp))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  newest <- function(pat) { d <- list.files(in_root, pattern = pat, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
    d <- d[dir.exists(d)]; if (length(d)) d[order(basename(d), decreasing = TRUE)] else character(0) }
  cfiles <- list.files(in_root, pattern = "^prep_census\\.csv$", full.names = TRUE, recursive = TRUE)
  cen <- if (length(cfiles)) do.call(rbind, lapply(cfiles, utils::read.csv, stringsAsFactors = FALSE)) else data.frame()
  problems <- character(0)
  want <- do.call(rbind, lapply(names(PROD_POPS), function(m) data.frame(master = m, pop_id = PROD_POPS[[m]], stringsAsFactors = FALSE)))
  if (nrow(cen)) {
    key <- paste(cen$master, cen$pop_id); dup <- key[duplicated(key)]
    # A rerun task may repeat populations; identical hashes are required, and one row per population is kept.
    if (length(dup)) {
      for (k in unique(dup)) { d <- cen[key == k, ]; if (length(unique(d$identity_sha256)) > 1L) problems <- c(problems, paste("population", k, "has two different identity hashes")) }
      cen <- cen[!duplicated(key), ]
    }
  }
  miss <- setdiff(paste(want$master, want$pop_id), paste(cen$master, cen$pop_id))
  if (length(miss)) problems <- c(problems, sprintf("%d production population(s) missing from the census outputs (first: %s)", length(miss), paste(head(miss, 3), collapse = ", ")))
  clip_fail <- if (nrow(cen)) sum(!as.logical(cen$clip_pass)) else NA
  info_fail <- if (nrow(cen)) sum(!as.logical(cen$info_pass)) else NA
  pdir <- newest("^prep-panel-"); gdir <- newest("^prep-g1g4-")
  fz <- file.path(engine_dir, AE_FROZEN_FILE)
  rec <- tryCatch(ae_frozen_load(fz, spec_sha = spec$actual), error = function(e) { problems <<- c(problems, conditionMessage(e)); NULL })
  if (!length(pdir)) problems <- c(problems, "no prep-panel output")
  if (!length(gdir)) problems <- c(problems, "no prep-g1g4 output")
  g6 <- if (length(pdir) && file.exists(file.path(pdir[1], "g6.csv"))) utils::read.csv(file.path(pdir[1], "g6.csv"), stringsAsFactors = FALSE) else data.frame()
  g1 <- if (length(gdir) && file.exists(file.path(gdir[1], "g1.csv"))) utils::read.csv(file.path(gdir[1], "g1.csv"), stringsAsFactors = FALSE) else data.frame()
  g4 <- if (length(gdir) && file.exists(file.path(gdir[1], "g4.csv"))) utils::read.csv(file.path(gdir[1], "g4.csv"), stringsAsFactors = FALSE) else data.frame()
  stop_ok <- !length(problems) && isTRUE(clip_fail == 0) && isTRUE(info_fail == 0)
  status <- if (stop_ok) 0L else 1L
  if (nrow(cen)) w17(cen[order(cen$master, cen$pop_id), c("master", "pop_id", "N", "identity_sha256", "derived_sha256", "size_ties")],
                     file.path(out_dir, "production_population_store.csv"))
  if (nrow(cen)) w17(cen, file.path(out_dir, "prep_census_all.csv"))
  by_master <- if (nrow(cen)) do.call(rbind, lapply(split(cen, cen$master), function(d) data.frame(master = d$master[1], n = nrow(d),
    clip_fail = sum(!as.logical(d$clip_pass)), info_min = min(d$info_pooled), info_K1_min = min(d$info_K1), info_K2_min = min(d$info_K2),
    info_fail = sum(!as.logical(d$info_pass)), deff_min = min(d$deff_aug_pooled), deff_med = stats::median(d$deff_aug_pooled),
    deff_max = max(d$deff_aug_pooled), deff_K1_med = stats::median(d$deff_aug_K1), deff_K2_med = stats::median(d$deff_aug_K2),
    deff_met = sum(as.logical(d$deff_target_met)), ties = sum(d$size_ties), sec_med = stats::median(d$seconds), stringsAsFactors = FALSE))) else data.frame()
  md <- c(sprintf("# Round-3 preparation gates (run %s)", stamp), "",
    sprintf("> Spec v9.8 (frozen SHA-256 %s verified). Stopping gates before any fitted cell of the dry run: the clip check and the informativeness requirement on every production population (R1; A3 §A3.2), every preparation output present, and the frozen kappa_AE record verified (R23). DEFF_aug is reported against its original, unmet 1.50 target (R24). G1, G4 and G6 are reported here and enter the execution sign-off.", substr(spec$actual, 1, 12)), "",
    sprintf("**Stopping gates: %s.** Populations: %d of %d present; clip failures %s; informativeness below 0.30: %s.%s",
            if (stop_ok) "PASS -- the dry run's cells may start" else "NOT PASSED -- the dry run's cells do not start", nrow(cen), nrow(want),
            clip_fail, info_fail, if (length(problems)) paste0(" Problems: ", paste(problems, collapse = "; ")) else ""), "",
    "| master | populations | clip failures | informativeness min (K1 / K2 min) | below 0.30 | DEFF_aug pooled min / median / max | per class median K1 / K2 | target 1.50 met | size ties | s per population |",
    "|---|---|---|---|---|---|---|---|---|---|",
    if (nrow(by_master)) sprintf("| %s | %d | %d | %.4f (%.4f / %.4f) | %d | %.4f / %.4f / %.4f | %.4f / %.4f | %d of %d (reported; the target is not met by the frozen DGP, R24) | %d | %.1f |",
      by_master$master, by_master$n, by_master$clip_fail, by_master$info_min, by_master$info_K1_min, by_master$info_K2_min, by_master$info_fail,
      by_master$deff_min, by_master$deff_med, by_master$deff_max, by_master$deff_K1_med, by_master$deff_K2_med, by_master$deff_met, by_master$n,
      by_master$ties, by_master$sec_med) else "| (none) |", "",
    sprintf("- **kappa_AE (C-AE-3):** %s", if (is.null(rec)) "MISSING or unverified" else if (isTRUE(rec$buildable))
      sprintf("frozen record verified: kappa_AE = %.10f, Q_0.9 = %.6f, content sha256 %s (FS populations 1..160)", rec$kappa, rec$stat, substr(rec$content_sha256, 1, 16)) else
      paste("NOT BUILDABLE (R10: recorded, the AE control is not run):", rec$reason)),
    if (nrow(g1)) sprintf("- **G1 (C3a host, binding):** %d of %d populations clear 3 MCSE at R = 10,000.", sum(g1$clears_3mcse[g1$binding]), sum(g1$binding)),
    if (nrow(g6)) c("- **G6 limb (a):**", sprintf("  - %s: predicted %.4f vs [%.2f, %.2f], distance %.4f, 3 MCSE %.4f -- %s", g6$control, g6$predicted, g6$lo, g6$hi,
                                                g6$distance, 3 * g6$mcse, g6$status)),
    if (nrow(g4)) sprintf("- **G4 (reported):** single-linear-index R^2 %s.", paste(sprintf("%s at H = %d: %.3f", g4$nuisance, g4$H, g4$R2_single_linear_index), collapse = "; ")))
  jsonlite::write_json(list(stopping_gates = if (stop_ok) "PASS" else "NOT PASSED", problems = problems, clip_failures = clip_fail,
                            informativeness_failures = info_fail, populations_present = nrow(cen), populations_expected = nrow(want),
                            kappa_ae = if (is.null(rec)) NULL else rec[c("kappa", "T2", "stat", "buildable", "content_sha256")],
                            kappa_ae_file_sha256 = sha(fz),
                            g1_binding = if (nrow(g1)) list(clear = sum(g1$clears_3mcse[g1$binding]), n = sum(g1$binding)) else NULL,
                            g6 = if (nrow(g6)) g6[, c("cell", "control", "predicted", "lo", "hi", "distance", "mcse", "graded", "inside", "status")] else NULL,
                            by_master = by_master),
                       file.path(out_dir, "prep_gates.json"), auto_unbox = TRUE, pretty = TRUE, digits = I(17), null = "null", na = "null")
  finish(out_dir, md, status, list(stopping_gates = if (stop_ok) "PASS" else "NOT PASSED"))
}
