# CENSUS LIMBS of the E2 constructions: the strict trajectory, the vanishing-fraction cell and the
# unequal-allocation pair.
#
#   Rscript src/simulation/env/run.R src/simulation/runs/run_e2_census.R [--traj strict|vanish|unequal] [--B 2] [--R 5] [--out-root <dir>]
#
# --traj strict (default): the E2 genuinely-strict trajectory's census limb
#   (the strictness witness is census algebra, not fitted simulation), H in {25, 50,
#   100, 200, 400}, m_h = 2, J = 4, f = 0.5, the first H strata and the first 4 PSUs of the
#   many-strata master [forced]. Delta_m/sigma2_m reported, never graded here.
# --traj vanish: the E2 vanishing-fraction cell, the ONLY Psi-bridge PASS
#   cell (census limb 120 x 0). (H, J, m_h) in {(25, 10, 2), (50, 20, 2),
#   (100, 40, 2), (200, 80, 2)}, f = 0.20 ... 0.025 (the reported trajectory), plus the bridge-only
#   terminal rung (200, 800, 2), f = 0.0025, all prefixes of master FPC (200 x 800).
#   Per rung: Psi of P_{0,H} on the saved-rule P0 grid (hash = grids/MANIFEST.json; U9's P0-grid limb),
#   and the bridge |Psi_N - Psi| / sigma_m,wor -- graded only at the terminal rung (200, 800) in
#   production (Q_0.9 <= 0.10); reported here. Memory: the full 200 x 800 frame is built
#   once per population and serves every rung; the terminal rung's nesting is checked on a sub-block
#   (the first NEST_SUB_H strata, all 800 PSUs) instead of a second full frame; the job requests 16 GiB.
# --traj unequal: the E2 unequal-allocation invariance pair:
#   H = 50, J = 20, prefix of the FPC master (its 200 x 100 block); the SAME population under
#   "high-on-high" (K_B = 2 strata m_h = 6, K_B = 1 m_h = 2) and "high-on-low" (swapped), m = 200 in
#   both. Census limb: Delta_m equal to floating point (REQUIRED, an implementation check), the
#   contrast in sigma2_m reported, G1's "the intended variance contrast exists" reported. Plus a
#   small sampled oracle limb (B x R draws per allocation; class-indexed m_h in the sample-stream
#   key) exercising the sampler, U1 and U8 on both allocations.
#
# Run modes: --mode dev|dry|prod (R/run_mode.R) and --pops a:b. dry: these are
# census-only cells (R = 0), so the dry run uses the reserved population IDs 1000001-1000005 and computes no
# confirmatory census output; prod: production populations (1..120). The unequal pair's sampled
# oracle limb is a dry-run implementation check (rep IDs 1000001..); production runs its census limb only.
# Populations: strict builds the (400, 4) block of master MS once per population (every rung is a prefix);
# unequal the (50, 20) block of master FPC; vanish the full 200 x 800 frame. Nesting is checked against a
# direct build (a sub-block at the terminal rung). No fitted estimator, no grade.
#
# Writes runs/out/dev-e2e-<strict-census|vanish-census|unequal-pair>-<UTC>/: census_trajectory.csv
# (or census_pair.csv and results.csv), census_strata.csv, checks.csv, manifest.json, SUMMARY.md.
# Required checks: master-prefix nesting, census FPC algebra, G2 feasibility, finite targets, the
# ambient RNG; vanish adds the exact f per rung, the P0 grid hashes and U9's P0-grid limb; unequal
# adds the Delta_m invariance, equal total m, and U1/U8/finite outputs on every sampled draw.

via_wrapper <- endsWith(getOption("sim.engine.script", ""), "/runs/run_e2_census.R")
if (!via_wrapper) stop("run through the locked wrapper: Rscript simulation/engine/env/run.R simulation/engine/runs/run_e2_census.R")
now <- function() proc.time()[["elapsed"]]
t_start <- now()
utc_start <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
genv <- globalenv()
ambient_rng <- function() list(kind = RNGkind(),
  seed = if (exists(".Random.seed", envir = genv, inherits = FALSE)) get(".Random.seed", envir = genv) else NULL)
rng_start <- ambient_rng()
engine_dir <- getOption("sim.engine.dir")
args <- getOption("sim.engine.args")
arg_of <- function(name, default) { i <- which(args == name); if (length(i) && length(args) > i[1]) args[i[1] + 1L] else default }
ENGINE_FILES <- c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R", "dgp.R", "census.R",
                  "quad_rules_io.R", "folds.R", "sampling.R", "variance.R", "estimators.R", "run_mode.R")
for (f in ENGINE_FILES) source(file.path(engine_dir, "R", f), local = TRUE)

A_CPSU <- 1.25                                  # frozen
N0 <- 20L
TRAJ <- arg_of("--traj", "strict")
if (!TRAJ %in% c("strict", "vanish", "unequal")) stop("--traj must be strict, vanish or unequal", call. = FALSE)
if (TRAJ == "strict") {
  MASTER_ID <- "MS"; MASTER_H <- 400L; MASTER_J <- 4L                           # the (400, 4) block of MS 800 x 20
  RUNGS <- data.frame(H = c(25L, 50L, 100L, 200L, 400L), J = 4L)
} else if (TRAJ == "vanish") {
  MASTER_ID <- "FPC"; MASTER_H <- 200L; MASTER_J <- 800L                        # J = 800 for the bridge-only rung
  RUNGS <- data.frame(H = c(25L, 50L, 100L, 200L, 200L), J = c(10L, 20L, 40L, 80L, 800L))
} else {
  # The 200 x 100 block of master FPC (200 x 800): entity addressing makes its prefixes identical
  # to the full master's, so only the block a cell needs is built (the unequal pair: H = 50, J = 20).
  MASTER_ID <- "FPC"; MASTER_H <- 50L; MASTER_J <- 20L                          # the (50, 20) block of FPC 200 x 800
  RUNGS <- data.frame(H = 50L, J = 20L)
}
TERMINAL <- c(H = 200L, J = 800L)        # vanish: the bridge-only terminal rung
NEST_SUB_H <- 10L                        # vanish: strata in the terminal rung's sub-block nesting check (A2)
MASTER_LABEL <- c(strict = "the (400, 4) block of master MS (800 x 20)", vanish = "master FPC (200 x 800, Amendment A2)",
                  unequal = "the (50, 20) block of master FPC (200 x 800)")[[TRAJ]]
M_H <- 2L
repo_root <- normalizePath(file.path(engine_dir, "..", ".."), winslash = "/", mustWork = TRUE)
RM <- run_mode_setup(arg_of, repo_root, B_default = 2L, R_default = 5L, census_only = TRUE, master_id = MASTER_ID)
MODE <- RM$mode; ROOT <- RM$root; POP_IDS <- RM$pops; B <- RM$B; RUN_LABEL <- RM$label
R <- as.integer(arg_of("--R", "5"))
# The unequal pair's sampled oracle limb: a dry-run (and dev) implementation check; production is census only.
REP_IDS <- if (TRAJ != "unequal" || MODE == "prod") integer(0) else if (MODE == "dry") DRY_ID_BASE + seq_len(R) else seq_len(R)
if (TRAJ == "unequal" && MODE != "prod" && (is.na(R) || R < 2L)) stop("--R must be >= 2 for the sampled oracle limb", call. = FALSE)
TAG <- c(strict = "strict-census", vanish = "vanish-census", unequal = "unequal-pair")[[TRAJ]]
if (RM$task > 0L) TAG <- sprintf("%s-t%03d", TAG, RM$task)
out_dir <- file.path(arg_of("--out-root", file.path(engine_dir, "runs", "out")), paste0(MODE, "-e2e-", TAG, "-", stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
co <- dgp_coefficients(A_CPSU)
cat(sprintf("[census] %s  mode %s  E2 %s: rungs (H, J) = %s, m_h = %d, B = %d (pops %s)%s  out=%s\n", utc_start, MODE, TAG,
            paste(sprintf("(%d,%d)", RUNGS$H, RUNGS$J), collapse = " "), M_H, B, paste(unique(range(POP_IDS)), collapse = ".."),
            if (length(REP_IDS)) sprintf(", R = %d per allocation (sampled oracle limb)", length(REP_IDS)) else "", out_dir))

checks <- list(); traj <- list(); strata <- list(); pairs <- list(); results <- list()
add <- function(pop, H, check, value, reference, pass, required = TRUE, alloc = NA_character_)
  checks[[length(checks) + 1L]] <<- data.frame(pop_id = pop, H = H, alloc = alloc, check = check, value = value,
                                               reference = reference, pass = isTRUE(pass), required = required,
                                               stringsAsFactors = FALSE)
T_PH <- c(p0 = 0, populations = 0, census = 0)   # census includes the unequal pair's sampled limb

# ---- vanish: Psi of P_{0,H} per rung on the saved-rule grid (U9 P0-grid limb; hashes) ----------
PSI_H <- list(); GRID_HASH <- list()
if (TRAJ == "vanish") {
  t0 <- now()
  QR <- quad_rules_dev(engine_dir)
  quad_rules_activate(QR, required = TRUE)
  ctx0 <- p0_context(RUNGS$H[1], co)
  gh <- QR$manifest$expected_grid_hashes
  # Each rung's grids are built FRESH (p0_grid at that H), so their hashes are comparable with the
  # MANIFEST records: reweighting another H's grid is equal in algebra but not bit for bit. eta0 at
  # the nodes does not depend on H (only the K_B class weights do), so ctx0's truth is reused after
  # asserting the nodes are identical.
  fresh_ctx <- function(H) {
    if (H == RUNGS$H[1]) return(ctx0)
    d1 <- attr(ctx0$primary, "definition"); d2 <- attr(ctx0$second, "definition")
    gp <- p0_grid(H, d1$rule, d1$n_c, d1$n_psu); gs <- p0_grid(H, d2$rule, d2$n_c, d2$n_psu)
    nodes <- c("C1", "C2", "C3", "C4", "Cpsu", "KB")
    if (!identical(gp[nodes], ctx0$primary[nodes]) || !identical(gs[nodes], ctx0$second[nodes]))
      stop("P0 grid nodes differ across H; the truth cannot be reused", call. = FALSE)
    ctx <- ctx0; ctx$H <- H; ctx$primary <- gp; ctx$second <- gs
    ctx$primary_hash <- p0_grid_hash(gp); ctx$second_hash <- p0_grid_hash(gs)
    ctx
  }
  for (H in RUNGS$H) {
    ctx <- fresh_ctx(H)
    gi <- if (is.data.frame(gh)) which(gh$H == H) else integer(0)
    exp_hash <- if (length(gi) == 1L) c(gh$primary_sha256[gi], gh$second_sha256[gi]) else c(NA_character_, NA_character_)
    for (j in 1:2) add(NA, H, sprintf("P0 grid hash (%s, H = %d) = grids/MANIFEST.json record (%s)", c("primary", "second")[j], H,
                                      ifelse(is.na(exp_hash[j]), "NO RECORD", substr(exp_hash[j], 1, 12))),
                       NA_real_, NA_real_, identical(c(ctx$primary_hash, ctx$second_hash)[j], exp_hash[j]))
    u9 <- u9_psi(ctx)
    for (i in seq_len(nrow(u9))) add(NA, H, paste("U9 P0-grid limb:", u9$check[i]), u9$value[i], u9$reference[i], u9$pass[i])
    PSI_H[[as.character(H)]] <- sum(ctx$primary$w * (ctx$primary_truth$Q1 - ctx$primary_truth$Q0))
    GRID_HASH[[as.character(H)]] <- ctx$primary_hash
  }
  rm(ctx0, ctx)
  T_PH[["p0"]] <- now() - t0
  cat(sprintf("[census] P0 grids: Psi(P0,H) = %s (%.1f s)\n", paste(sprintf("%.8f", unlist(PSI_H)), collapse = " / "), T_PH[["p0"]]))
}

census_row <- function(b, H, J, cen, fs, pop) data.frame(
  pop_id = b, H = H, J = J, m_h = M_H, m = H * M_H, f = fs$f_bar, fpc = fs$fpc_bar, N = cen$N, M_pop = cen$M_pop,
  Psi_N = cen$Psi_N, sigma2_m = cen$sigma2_m, sigma2_m_wor = cen$sigma2_m_wor, Delta_m = cen$Delta_m,
  Delta_over_sigma2_m = fs$Delta_over_sigma2_m, sigma_ratio = fs$sigma_ratio, beta_pooled = fs$beta_pooled,
  mean_beta_h = mean(fs$per_stratum$beta_h), single_stage_fpc_deficit = fs$single_stage_fpc_deficit,
  c3b_overlay_080 = fs$overlay[["0.80"]], c3b_overlay_090 = fs$overlay[["0.90"]], c3b_overlay_095 = fs$overlay[["0.95"]],
  min_M = min(pop$psu$M), max_M = max(pop$psu$M), stringsAsFactors = FALSE)

common_checks <- function(b, H, J, cen, fs, pop, alloc = NA_character_) {
  tol <- 64 * .Machine$double.eps * cen$sigma2_m
  add(b, H, "census FPC algebra: sigma2_m - sigma2_m_wor = Delta_m = sum_h f_h B_h (64 eps sigma2_m)",
      cen$sigma2_m - cen$sigma2_m_wor, cen$Delta_m, fs$identity_gap <= tol, alloc = alloc)
  add(b, H, "G2 feasibility: n0/M_k >= 0.40 for every PSU (min n0/M_k in value)", N0 / max(pop$psu$M), 0.40,
      N0 / max(pop$psu$M) >= 0.40, alloc = alloc)
  fin <- all(is.finite(c(cen$Psi_N, cen$sigma2_m, cen$sigma2_m_wor, cen$Delta_m)))
  add(b, H, "finite census targets (Psi_N, sigma2_m, sigma2_m_wor, Delta_m)", as.numeric(fin), 1, fin, alloc = alloc)
}

for (b in POP_IDS) {
  t_b <- now()
  master <- build_population(MASTER_H, MASTER_J, co, MASTER_ID, b, ROOT)
  T_PH[["populations"]] <- T_PH[["populations"]] + (now() - t_b)
  for (ri in seq_len(nrow(RUNGS))) {
    H <- RUNGS$H[ri]; J <- RUNGS$J[ri]
    t_p <- now()
    full <- H == MASTER_H && J == MASTER_J
    pop <- if (full) master else population_prefix(master, H, J)
    if (full && TRAJ == "vanish") {
      # Never a second full 200 x 800 frame. The sub-block (first NEST_SUB_H strata, all J PSUs)
      # of the frame must equal its direct build -- the same entity-addressing identity, at the full J.
      ok <- identical(population_prefix(pop, NEST_SUB_H, J), build_population(NEST_SUB_H, J, co, MASTER_ID, b, ROOT))
      nest_lab <- sprintf("prefix view (%d,%d) = %s: sub-block (%d,%d) identical() to direct build (nesting; Amendment A2)",
                          H, J, MASTER_LABEL, NEST_SUB_H, J)
    } else {
      ok <- identical(pop, build_population(H, J, co, MASTER_ID, b, ROOT))
      nest_lab <- sprintf("prefix view (%d,%d) of %s identical() to direct build (nesting)", H, J, MASTER_LABEL)
      if (full) nest_lab <- sprintf("the block (%d,%d) of %s identical() to a direct build (nesting)", H, J, MASTER_LABEL)
    }
    T_PH[["populations"]] <- T_PH[["populations"]] + (now() - t_p)
    add(b, H, nest_lab, as.numeric(ok), 1, ok)
    t_c <- now()
    if (TRAJ != "unequal") {
      cen <- census_targets(pop, m_h = M_H, n0 = N0)
      fs <- census_fpc_summary(cen)
      common_checks(b, H, J, cen, fs, pop)
      row <- census_row(b, H, J, cen, fs, pop)
      if (TRAJ == "strict") {
        add(b, H, "strictness witness (report): Delta_m / sigma2_m > 0 at f = 0.5", fs$Delta_over_sigma2_m, 0,
            fs$Delta_over_sigma2_m > 0, required = FALSE)
      } else {
        f_exact <- all(cen$per_stratum$f_h == M_H / J)
        add(b, H, sprintf("vanishing rung: f_h = m_h/J = %d/%d exactly in every stratum", M_H, J), cen$per_stratum$f_h[1], M_H / J, f_exact)
        Psi <- PSI_H[[as.character(H)]]
        row$Psi <- Psi; row$Psi_N_minus_Psi <- cen$Psi_N - Psi
        row$bridge <- census_bridge(cen, Psi)
        row$terminal <- H == TERMINAL[["H"]] && J == TERMINAL[["J"]]
        add(b, H, sprintf("Psi-bridge (report): |Psi_N - Psi| / sigma_m,wor%s", if (row$terminal) " at the TERMINAL rung (200, 800) (graded after production by the population-only percentile bootstrap of Q_0.9, R6; margin 0.10, Amendment A2)" else sprintf(" at (%d, %d), the reported trajectory", H, J)),
            row$bridge, 0.10, is.finite(row$bridge), required = FALSE)
      }
      traj[[length(traj) + 1L]] <- row
      strata[[length(strata) + 1L]] <- cbind(data.frame(pop_id = b, H = H, KB = kb_of_stratum(fs$per_stratum$h)), fs$per_stratum)
    } else {
      up <- census_unequal_pair(pop, 6L, 2L, N0)
      for (al in c("high_on_high", "high_on_low")) {
        cen <- up[[al]]; fs <- census_fpc_summary(cen)
        common_checks(b, H, J, cen, fs, pop, alloc = al)
        r <- census_row(b, H, J, cen, fs, pop); r$alloc <- al
        r$m_high_class <- if (al == "high_on_high") 6L else 2L
        r$m <- cen$m
        pairs[[length(pairs) + 1L]] <- r
        strata[[length(strata) + 1L]] <- cbind(data.frame(pop_id = b, H = H, alloc = al, KB = up$KB), fs$per_stratum)
      }
      add(b, H, "unequal-allocation invariance: |Delta_m(high-on-high) - Delta_m(high-on-low)| <= 64 eps Delta_m (implementation check, §12.4)",
          up$delta_gap, 0, up$delta_gap <= 64 * .Machine$double.eps * up$high_on_high$Delta_m)
      add(b, H, sprintf("unequal-allocation pair: identical total m in both allocations (m = %d)", up$high_on_high$m),
          up$high_on_low$m, up$high_on_high$m, up$m_equal)
      add(b, H, "unequal-allocation pair: one population, run twice (identical Psi_N and N)", as.numeric(identical(up$high_on_high$Psi_N, up$high_on_low$Psi_N)),
          1, identical(up$high_on_high$Psi_N, up$high_on_low$Psi_N) && up$high_on_high$N == up$high_on_low$N)
      add(b, H, "G1 (report): the intended variance contrast exists -- mean census S2_t,h of the K_B = 2 strata exceeds K_B = 1's",
          up$contrast[["S2t_K2"]], up$contrast[["S2t_K1"]], up$contrast[["S2t_K2"]] > up$contrast[["S2t_K1"]], required = FALSE)
      # Sampled oracle limb: the same population under both allocations (distinct cell keys).
      frame <- sampling_frame(pop, H, J)
      for (al in c("high_on_high", "high_on_low")) {
        mh_cls <- if (al == "high_on_high") c(K1 = 2L, K2 = 6L) else c(K1 = 6L, K2 = 2L)
        f_cls <- if (al == "high_on_high") c(K1 = "1/10", K2 = "3/10") else c(K1 = "3/10", K2 = "1/10")
        cell <- list(regime = "many", H = H, m_h = mh_cls, J = J, m = 200L, f = f_cls, e2c = "unequal")
        cen <- up[[al]]
        for (r in REP_IDS) {
          rng_r <- ambient_rng()
          s <- draw_sample(frame, mh_cls, N0, list(cell = cell, population_id = b, rep_id = r), ROOT)
          u1 <- check_u1(s, pop)
          for (i in seq_len(nrow(u1$checks))) add(b, H, paste("U1", u1$checks$check[i]), as.numeric(u1$checks$pass[i]), 1, u1$checks$pass[i], alloc = al)
          or <- oracle_aipw(s)
          u8 <- check_u8(s$phi0, s$w, s$h, s$k)
          add(b, H, "U8 se2 = survey::svymean, fpc not declared (rel 1e-12)", u8$se2, u8$se2_survey, u8$pass, alloc = al)
          mh_real <- tapply(paste(s$h, s$k), s$h, function(x) length(unique(x)))
          want <- ifelse(up$KB == 2L, mh_cls[["K2"]], mh_cls[["K1"]])
          add(b, H, "sampled allocation matches the designated class allocation in every stratum",
              as.numeric(identical(as.integer(mh_real), as.integer(want))), 1, identical(as.integer(mh_real), as.integer(want)), alloc = al)
          Z <- (or$psi - cen$Psi_N) / or$se
          fin <- all(is.finite(c(or$psi, or$se, Z)))
          add(b, H, "finite outputs: oracle psi, se, Z", as.numeric(fin), 1, fin, alloc = al)
          rng_ok_r <- identical(rng_r, ambient_rng())
          add(b, H, "ambient RNG state and kind unchanged by the draw", as.numeric(rng_ok_r), 1, rng_ok_r, alloc = al)
          results[[length(results) + 1L]] <- data.frame(pop_id = b, rep_id = r, alloc = al, arm = "oracle", psi_hat = or$psi,
            se_hat = or$se, Psi_N = cen$Psi_N, sigma_m = sqrt(cen$sigma2_m), sigma_m_wor = sqrt(cen$sigma2_m_wor), Z = Z,
            R_var = or$se^2 / cen$sigma2_m, N_hat_over_N_minus_1 = sum(s$w) / cen$N - 1, n_units = nrow(s),
            n_psu = sum(attr(s, "design")$m_h), sample_key = attr(s, "design")$sample_key, stringsAsFactors = FALSE)
        }
      }
    }
    T_PH[["census"]] <- T_PH[["census"]] + (now() - t_c)
  }
  rm(master)
  cat(sprintf("[census] population %d: %d rung(s) in %.1f s\n", b, nrow(RUNGS), now() - t_b))
}
if (TRAJ == "vanish") quad_rules_deactivate()
rng_ok <- identical(rng_start, ambient_rng())
add(NA, NA, "ambient RNG state and kind unchanged by the whole run", as.numeric(rng_ok), 1, rng_ok)

checks_df <- do.call(rbind, checks)
traj_df <- if (length(traj)) do.call(rbind, traj) else data.frame()
pairs_df <- if (length(pairs)) do.call(rbind, pairs) else data.frame()
results_df <- if (length(results)) do.call(rbind, results) else data.frame()
strata_df <- if (length(strata)) do.call(rbind, strata) else data.frame()
w17 <- function(df, path) {
  for (cl in names(df)) if (is.double(df[[cl]])) df[[cl]] <- ifelse(is.na(df[[cl]]), NA_character_, sprintf("%.17g", df[[cl]]))
  utils::write.csv(df, path, row.names = FALSE, na = "NA")
}
if (nrow(traj_df)) w17(traj_df, file.path(out_dir, "census_trajectory.csv"))
if (nrow(pairs_df)) w17(pairs_df, file.path(out_dir, "census_pair.csv"))
if (nrow(results_df)) w17(results_df, file.path(out_dir, "results.csv"))
w17(strata_df, file.path(out_dir, "census_strata.csv"))
w17(checks_df, file.path(out_dir, "checks.csv"))
n_req <- sum(checks_df$required); n_fail <- sum(checks_df$required & !checks_df$pass)
status <- if (n_fail) 1L else 0L
peak_kb <- local({ f <- "/proc/self/status"; l <- if (file.exists(f)) grep("^VmHWM:", readLines(f, warn = FALSE), value = TRUE) else character(0)
  if (length(l)) as.numeric(gsub("[^0-9]", "", l[1])) else NA_real_ })
wall <- now() - t_start
T_PH <- c(T_PH, other = wall - sum(T_PH), wall = wall)
sha <- function(p) if (file.exists(p)) digest::digest(file = p, algo = "sha256") else NA_character_
repo_root <- normalizePath(file.path(engine_dir, "..", ".."), winslash = "/", mustWork = TRUE)
jsonlite::write_json(list(
  label = sprintf("%s -- E2 %s", RUN_LABEL, TAG), mode = MODE, spec_freeze = RM$spec,
  created_utc = utc_start, wall_seconds = wall, phase_seconds = as.list(T_PH), memory = list(peak_rss_kb = peak_kb),
  package_commit = if (file.exists(file.path(repo_root, "PACKAGE-COMMIT.txt"))) readLines(file.path(repo_root, "PACKAGE-COMMIT.txt"), warn = FALSE)[1] else NA,
  config = list(traj = TRAJ, root = ROOT, master_id = MASTER_ID, master_dims = c(H = MASTER_H, J = MASTER_J), a_cpsu = A_CPSU,
                rungs = RUNGS, m_h = M_H, n0 = N0, B = B, R = length(REP_IDS), population_ids = POP_IDS, rep_ids = REP_IDS,
                task = RM$task, block = RM$block,
                rulings = c(strict = "none (round 1)", vanish = "R12 (rungs), FPC master", unequal = "R13 (FPC master); allocation: author decision 2026-09-16")[[TRAJ]],
                psi_by_H = if (TRAJ == "vanish") PSI_H else NULL, grid_primary_hash = if (TRAJ == "vanish") GRID_HASH else NULL),
  keys = list(population = vapply(POP_IDS, function(b) rng_key("population", list(root = ROOT, master_id = MASTER_ID, population_id = b)), ""),
              sample = if (nrow(results_df)) unique(results_df$sample_key) else NULL),
  engine_files_sha256 = c(stats::setNames(lapply(file.path(engine_dir, "R", ENGINE_FILES), sha), paste0("R/", ENGINE_FILES)),
                          list("runs/run_e2_census.R" = sha(getOption("sim.engine.script")))),
  runtime = list(R = R.version.string, platform = R.version$platform),
  checks = list(n = nrow(checks_df), n_required = n_req, n_required_failed = n_fail, exit_status = status)),
  file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = I(17), null = "null", na = "null")

hdr <- c(strict = "E2 genuinely-strict trajectory: census limb", vanish = "E2 vanishing-fraction cell: census limb and the Psi-bridge",
         unequal = "E2 unequal-allocation invariance pair")[[TRAJ]]
md <- c(sprintf("# %s (%s), run %s", hdr, toupper(MODE), stamp), "",
  sprintf("> **%s.** Census algebra (spec %s, frozen SHA-256 verified: %s; adopted RNG contract, root `%s`; a_cpsu = 1.25, R1). Nothing here is graded.",
          RUN_LABEL, RM$spec$version, RM$spec$ok, ROOT), "",
  sprintf("- Rungs (H, J): %s; m_h = %d; prefixes of %s, population IDs %s%s.",
          paste(sprintf("(%d, %d)", RUNGS$H, RUNGS$J), collapse = ", "), M_H, MASTER_LABEL, paste(unique(range(POP_IDS)), collapse = ".."),
          if (length(REP_IDS)) sprintf("; sampled oracle limb: %d draws per allocation and population", length(REP_IDS)) else ""),
  sprintf("- Wall time %.1f s (P0 grids %.1f, populations %.1f, census and sampling %.1f, other %.1f); peak RSS %s MB.", wall,
          T_PH[["p0"]], T_PH[["populations"]], T_PH[["census"]], T_PH[["other"]], if (is.na(peak_kb)) "NA" else sprintf("%.0f", peak_kb / 1024)), "",
  sprintf("**Required checks: %d of %d passed. Exit status %d.**", n_req - n_fail, n_req, status))
fam <- sub("^(U[0-9]+ P0-grid limb|U1|U8|P0 grid hash|prefix view|the block|census FPC algebra|G2 feasibility|finite|ambient RNG|vanishing rung|Psi-bridge|strictness witness|unequal-allocation invariance|unequal-allocation pair|G1|sampled allocation).*$", "\\1", checks_df$check)
tab <- do.call(rbind, lapply(split(seq_len(nrow(checks_df)), list(fam, checks_df$required), drop = TRUE), function(i)
  data.frame(family = fam[i[1]], required = checks_df$required[i[1]], rows = length(i), passed = sum(checks_df$pass[i]))))
tab <- tab[order(!tab$required, tab$family), ]
md <- c(md, "", "| check family | required | rows | passed |", "|---|---|---|---|",
        sprintf("| %s | %s | %d | %d |", tab$family, ifelse(tab$required, "yes", "no (report)"), tab$rows, tab$passed))
if (TRAJ == "strict") md <- c(md, "",
  "| pop | H | m | N | Delta_m / sigma2_m | sigma_m / sigma_m_wor | pooled beta | mean beta_h | overlay 0.95 |",
  "|---|---|---|---|---|---|---|---|---|",
  sprintf("| %d | %d | %d | %d | %.4f | %.4f | %.4f | %.4f | %.4f |", traj_df$pop_id, traj_df$H, traj_df$m, traj_df$N,
          traj_df$Delta_over_sigma2_m, traj_df$sigma_ratio, traj_df$beta_pooled, traj_df$mean_beta_h, traj_df$c3b_overlay_095))
if (TRAJ == "vanish" && nrow(traj_df)) {
  term <- traj_df[traj_df$terminal, , drop = FALSE]
  md <- c(md, "", "## Psi-bridge along the vanishing-fraction rungs (report)", "",
    "> The bridge |Psi_N - Psi| / sigma_m,wor is graded only at the bridge-only terminal rung (200, 800) (Amendment A2, R22) and only in production, as Q_0.9 over the cell's populations against the unchanged 0.10 (§12.4, secondary Psi); the four R12 rungs are the reported trajectory. Psi is the generating-law parameter of P_{0,H} on the saved-rule P0 grid. Leading-order sizing forecast at the terminal rung (conditional, not a probability implied by a correct theorem): Q_0.9 ~ 0.07.", "",
    "| pop | H | J | f | m | N | Psi_N | Psi | Psi_N - Psi | sigma_m,wor | bridge | Delta_m / sigma2_m |",
    "|---|---|---|---|---|---|---|---|---|---|---|---|",
    sprintf("| %d | %d | %d | %.3f | %d | %d | %.6f | %.6f | %+.2e | %.6f | %.4f | %.4f |", traj_df$pop_id, traj_df$H, traj_df$J, traj_df$f,
            traj_df$m, traj_df$N, traj_df$Psi_N, traj_df$Psi, traj_df$Psi_N_minus_Psi, sqrt(traj_df$sigma2_m_wor), traj_df$bridge,
            traj_df$Delta_over_sigma2_m), "",
    sprintf("At the terminal rung over these %d populations: max bridge %.4f; Q_0.9 (type 7) %.4f. Descriptive only (a dry run computes no confirmatory output; production grades it, R6).",
            nrow(term), max(term$bridge), unname(stats::quantile(term$bridge, 0.9, type = 7))))
}
if (TRAJ == "unequal" && nrow(pairs_df)) {
  md <- c(md, "", "## Census limb: the same population under both allocations", "",
    "> Delta_m = N^-2 sum_h M_pop,h S2_t,h does not depend on the allocation (an implementation check at floating-point tolerance); the contrast lives in sigma2_m. \"High\" is the frozen class label K_B = 2, designated before any realised B_h is seen.", "",
    "| pop | allocation | m_h (K_B = 2 / 1) | m | sigma2_m | sigma2_m,wor | Delta_m | Delta_m / sigma2_m |",
    "|---|---|---|---|---|---|---|---|",
    sprintf("| %d | %s | %d / %d | %d | %.6e | %.6e | %.10e | %.4f |", pairs_df$pop_id, pairs_df$alloc, pairs_df$m_high_class,
            ifelse(pairs_df$m_high_class == 6L, 2L, 6L), pairs_df$m, pairs_df$sigma2_m, pairs_df$sigma2_m_wor, pairs_df$Delta_m,
            pairs_df$Delta_over_sigma2_m))
  if (nrow(results_df)) {
    per <- do.call(rbind, lapply(split(results_df, list(results_df$alloc, results_df$pop_id), drop = TRUE), function(d)
      data.frame(pop = d$pop_id[1], alloc = d$alloc[1], n = nrow(d), mean_Z = mean(d$Z), R_var = mean(d$R_var),
                 n_units = d$n_units[1], stringsAsFactors = FALSE)))
    md <- c(md, "", sprintf("Sampled oracle limb (%d draws per allocation and population; descriptive):", length(REP_IDS)), "",
            "| pop | allocation | draws | units per draw | mean Z | mean R_var |", "|---|---|---|---|---|---|",
            sprintf("| %d | %s | %d | %d | %+.3f | %.3f |", per$pop, per$alloc, per$n, per$n_units, per$mean_Z, per$R_var))
  }
}
if (n_fail) md <- c(md, "", "Failed required checks:", "", sprintf("- pop %s, H %s, %s: %s", checks_df$pop_id, checks_df$H,
                                                                 checks_df$alloc, checks_df$check)[checks_df$required & !checks_df$pass])
writeLines(md, file.path(out_dir, "SUMMARY.md"), useBytes = TRUE)
cat(sprintf("[census] done in %.1f s: %d checks, %d required, %d failed; exit %d\n[census] wrote %s\n", wall, nrow(checks_df), n_req, n_fail, status, out_dir))
if (status != 0L) quit(save = "no", status = status)
