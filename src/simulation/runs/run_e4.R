# E4 -- the model-layer-free deterministic array (claim C7: an
# IMPLEMENTATION CHECK, never PASS or RED-FLAG evidence).
#
#   Rscript src/simulation/env/run.R src/simulation/runs/run_e4.R [--R 2000] [--u8-first 20] [--out-root <dir>]
#
# Per rung H in {25, 50, 100, 200, 400} (m_h = 2, J = 20, f = 0.1, m = 2H):
#   census (free): the regularity gate on the regular array -- uniform bound C_t = 1.25, non-degeneracy
#     0.01 <= m sigma2_m,wor <= 100 (both REQUIRED), the Lindeberg ratio reported (a limit) -- and
#     G6 limb (b) on the dominant-stratum control: ratio >= delta_L = 0.50 (REQUIRED, every rung),
#     with its breach of the uniform bound stated (REQUIRED);
#   sampled oracle: R SRS-WOR draws of m_h = 2 of J = 20 PSUs per stratum through the engine's
#     sampler (one unit per PSU, n0 = 1, w = 10, N_hat = N exactly). REQUIRED: U1 on every
#     draw, N_hat == N bit-exact on every draw, U8 (survey::svymean, fpc not declared) on the first
#     --u8-first draws, finite outputs, ambient RNG. U2 is recorded as VACUOUS, never as a
#     pass. C1's oracle-limb margins on Z_oracle (mean, SD, Kolmogorov distance), coverage and the
#     first-stage UC/WOR variance behaviour are REPORTED (no grade).
# No model is posited or invoked: no DGP, no eta0, no learner, no P0 grid. E4 is not a two-stage
# variance check and is never graded on C3a.
#
# Run modes: --mode dev|dry|prod (R/run_mode.R) with --pops 1 (population id 1 names the single
# fixed array of each rung; the cell e2c label "e4" keeps E4's sample streams apart from the E1 many-strata
# cells at the same (H, m_h, J)); dry draws use rep IDs 1000001.., production 1..2000. REQUIRED:
# every rung's regular array and dominant-stratum control reproduce the frozen hashes of
# manifest/E4-ARRAY-HASHES.csv.
#
# Writes runs/out/dev-e4-<UTC>/: e4_census.csv, results.csv, checks.csv, manifest.json, SUMMARY.md.

via_wrapper <- endsWith(getOption("sim.engine.script", ""), "/runs/run_e4.R")
if (!via_wrapper) stop("run through the locked wrapper: Rscript simulation/engine/env/run.R simulation/engine/runs/run_e4.R")
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
ENGINE_FILES <- c("rng_streams.R", "census.R", "sampling.R", "variance.R", "e4.R", "run_mode.R")
for (f in ENGINE_FILES) source(file.path(engine_dir, "R", f), local = TRUE)

repo_root <- normalizePath(file.path(engine_dir, "..", ".."), winslash = "/", mustWork = TRUE)
RM <- run_mode_setup(arg_of, repo_root, B_default = 1L, R_default = 2000L)
MODE <- RM$mode; ROOT <- RM$root; REP_IDS <- RM$reps; R <- RM$R; RUN_LABEL <- RM$label
if (!identical(RM$pops, 1L)) stop("E4 has one fixed array per rung: --pops 1", call. = FALSE)
U8_FIRST <- as.integer(arg_of("--u8-first", "20"))
if (is.na(R) || R < 2L) stop("--R must be >= 2", call. = FALSE)
E4_FROZEN <- utils::read.csv(file.path(engine_dir, "manifest", "E4-ARRAY-HASHES.csv"), stringsAsFactors = FALSE)
if (is.na(U8_FIRST) || U8_FIRST < 1L) stop("--u8-first must be >= 1", call. = FALSE)
out_dir <- file.path(arg_of("--out-root", file.path(engine_dir, "runs", "out")), paste0(MODE, "-e4-", stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("[e4] %s  E4 deterministic array: H = %s, m_h = %d, J = %d, R = %d draws per rung  out=%s\n",
            utc_start, paste(E4_H, collapse = "/"), E4_MH, E4_J, R, out_dir))

checks <- list(); cens <- list(); res <- list(); T_PH <- c(census = 0, sampling = 0)
add <- function(H, check, value, reference, pass, required = TRUE, status = NULL)
  checks[[length(checks) + 1L]] <<- data.frame(H = H, check = check, value = value, reference = reference,
    pass = pass, required = required,
    status = if (!is.null(status)) status else if (!required) "REPORT" else if (isTRUE(pass)) "PASS" else "FAIL",
    stringsAsFactors = FALSE)

for (H in E4_H) {
  t_c <- now()
  pop <- e4_population(e4_array(H))
  pop_dom <- e4_population(e4_array(H, dominant = TRUE))
  reg <- e4_census(pop); dom <- e4_census(pop_dom)
  g <- e4_gate_rows(H, reg, dom)
  for (i in seq_len(nrow(g))) add(H, g$check[i], g$value[i], g$reference[i], g$pass[i], g$required[i])
  for (a in list(list(h = pop$meta$hash, lab = "regular"), list(h = pop_dom$meta$hash, lab = "dominant-stratum control"))) {
    fz <- E4_FROZEN$array_sha256[E4_FROZEN$H == H & E4_FROZEN$array == a$lab]
    add(H, sprintf("E4 %s array reproduces its frozen hash (A3 §A3.11; %s)", a$lab, if (length(fz) == 1L) substr(fz, 1, 12) else "NO RECORD"),
        NA_real_, NA_real_, length(fz) == 1L && identical(a$h, fz))
  }
  add(H, "E4 Lindeberg ratio x H (report; ~constant when the ratio falls like 1/H)", reg$lindeberg * H, NA_real_, TRUE, FALSE)
  for (a in list(list(x = reg, lab = "regular"), list(x = dom, lab = "dominant-stratum control")))
    cens[[length(cens) + 1L]] <- data.frame(H = H, array = a$lab, N = a$x$N, m = a$x$m, Psi_N = a$x$Psi_N,
      sigma2_m = a$x$sigma2_m, sigma2_m_wor = a$x$sigma2_m_wor, Delta_m = a$x$Delta_m, lindeberg = a$x$lindeberg,
      lindeberg_times_H = a$x$lindeberg * H, m_sigma2_wor = a$x$m_sigma2_wor, max_abs_t = a$x$max_abs_t,
      array_sha256 = if (a$lab == "regular") pop$meta$hash else pop_dom$meta$hash, stringsAsFactors = FALSE)
  T_PH[["census"]] <- T_PH[["census"]] + (now() - t_c)

  # Sampled oracle draws on the regular array (one unit per PSU).
  t_s <- now()
  cell <- list(regime = "many", H = H, m_h = E4_MH, J = E4_J, m = E4_MH * H, f = "1/10", e2c = "e4")
  frame <- sampling_frame(pop, H, E4_J, psu_cols = character(0))
  u1_fail <- 0L; nhat_bad <- 0L; rng_bad <- 0L; fin_bad <- 0L; u8_fail <- 0L
  psi <- se2 <- numeric(R)
  for (ri in seq_len(R)) {
    r <- REP_IDS[ri]
    rng_r <- ambient_rng()
    s <- draw_sample(frame, E4_MH, 1L, list(cell = cell, population_id = 1L, rep_id = r), ROOT, psu_cols = character(0))
    u1 <- check_u1(s, pop)
    u1_fail <- u1_fail + !u1$pass
    nhat_bad <- nhat_bad + !identical(sum(s$w), as.numeric(reg$N))
    uc <- hajek_uc(s$phi0, s$w, s$h, s$k)
    psi[ri] <- uc$psi; se2[ri] <- uc$se2
    fin_bad <- fin_bad + !all(is.finite(c(uc$psi, uc$se2))) + !(uc$se2 > 0)
    if (ri <= U8_FIRST) {
      u8 <- check_u8(s$phi0, s$w, s$h, s$k)
      u8_fail <- u8_fail + !u8$pass
    }
    rng_bad <- rng_bad + !identical(rng_r, ambient_rng())
  }
  Z <- (psi - reg$Psi_N) / sqrt(se2)
  res[[length(res) + 1L]] <- data.frame(H = H, rep_id = REP_IDS, psi_hat = psi, se2_hat = se2, Z = Z)
  add(H, sprintf("U1 (bit-exact weights, design counts) held on all %d draws (failures in value)", R), u1_fail, 0, u1_fail == 0L)
  add(H, sprintf("R11: N_hat == N bit-exact on all %d draws (w = M_pop,h/m_h = 10; draws off in value)", R), nhat_bad, 0, nhat_bad == 0L)
  add(H, sprintf("U8 se2 = survey::svymean, fpc not declared (rel 1e-12) on the first %d draws (failures in value)", min(R, U8_FIRST)),
      u8_fail, 0, u8_fail == 0L)
  add(H, sprintf("finite outputs: oracle psi, se2 > 0 on all %d draws (failures in value)", R), fin_bad, 0, fin_bad == 0L)
  add(H, sprintf("ambient RNG state and kind unchanged by every draw (%d draws; moved in value)", R), rng_bad, 0, rng_bad == 0L)
  add(H, "U2: N_hat/N - 1 is identically 0 under R11's unit layer -- VACUOUS, never a pass", 0, 0, NA, required = FALSE, status = "VACUOUS")
  ks <- ks_normal(Z)
  rep_rows <- data.frame(
    check = c(sprintf("C1 oracle limb (report): mean Z_oracle over %d draws", R),
              sprintf("C1 oracle limb (report): SD Z_oracle over %d draws", R),
              sprintf("C1 oracle limb (report): Kolmogorov distance of Z_oracle from N(0,1) over %d draws", R),
              "coverage of the 95% interval (report)",
              "first stage (report): mean UC se2 / sigma2_m",
              "first stage (report): Var(psi_hat) / sigma2_m,wor (the WOR design variance)"),
    value = c(mean(Z), stats::sd(Z), ks, mean(abs(Z) <= stats::qnorm(0.975)), mean(se2) / reg$sigma2_m,
              stats::var(psi) / reg$sigma2_m_wor),
    reference = c(0, 1, 0, 0.95, 1, 1))
  for (i in seq_len(nrow(rep_rows))) add(H, rep_rows$check[i], rep_rows$value[i], rep_rows$reference[i], TRUE, FALSE)
  T_PH[["sampling"]] <- T_PH[["sampling"]] + (now() - t_s)
  cat(sprintf("[e4] H = %3d: L = %.4f (control %.4f), m sigma2_wor = %.4f; %d draws: mean Z %+.3f, SD Z %.3f, KS %.4f (%.1f s)\n",
              H, reg$lindeberg, dom$lindeberg, reg$m_sigma2_wor, R, mean(Z), stats::sd(Z), ks, now() - t_c))
  flush(stdout())
}
rng_ok <- identical(rng_start, ambient_rng())
add(NA, "ambient RNG state and kind unchanged by the whole run", as.numeric(rng_ok), 1, rng_ok)

checks_df <- do.call(rbind, checks); cens_df <- do.call(rbind, cens); res_df <- do.call(rbind, res)
w17 <- function(df, path) {
  for (cl in names(df)) if (is.double(df[[cl]])) df[[cl]] <- ifelse(is.na(df[[cl]]), NA_character_, sprintf("%.17g", df[[cl]]))
  utils::write.csv(df, path, row.names = FALSE, na = "NA")
}
w17(cens_df, file.path(out_dir, "e4_census.csv")); w17(res_df, file.path(out_dir, "results.csv"))
w17(checks_df, file.path(out_dir, "checks.csv"))
req <- checks_df$required
n_req <- sum(req); n_fail <- sum(req & !(checks_df$pass %in% TRUE))
status <- if (n_fail) 1L else 0L
peak_kb <- local({ f <- "/proc/self/status"; l <- if (file.exists(f)) grep("^VmHWM:", readLines(f, warn = FALSE), value = TRUE) else character(0)
  if (length(l)) as.numeric(gsub("[^0-9]", "", l[1])) else NA_real_ })
wall <- now() - t_start
T_PH <- c(T_PH, other = wall - sum(T_PH), wall = wall)
sha <- function(p) if (file.exists(p)) digest::digest(file = p, algo = "sha256") else NA_character_
repo_root <- normalizePath(file.path(engine_dir, "..", ".."), winslash = "/", mustWork = TRUE)
jsonlite::write_json(list(
  label = paste(RUN_LABEL, "-- E4 deterministic array (C7 implementation check)"), mode = MODE, spec_freeze = RM$spec,
  created_utc = utc_start, wall_seconds = wall, phase_seconds = as.list(T_PH), memory = list(peak_rss_kb = peak_kb),
  package_commit = if (file.exists(file.path(repo_root, "PACKAGE-COMMIT.txt"))) readLines(file.path(repo_root, "PACKAGE-COMMIT.txt"), warn = FALSE)[1] else NA,
  config = list(root = ROOT, H = E4_H, J = E4_J, m_h = E4_MH, n0 = 1L, R = R, rep_ids = c(first = REP_IDS[1], last = REP_IDS[R]), u8_first = U8_FIRST,
                recipe = "§12.7: t_hk = A s(h) cos(2 pi k/J) + A a(h), A = 1; control s(1) = sqrt(H)",
                constants = list(C_t = E4_CT, delta_L = E4_DELTA_L, cond_nd = E4_ND_C),
                unit_layer = "R11 (DECIDED 2026-10-11): M_hk = 1, w = M_pop,h/m_h, N_hat = N; U2 VACUOUS",
                array_sha256 = stats::setNames(as.list(cens_df$array_sha256), paste(cens_df$array, cens_df$H))),
  keys = list(sample_example = rng_key("sample", list(root = ROOT, cell = list(regime = "many", H = 25L, m_h = 2L, J = 20L, m = 50L,
                                                                                f = "1/10", e2c = "e4"), population_id = 1L, rep_id = 1L))),
  engine_files_sha256 = c(stats::setNames(lapply(file.path(engine_dir, "R", ENGINE_FILES), sha), paste0("R/", ENGINE_FILES)),
                          list("runs/run_e4.R" = sha(getOption("sim.engine.script")))),
  runtime = list(R = R.version.string, platform = R.version$platform),
  checks = list(n = nrow(checks_df), n_required = n_req, n_required_failed = n_fail, exit_status = status)),
  file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = I(17), null = "null", na = "null")

rep_v <- function(H, pat) { i <- checks_df$H %in% H & grepl(pat, checks_df$check, fixed = TRUE); if (any(i)) checks_df$value[i][1] else NA }
md <- c(sprintf("# E4 deterministic array (%s), run %s", toupper(MODE), stamp), "",
  sprintf("> **%s.**", RUN_LABEL),
  "> **An implementation check (C7), never PASS or RED-FLAG evidence about theorem truth.** No model is posited or invoked: the census PSU totals are written down directly (§12.7). R11's unit layer: one unit per PSU, so N_hat = N exactly and U2 is VACUOUS. Every Z / coverage / variance number is reported, not graded.", "",
  sprintf("- Rungs H = %s, m_h = %d, J = %d (f = 0.1, m = 2H); %d SRS-WOR draws per rung; U8 on the first %d draws.", paste(E4_H, collapse = ", "), E4_MH, E4_J, R, min(R, U8_FIRST)),
  sprintf("- Wall time %.1f s (census %.1f, sampling and oracle %.1f, other %.1f); peak RSS %s MB.", wall, T_PH[["census"]], T_PH[["sampling"]],
          T_PH[["other"]], if (is.na(peak_kb)) "NA" else sprintf("%.0f", peak_kb / 1024)), "",
  sprintf("**Required checks: %d of %d passed. Exit status %d.** U2: VACUOUS at every rung (R11), not counted as a pass.", n_req - n_fail, n_req, status), "",
  "## Regularity gate and the dominant-stratum control (census algebra)", "",
  "| H | m sigma2_m,wor | max abs t | Lindeberg ratio | ratio x H | control ratio (>= 0.50) | control max abs t (> 1.25) |",
  "|---|---|---|---|---|---|---|")
reg_df <- cens_df[cens_df$array == "regular", ]; dom_df <- cens_df[cens_df$array != "regular", ]
md <- c(md, sprintf("| %d | %.4f | %.4f | %.5f | %.3f | %.4f | %.2f |", reg_df$H, reg_df$m_sigma2_wor, reg_df$max_abs_t,
                    reg_df$lindeberg, reg_df$lindeberg_times_H, dom_df$lindeberg, dom_df$max_abs_t), "",
        "## Sampled oracle (report only)", "",
        "| H | mean Z | SD Z | Kolmogorov | cover95 | mean se2 / sigma2_m | Var(psi) / sigma2_m,wor |", "|---|---|---|---|---|---|---|",
        vapply(E4_H, function(H) sprintf("| %d | %+.4f | %.4f | %.4f | %.4f | %.4f | %.4f |", H,
          rep_v(H, "mean Z_oracle"), rep_v(H, "SD Z_oracle"), rep_v(H, "Kolmogorov"), rep_v(H, "coverage of the 95%"),
          rep_v(H, "mean UC se2"), rep_v(H, "Var(psi_hat)")), ""))
if (n_fail) md <- c(md, "", "Failed required checks:", "", sprintf("- H %s: %s", checks_df$H, checks_df$check)[req & !(checks_df$pass %in% TRUE)])
writeLines(md, file.path(out_dir, "SUMMARY.md"), useBytes = TRUE)
cat(sprintf("[e4] done in %.1f s: %d checks, %d required, %d failed; exit %d\n[e4] wrote %s\n", wall, nrow(checks_df), n_req, n_fail, status, out_dir))
if (status != 0L) quit(save = "no", status = status)
