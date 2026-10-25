# The grading pass over saved PRODUCTION output. One function,
# grade_run(), called by runs/run_grade.R and by tests/testthat/test-grading.R on synthetic trees.
# Requires R/rng_streams.R, R/run_mode.R, R/learner_rd.R (RD_ALPHAS, rd_arm), R/grading.R, R/grade_plan.R,
# R/grade_io.R. Post-processing only.
#
# What is graded, where (terminal rungs; regime scope):
#   C1   Z_CF = (psi_hat - Psi_N)/sigma_m,wor (column Z_wor), three margins, composed per cell:
#        many-strata Ra-AIPW at m = 1600 and Rb-AIPW at m = 800; the E2 fitted cells (FPC f = 0.5, strict
#        m = 200), Ra-AIPW and Rb-AIPW. R-b needs its quarter-rate verdict.
#        Implementation check (population-indexed): Z_oracle at the oracle terminal rung m = 1600, both regimes.
#   C2   many-strata Ra-AIPW (m = 1600) and Rb-AIPW (m = 800): fitted mean, fitted SD, paired transfer (R-b:
#        its L2 verdict). Implementation check: the oracle limb's mean and SD per population (many-strata m = 1600).
#   C3a  FPC f = 0.50 oracle, populations 1-12 x 10,000: signed Var(S_des)/sigma2_m - (1 - Delta_m/sigma2_m).
#   C3b  each FPC oracle cell, per population and level 0.80/0.90/0.95: coverage of oracle +- z sigma_m minus
#        the level (f = 0.50: the first 2000 redraws), only where that population's own Z_oracle passes all
#        three C1 margins.
#   C4a  Q_0.9(|T1|/sigma_m,wor): fixed-strata Ra-AIPW at m = 1600; each R-d alpha at m = 3200 (G5 must PASS).
#   C4b  Q_0.9(|T2|/sigma_m,wor) at m = 3200: alpha .30, .35 graded; .25 BOUNDARY; .20 STRESS (reported only).
#        alpha .20 and .35 were STOPPED after G5 -- UNGRADABLE (premise), not run in production.
#   C5   PASS iff, at a premise-violating rung alpha < 1/4 (.25 is BOUNDARY), G5 established
#        r_m -> 0 and C4a PASSes there; never RED-FLAG. The only such rung (.20) is stopped: UNGRADABLE.
#   A remainder grade (C4a, C4b) on an arm with an UNRESOLVED T2 in that cell is
#        UNGRADABLE (unresolved T2), the survivors' statistic descriptive; no value is substituted, no replicate dropped.
#   Ra-TMLE's T1 and T2 enter no grade and no display (e1_fixed_strata_by_arm carries a note
#        instead); they stay archived in the runs' results.csv.
#   C6   many-strata m = 800, Ra-TMLE and Rb-TMLE: Q_0.9 of |psi_TMLE - psi_AIPW|/sigma_m,wor and of
#        |se2_TMLE - se2_AIPW|/sigma2_m; only if the matched AIPW cell's C1 and C2 PASS, R-b's
#        quarter-rate holds and lambda_min >= c_H in every retained replicate.
#   sec. Psi  Q_0.9(|Psi_N - Psi|/sigma_m,wor) at the bridge-only terminal rung (200, 800), 120 populations.
#   C7   E4, per rung: Z_oracle margins (implementation check) and E4's required rows.
#   Power controls (E3 panel): graded only where G6 limb (a) said so; predicted to fail.
# Everything else is a display (trajectories, envelopes, the C3b overlay, the bridge trajectory, the
# Delta_m ~ Var(S_mod) relation, the strictness witness, the E3 foil and V-sensitivity tables).
#
# Every graded quantity passes ONE gate, in this order (grade_stat):
#   1. NOT GRADED  -- the cell is incomplete against the frozen design, or refused (design or provenance mismatch);
#   2. BLOCKED     -- a failed or missing required check affecting it, or a non-finite required value in a
#                     replicate marked OK;
#   3. UNGRADABLE (premise: ...) -- a premise the statistic needs was not established;
#   4. UNGRADABLE (fitting failures: ...) -- an arm it uses has a FAILED fit in that cell (any folder, superseded
#                     included); the survivors' statistic and interval are reported, labelled descriptive;
#      UNGRADABLE (unresolved T2: ...) -- a remainder grade whose arm has a replicate with an
#                     UNRESOLVED T2; the survivors' statistic and interval are reported, labelled descriptive;
#   5. the decision operator on the 95% interval.

GRADE_TMLE_C_H <- 0.01                  # c_H

# ---- preparation and execution-check inputs the grades depend on -------------------------------------------
newest_file <- function(root, pattern) {
  f <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(f)) f[order(basename(dirname(f)), decreasing = TRUE)][1] else NA_character_
}
# assembly_roots: further folders whose assemble/ may hold the execution-check summary (results/simulation/execution-check
# in this repository); the preparation outputs are under round3_root. The newest execution_gates.json wins.
read_round3 <- function(round3_root, assembly_roots = character(0)) {
  ex <- newest_file(file.path(c(assembly_roots, round3_root), "assemble"), "^execution_gates\\.json$")
  adir <- if (is.na(ex)) NA_character_ else dirname(ex)
  rd <- function(d, f) if (!is.na(d) && file.exists(file.path(d, f))) grade_read_csv(file.path(d, f)) else data.frame()
  g6f <- newest_file(file.path(round3_root, "prep-panel"), "^g6\\.csv$")
  pj <- newest_file(file.path(round3_root, "prep-summary"), "^prep_gates\\.json$")
  kj <- newest_file(file.path(round3_root, "prep-panel"), "^kappa_ae\\.json$")
  sha <- function(p) if (!is.na(p) && file.exists(p)) digest::digest(file = p, algo = "sha256") else NA_character_
  paths <- c(ex, if (is.na(adir)) rep(NA_character_, 3) else file.path(adir, c("g5.csv", "rb_gates.csv", "g7.csv")), g6f, pj, kj)
  list(exec = if (is.na(ex)) NULL else jsonlite::fromJSON(ex, simplifyVector = TRUE),
       g5 = rd(adir, "g5.csv"), rb = rd(adir, "rb_gates.csv"), g7 = rd(adir, "g7.csv"),
       g6 = if (is.na(g6f)) data.frame() else grade_read_csv(g6f),
       prep = if (is.na(pj)) NULL else jsonlite::fromJSON(pj, simplifyVector = TRUE),
       kappa = if (is.na(kj)) NULL else jsonlite::fromJSON(kj, simplifyVector = TRUE),
       files = data.frame(input = c("execution_gates.json", "g5.csv", "rb_gates.csv", "g7.csv", "g6.csv", "prep_gates.json", "kappa_ae.json"),
                          path = paths, sha256 = unname(vapply(paths, sha, "")), stringsAsFactors = FALSE))
}

# The nested rate verdicts for many-strata R-b; missing = unresolved.
rb_premise <- function(rb, level) {
  r <- if (nrow(rb)) rb[rb$regime %in% "many", , drop = FALSE] else rb
  if (!nrow(r) || !is.finite(r$upper_95[1])) return(list(ok = FALSE, why = "R-b's rate gate is missing from the round-3 assembly (unresolved)"))
  u <- r$upper_95[1]
  if (level == "L2") list(ok = u < 0, why = sprintf("R-b L2 gate: one-sided 95%% upper bound on the slope %.4f %s 0", u, if (u < 0) "<" else ">="))
  else list(ok = u < -0.25, why = sprintf("R-b quarter-rate gate: upper bound %.4f %s -1/4", u, if (u < -0.25) "<" else ">="))
}
g5_premise <- function(g5, alpha) {
  r <- if (nrow(g5)) g5[abs(g5$alpha - alpha) < 1e-9, , drop = FALSE] else g5
  if (!nrow(r)) return(list(ok = FALSE, why = sprintf("G5 missing for alpha %.2f (unresolved)", alpha)))
  ok <- isTRUE(as.logical(r$pass[1]))
  list(ok = ok, why = sprintf("G5 alpha %.2f: %s", alpha, if (ok) "PASS (r_m -> 0 at the engineered rate)" else "not PASS: the trajectory is UNGRADABLE and STOPPED (§1.7a)"))
}

# rbind two data frames over the union of their columns (a display row may lack a column another arm has).
rbind_cols_fill <- function(a, b) {
  if (is.null(a) || !nrow(a)) return(b)
  cols <- union(names(a), names(b))
  for (cl in setdiff(cols, names(a))) a[[cl]] <- NA
  for (cl in setdiff(cols, names(b))) b[[cl]] <- NA
  rbind(a[, cols, drop = FALSE], b[, cols, drop = FALSE])
}

frac_str <- function(a, b) { g <- gcd2(as.integer(a), as.integer(b)); sprintf("%d/%d", as.integer(a) %/% g, as.integer(b) %/% g) }
stat_label <- function(...) gsub("[^A-Za-z0-9._+-]", "", paste(..., sep = "."))

# Population-equal summaries of one arm's rows (displays).
arm_summary <- function(d, cols = c("Z_wor", "R_var")) {
  if (!nrow(d)) return(data.frame(rows = 0L, populations = 0L))
  f <- pop_freq(d$pop_id)
  out <- data.frame(rows = nrow(d), populations = length(unique(d$pop_id)))
  if ("Z_wor" %in% names(d) && "Z_wor" %in% cols) {
    out$mean_Z <- f_mean(d$Z_wor, f); out$sd_Z <- f_sd(d$Z_wor, f); out$ks_Z <- f_ks(d$Z_wor, f)
  }
  if ("cover95" %in% names(d)) out$cover95_studentized <- f_mean(as.numeric(as.logical(d$cover95)), f)
  if ("R_var" %in% names(d) && "R_var" %in% cols) { out$mean_R_var <- f_mean(d$R_var, f); out$sd_R_var <- f_sd(d$R_var, f) }
  out
}


# Displays carry their basis: every replicate, or the survivors of FAILED fits (descriptive).
basis_of <- function(failed) if (is.finite(failed) && failed > 0) sprintf("survivors only (%d FAILED; descriptive, R25)", as.integer(failed)) else "all replicates"

# Reproduction release: output made by the release carries the release's identity (config/execution-identity.json,
# written to PACKAGE-COMMIT.txt by its setup) instead of a commit. For such output the production guard holds when
# (a) every production run folder records this release's identity, (b) every engine file a run hashed has the SHA-256
# of this release's file, and (c) each of those files has the executable code of the file at the commit whose
# execution check passed -- the same bytes, or the same code with only comments edited (code_sha256) -- except the
# files the release lists as adapted. Output made at a commit keeps the original check.
release_root_default <- function() {
  e <- getOption("sim.engine.dir")
  if (is.null(e)) NULL else normalizePath(file.path(e, "..", ".."), winslash = "/", mustWork = FALSE)
}

# SHA-256 of an R file's executable code: the file with every comment removed (comments as R's parser finds them),
# trailing blanks removed and blank lines dropped. Comment-only edits keep it; any change to code or strings does not.
code_sha256 <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  pd <- utils::getParseData(parse(path, keep.source = TRUE, encoding = "UTF-8"), includeText = TRUE)
  cm <- pd[pd$token == "COMMENT", c("line1", "text")]
  for (k in seq_len(nrow(cm))) {
    l <- lines[cm$line1[k]]
    if (!endsWith(l, cm$text[k])) stop(path, ": comment not at the end of line ", cm$line1[k], call. = FALSE)
    lines[cm$line1[k]] <- substr(l, 1L, nchar(l) - nchar(cm$text[k]))
  }
  keep <- sub("[[:space:]]+$", "", lines)
  digest::digest(paste0(paste(keep[nzchar(keep)], collapse = "\n"), "\n"), algo = "sha256", serialize = FALSE)
}

release_identity_check <- function(commits, engine_files, release_root, r3_commit) {
  out <- list(applies = FALSE, ok = FALSE, check = "", detail = "")
  if (!length(commits) || !any(startsWith(commits, "release:"))) return(out)      # output made at a commit
  out$applies <- TRUE
  out$check <- "reproduction: every production run folder ran this release's code, which has the code the execution check PASSed except the release's adapted files"
  f <- if (is.null(release_root)) NA_character_ else file.path(release_root, "config", "execution-identity.json")
  if (is.na(f) || !file.exists(f)) { out$detail <- "no release record (config/execution-identity.json)"; return(out) }
  rec <- jsonlite::fromJSON(f, simplifyVector = FALSE)
  id <- paste0("release:", rec$release_id)
  ec <- rec$execution_check
  files <- unlist(ec$engine_files); code <- unlist(ec$engine_code); adapted <- as.character(unlist(ec$adapted))
  eng <- file.path(release_root, "src", "simulation")
  sha <- function(p) { q <- file.path(eng, p); if (file.exists(q)) digest::digest(file = q, algo = "sha256") else "missing" }
  now <- vapply(names(files), sha, "")
  same <- now == files
  code_now <- vapply(names(files), function(p) {
    q <- file.path(eng, p)
    if (!is.null(code) && p %in% names(code) && file.exists(q)) code_sha256(q) else NA_character_
  }, "")
  ref_code <- if (is.null(code)) rep(NA_character_, length(files)) else unname(code[names(files)])
  same_code <- !same & !is.na(code_now) & !is.na(ref_code) & code_now == ref_code
  prob <- character(0)
  if (!(length(commits) == 1L && identical(commits, id)))
    prob <- c(prob, sprintf("package identity %s, this release is %s", paste(commits, collapse = ","), id))
  if (!identical(ec$commit, r3_commit))
    prob <- c(prob, sprintf("execution-check commit %s in the release record, %s in the execution-check summary", ec$commit %||% "missing", r3_commit))
  bad_exec <- names(files)[!same & !same_code & !names(files) %in% adapted]
  if (length(bad_exec)) prob <- c(prob, paste("differ from the execution-check commit:", paste(bad_exec, collapse = ", ")))
  recorded <- engine_files[!vapply(engine_files, is.null, TRUE)]
  if (length(recorded) < length(engine_files))
    prob <- c(prob, sprintf("%d run folders recorded no engine-file hashes", length(engine_files) - length(recorded)))
  unknown <- character(0); bad_run <- character(0)
  for (ef in recorded) for (p in names(ef)) {
    if (!p %in% names(files)) unknown <- c(unknown, p) else if (!identical(as.character(ef[[p]]), now[[p]])) bad_run <- c(bad_run, p)
  }
  if (length(unknown)) prob <- c(prob, paste("hashed files not in the release record:", paste(unique(unknown), collapse = ", ")))
  if (length(bad_run)) prob <- c(prob, paste("run folders recorded other code for:", paste(unique(bad_run), collapse = ", ")))
  out$ok <- !length(prob) && length(recorded) > 0L
  keep <- !names(files) %in% adapted
  out$detail <- if (out$ok) sprintf("%s; %d run folders; %d engine files with the code of %s (%d with comments edited), adapted: %s",
                                    id, length(recorded), sum((same | same_code) & keep), substr(ec$commit, 1, 7),
                                    sum(same_code & keep), paste(adapted, collapse = ", "))
                else paste(prob, collapse = "; ")
  out
}

grade_run <- function(in_root, round3_root, cells, out_dir, frozen_sha, pass = GRADE_PASS, n_boot = GRADE_B_BOOT,
                      root = SIM_PROD_ROOT, frozen = NULL, assembly_roots = character(0), release_root = release_root_default(),
                      say = function(...) { cat(sprintf(...), "\n", sep = ""); flush(stdout()) }) {
  t0 <- proc.time()[["elapsed"]]
  dir.create(file.path(out_dir, "displays"), recursive = TRUE, showWarnings = FALSE)
  r3 <- read_round3(round3_root, assembly_roots)
  if (is.null(frozen)) frozen <- frozen_cells(ae_buildable = if (is.null(r3$kappa)) TRUE else isTRUE(r3$kappa$buildable))
  G <- list(); DSP <- list(); CMP <- list(); PROV <- list(); IMPL <- list(); FAILS <- list(); BLK <- list(); HOLD <- list()
  ENG <- list()                                       # engine-file SHA-256 recorded by each used run folder
  add_grade <- function(claim, role, cell, arm, population, statistic, res, premise = "", note = "", spec_ref = "", failed = NA_integer_) {
    G[[length(G) + 1L]] <<- cbind(data.frame(claim = claim, role = role, cell = cell, arm = arm, population = population,
                                             statistic = statistic, stringsAsFactors = FALSE), res,
                                  data.frame(failed_replicates = failed, premise = premise, note = note, spec_ref = spec_ref, stringsAsFactors = FALSE))
    utils::tail(res$disposition, 1)
  }
  add_display <- function(name, df) if (nrow(df)) DSP[[name]] <<- rbind_cols_fill(DSP[[name]], df)
  key <- function(cellkey, label) list(root = root, cell = cellkey, statistic = label, pass = pass)
  g_none <- function(claim, role, cellid, arm, pop, stat_name, disposition, estimate = NA_real_, region = NULL, premise = "", note = "", failed = NA_integer_)
    add_grade(claim, role, cellid, arm, pop, stat_name, ungraded_row(disposition, estimate, if (is.null(region)) c(NA, NA) else GRADE_MARGINS[[region]]),
              premise, note, if (is.null(region)) "" else GRADE_SPEC_REF[[region]], failed)

  # The gate of one cell: incomplete output, the required-check audit, the FAILED fits (every folder).
  make_gate <- function(cp, au, failed_by_arm, why_not) function(arms, extra_block = NULL, rung = NULL) {
    if (cp$status != "COMPLETE") return(list(state = "notgraded", disp = why_not))
    blocks <- unique(c(au$all, unlist(au$arms[intersect(arms, names(au$arms))], use.names = FALSE),
                       if (!is.null(rung)) au$rungs[[as.character(rung)]], extra_block))
    if (length(blocks)) return(list(state = "blocked", disp = paste0("BLOCKED (", paste(blocks, collapse = "; "), ")")))
    nf <- failed_by_arm[intersect(arms, names(failed_by_arm))]; nf <- nf[nf > 0]
    list(state = "ok", failed = sum(nf),
         fail_why = if (length(nf)) sprintf("fitting failures: %s FAILED (R4); survivors' statistic descriptive (author ruling R25)",
                                            paste(sprintf("%s %d", names(nf), as.integer(nf)), collapse = ", ")) else NULL)
  }

  # ONE graded quantity: gate -> premise -> fitting failures ->. kind: joint | draws | census.
  # disp_map (optional) rewrites a graded disposition (power controls, boundary/stress rows).
  grade_stat <- function(gate, claim, role, cellid, cellkey, arm_label, arms_used, d, stat_name, stat_fun, region,
                         kind = "joint", pop = NULL, premise = list(ok = TRUE, why = ""), note = "", nonfinite = 0L,
                         rung = NULL, disp_map = NULL, values = NULL, unresolved = 0L) {
    g <- gate(arms_used, extra_block = if (nonfinite > 0) sprintf("non-finite required value in %d replicate(s) marked OK", as.integer(nonfinite)), rung = rung)
    pid <- if (is.null(pop)) NA_integer_ else pop
    fl <- if (is.null(g$failed)) NA_integer_ else g$failed
    if (g$state != "ok") return(g_none(claim, role, cellid, arm_label, pid, stat_name, g$disp, region = region, premise = premise$why, note = note, failed = fl))
    if (!premise$ok) return(g_none(claim, role, cellid, arm_label, pid, stat_name,
                                   paste0("UNGRADABLE (", paste(c(paste("premise:", premise$why), g$fail_why), collapse = "; "), ")"),
                                   region = region, premise = premise$why, note = note, failed = fl))
    n_obs <- if (kind == "census") length(values) else nrow(d)
    if (!n_obs) return(g_none(claim, role, cellid, arm_label, pid, stat_name, "NOT GRADED (no retained replicates)", region = region, failed = fl))
    lab <- stat_label(claim, arm_label, stat_name, if (!is.null(pop)) paste0("pop", pop))
    b <- switch(kind, joint = boot_joint_stat(d, d$pop_id, stat_fun, key(cellkey, lab), n_boot),
                draws = boot_draws_stat(d, stat_fun, key(cellkey, lab), n_boot),
                census = boot_census_stat(values, stat_fun, key(cellkey, lab), n_boot))
    res <- grade_boot(b, GRADE_MARGINS[[region]])
    unres_why <- if (is.finite(unresolved) && unresolved > 0)
      sprintf("unresolved T2: %d replicate(s) of %s; survivors' statistic descriptive (round 3d rule, DECIDED 2026-10-16)", as.integer(unresolved), arm_label) else NULL
    if (!is.null(g$fail_why) || !is.null(unres_why)) {
      res$disposition <- paste0("UNGRADABLE (", paste(c(g$fail_why, unres_why), collapse = "; "), ")")
      res$basis <- paste0("descriptive (survivors only", if (!is.null(g$fail_why)) ", R25", if (!is.null(unres_why)) ", unresolved T2 excluded", ")")
    }
    else if (!is.null(disp_map)) { m <- disp_map(res$disposition); res$disposition <- m$disposition; if (!is.null(m$basis)) res$basis <- m$basis }
    add_grade(claim, role, cellid, arm_label, pid, stat_name, res, premise$why, note, GRADE_SPEC_REF[[region]], fl)
  }
  Z3 <- list(mean_Z = list(STAT$mean_z("Z_wor"), "z_location"), sd_Z_minus_1 = list(STAT$sd_z_minus1("Z_wor"), "z_scale"),
             KS_Z = list(STAT$ks_z("Z_wor"), "z_ks"))
  z_margins <- function(gate, claim, role, cellid, cellkey, arm, d, kind = "joint", pop = NULL, premise = list(ok = TRUE, why = ""),
                        note = "", nonfinite = 0L, rung = NULL)
    vapply(names(Z3), function(s) grade_stat(gate, claim, role, cellid, cellkey, arm, arm, d, s, Z3[[s]][[1]], Z3[[s]][[2]],
                                              kind = kind, pop = pop, premise = premise, note = note, nonfinite = nonfinite, rung = rung), "")

  rbL2 <- rb_premise(r3$rb, "L2"); rbQ <- rb_premise(r3$rb, "quarter")
  OKP <- list(ok = TRUE, why = "")
  arm_premise <- function(arm, what) {
    if (!startsWith(arm, "Rb")) return(list(ok = TRUE, why = "R-a: L2 and quarter-rate premises satisfied by construction (§12.11)"))
    if (what == "C2") rbL2 else if (rbL2$ok) rbQ else rbL2
  }
  ctrl_of <- function(pattern) { g <- r3$g6; if (!nrow(g)) return(NULL); i <- grep(pattern, g$control); if (length(i)) g[i[1], ] else NULL }

  cells <- cells[order(!grepl("^P-ms-H(800|400)$", cells$id)), , drop = FALSE]
  for (ci in seq_len(nrow(cells))) {
    row <- cells[ci, ]; id <- row$id
    t_c <- proc.time()[["elapsed"]]
    x <- tryCatch(read_prod_cell(in_root, row, frozen_sha, root, frozen), error = function(e) e)
    if (inherits(x, "error")) { CMP[[id]] <- data.frame(cell = id, status = "REFUSED", detail = paste("read error:", conditionMessage(x)), missing_tasks = "", runs = 0L, superseded = 0L); next }
    cp <- cell_completeness(x)
    au <- if (cp$status == "COMPLETE") required_check_audit(x) else list(all = character(0), arms = list(), rungs = list())
    failed_by_arm <- if (nrow(x$failed)) table(x$failed$arm) else table(character(0))
    failed_by_arm <- stats::setNames(as.integer(failed_by_arm), names(failed_by_arm))
    CMP[[id]] <- data.frame(cell = id, status = cp$status, detail = cp$detail, missing_tasks = paste(cp$missing_tasks, collapse = ","),
                            runs = length(x$folders$used), superseded = length(x$folders$superseded), stringsAsFactors = FALSE)
    if (nrow(x$provenance)) PROV[[id]] <- cbind(cell = id, x$provenance)
    ENG[[id]] <- lapply(x$manifests, function(m) m$engine_files_sha256)
    blk <- c(au$all, unlist(lapply(names(au$arms), function(a) paste0(a, ": ", au$arms[[a]]))), unlist(lapply(names(au$rungs), function(h) paste0("H = ", h, ": ", au$rungs[[h]]))))
    if (length(blk)) BLK[[id]] <- data.frame(cell = id, block = blk, stringsAsFactors = FALSE)
    if (nrow(x$checks) && "required" %in% names(x$checks)) {
      req <- as.logical(x$checks$required) %in% TRUE; ok <- as.logical(x$checks$pass) %in% TRUE
      IMPL[[id]] <- data.frame(cell = id, required_checks = sum(req), failed = sum(req & !ok),
                               declared = sum(x$required_declared, na.rm = TRUE),
                               first_failure = if (any(req & !ok)) x$checks$check[which(req & !ok)[1]] else "", stringsAsFactors = FALSE)
    }
    if (nrow(x$results) && "arm" %in% names(x$results)) {
      arms_all <- sort(unique(x$results$arm))
      FAILS[[id]] <- data.frame(cell = id, arm = arms_all, replicates = as.integer(table(factor(x$results$arm, levels = arms_all))),
                                failed = as.integer(vapply(arms_all, function(a) if (a %in% names(failed_by_arm)) failed_by_arm[[a]] else 0L, 0L)),
                                failed_in_superseded_runs_only = if (x$failed_superseded_only > 0) x$failed_superseded_only else 0L, stringsAsFactors = FALSE)
    }
    complete <- cp$status == "COMPLETE"
    why_not <- if (complete) "" else sprintf("NOT GRADED (%s output: %s)", tolower(cp$status), cp$detail)
    gate <- make_gate(cp, au, failed_by_arm, why_not)
    cellkey <- if (length(x$manifests) && !is.null(x$manifests[[1]]$config$cell)) x$manifests[[1]]$config$cell else NULL
    r <- x$results
    fa <- function(a) if (a %in% names(failed_by_arm)) failed_by_arm[[a]] else 0L
    say("[grade] %s: %s (%s)%s", id, cp$status, cp$detail, if (length(blk)) sprintf("; %d block(s)", length(blk)) else "")
    m_ms <- if (grepl("^P-ms-H", id)) 2L * as.integer(sub("^P-ms-H", "", id)) else NA_integer_
    m_fs <- if (grepl("^P-fs-m", id)) as.integer(sub("^P-fs-m", "", id)) else NA_integer_

    err <- tryCatch({
    # ---------------- E1 oracle cells (design-conditional, 12 x 4000) ----------------
    if (grepl("^P-(ms|fs)-or-", id)) {
      regime <- if (grepl("^P-ms", id)) "many-strata" else "fixed-strata"
      terminal <- id %in% c("P-ms-or-H800", "P-fs-or-m1600")
      if (nrow(r)) for (p in sort(unique(r$pop_id))) {
        a <- arm_rows(r[r$pop_id == p, , drop = FALSE], ARM_ORACLE, c("Z_wor", "R_var", "psi_hat"))
        d <- a$rows
        if (nrow(d) >= 2L) add_display("e1_oracle_by_population", data.frame(cell = id, regime = regime, population = p, rows = nrow(d),
          mean_Z = mean(d$Z_wor), sd_Z = stats::sd(d$Z_wor), ks_Z = f_ks(d$Z_wor, rep(1, nrow(d))), mean_R_var = mean(d$R_var), sd_R_var = stats::sd(d$R_var),
          linearization_monitor = stats::var(d$psi_hat) / d$sigma_m_wor[1]^2, cover95_studentized = mean(as.logical(d$cover95))))
        if (!terminal) next
        z_margins(gate, "C1 Z_oracle", "implementation check", id, cellkey, ARM_ORACLE, d, kind = "draws", pop = p, nonfinite = a$nonfinite_ok,
                  note = sprintf("%s oracle at the terminal rung m = 1600; population-indexed (R20), never composed", regime))
        if (regime == "many-strata") {
          grade_stat(gate, "C2 oracle limb", "implementation check", id, cellkey, ARM_ORACLE, ARM_ORACLE, d, "mean_R_var_minus_1", STAT$mean_minus1("R_var"),
                     "rvar_oracle_mean", kind = "draws", pop = p, nonfinite = a$nonfinite_ok,
                     note = "analytic-vs-MC normalisation (E[R_var] = 1 by definition); a breach is a coding fault, never C2 evidence")
          grade_stat(gate, "C2 oracle limb", "implementation check", id, cellkey, ARM_ORACLE, ARM_ORACLE, d, "sd_R_var", STAT$sd("R_var"),
                     "rvar_oracle_sd", kind = "draws", pop = p, nonfinite = a$nonfinite_ok, note = "oracle-tier concentration within a fixed population (R20)")
        }
      }
    }

    # ---------------- E1 many-strata fitted cells (joint, 400 x 5) ----------------
    if (grepl("^P-ms-H", id)) {
      orc <- if (nrow(r)) arm_rows(r, ARM_ORACLE, c("R_var"))$rows else r
      for (a in unique(r$arm)) {
        ar <- arm_rows(r, a, c("Z_wor", "R_var")); d <- ar$rows
        s <- arm_summary(d)
        if (a != ARM_ORACLE && nrow(d) && nrow(orc)) {
          j <- merge(d[, c("pop_id", "rep_id", "R_var")], orc[, c("pop_id", "rep_id", "R_var")], by = c("pop_id", "rep_id"), suffixes = c("", "_or"))
          s$transfer_mean <- f_mean(abs(j$R_var - j$R_var_or), pop_freq(j$pop_id))
        }
        add_display("e1_many_strata_by_arm", cbind(data.frame(cell = id, m = m_ms, arm = a, failed = fa(a), basis = basis_of(fa(a))), s))
      }
      grade_c1c2 <- function(a) {
        ar <- arm_rows(r, a, c("Z_wor", "R_var")); d <- ar$rows
        or_a <- arm_rows(r, ARM_ORACLE, "R_var")
        z_margins(gate, "C1", "scientific", id, cellkey, a, d, premise = arm_premise(a, "C1"), nonfinite = ar$nonfinite_ok,
                  note = sprintf("terminal rung m = %d (§3/E1)", m_ms))
        pr2 <- arm_premise(a, "C2")
        grade_stat(gate, "C2", "scientific", id, cellkey, a, a, d, "mean_R_var_minus_1", STAT$mean_minus1("R_var"), "rvar_fit_mean", premise = pr2, nonfinite = ar$nonfinite_ok)
        grade_stat(gate, "C2", "scientific", id, cellkey, a, a, d, "sd_R_var", STAT$sd("R_var"), "rvar_fit_sd", premise = pr2, nonfinite = ar$nonfinite_ok)
        j <- merge(d[, c("pop_id", "rep_id", "R_var")], or_a$rows[, c("pop_id", "rep_id", "R_var")], by = c("pop_id", "rep_id"), suffixes = c("", "_or"))
        j$transfer <- j$R_var - j$R_var_or         # (se2_fit - se2_or)/sigma2_m,b, same replicate
        grade_stat(gate, "C2", "scientific", id, cellkey, a, c(a, ARM_ORACLE), j, "transfer", STAT$mean_abs("transfer"), "transfer", premise = pr2,
                   nonfinite = ar$nonfinite_ok + or_a$nonfinite_ok, note = "per-replicate |R_var,fit - R_var,oracle| averaged (§12.4 :1291)")
      }
      if (id == "P-ms-H800") grade_c1c2(ARM_RA_AIPW)
      if (id == "P-ms-H400") {
        grade_c1c2(ARM_RB_AIPW)
        HOLD$tmle <- r[r$arm %in% c(ARM_RA_TMLE, ARM_RB_TMLE), , drop = FALSE]; HOLD$tmle_key <- cellkey; HOLD$gate400 <- gate
      }
    }

    # ---------------- E1 fixed-strata fitted cells (joint; C4 and displays) ----------------
    if (grepl("^P-fs-m", id)) {
      for (a in unique(r$arm)) {
        ar <- arm_rows(r, a, c("psi_hat")); d <- ar$rows
        s <- arm_summary(d)
        s$t2_unresolved <- if ("T2_status" %in% names(r)) sum(r$arm == a & r$T2_status %in% "unresolved") else 0L
        if (grade_remainder_excluded(a)) {                  # no T1/T2 column for this arm, a note instead
          s$t2_unresolved <- NA_integer_; s$remainder_note <- GRADE_R28_NOTE
        } else if (nrow(d) && all(c("T1_over_sigma_wor", "T2_over_sigma_wor") %in% names(d)) && any(is.finite(d$T1_over_sigma_wor))) {
          dd <- d[is.finite(d$T1_over_sigma_wor) & is.finite(d$T2_over_sigma_wor), , drop = FALSE]
          f <- pop_freq(dd$pop_id)
          s$q50_abs_T1_std <- f_q7(abs(dd$T1_over_sigma_wor), f, 0.5); s$q90_abs_T1_std <- f_q7(abs(dd$T1_over_sigma_wor), f, 0.9)
          s$q50_abs_T2_std <- f_q7(abs(dd$T2_over_sigma_wor), f, 0.5); s$q90_abs_T2_std <- f_q7(abs(dd$T2_over_sigma_wor), f, 0.9)
          s$q50_sqrt_m_abs_T1 <- f_q7(sqrt(m_fs) * abs(dd$T1), f, 0.5); s$q50_sqrt_m_abs_T2 <- f_q7(sqrt(m_fs) * abs(dd$T2), f, 0.5)
          if ("r_hat" %in% names(dd) && any(is.finite(dd$r_hat))) {
            k <- is.finite(dd$r_hat) & dd$r_hat > 0
            if (any(k)) s$q90_sqrt_m_T1_over_rhat <- f_q7(sqrt(m_fs) * abs(dd$T1[k]) / dd$r_hat[k], pop_freq(dd$pop_id[k]), 0.9)
            k2 <- is.finite(dd$r_Q) & is.finite(dd$r_g) & dd$r_Q * dd$r_g > 0
            if (any(k2)) s$q90_T2_over_rQ_rg <- f_q7(abs(dd$T2[k2]) / (dd$r_Q[k2] * dd$r_g[k2]), pop_freq(dd$pop_id[k2]), 0.9)
            if (any(k)) s$mean_r_hat <- f_mean(dd$r_hat[k], pop_freq(dd$pop_id[k]))
          }
        }
        add_display("e1_fixed_strata_by_arm", cbind(data.frame(cell = id, m = m_fs, arm = a, failed = fa(a), basis = basis_of(fa(a))), s))
      }
      if (id == "P-fs-m1600") {
        ar <- arm_rows(r, ARM_RA_AIPW, "T1_over_sigma_wor")
        grade_stat(gate, "C4a", "scientific", id, cellkey, ARM_RA_AIPW, ARM_RA_AIPW, ar$rows, "Q90_abs_T1_std", STAT$q90_abs("T1_over_sigma_wor"), "negligible",
                   premise = list(ok = TRUE, why = "R-a: r_m -> 0 by construction (§12.11)"), note = "terminal rung m = 1600", nonfinite = ar$nonfinite_ok,
                   unresolved = ar$unresolved)
      }
      if (id == "P-fs-m3200") for (al in RD_ALPHAS) {
        a <- rd_arm(al)
        role_b <- if (al > 0.25) "scientific" else if (abs(al - 0.25) < 1e-9) "boundary" else "stress"
        if (grade_rd_stopped(al)) {                     # STOPPED after G5; not run in production
          why <- sprintf("UNGRADABLE (premise: G5 alpha %.2f STOPPED after G5 -- failed engineering, recorded in round 3c; ruling R27; not run in production)", al)
          g_none("C4a", "scientific", id, a, NA_integer_, "Q90_abs_T1_std", why, region = "negligible", premise = g5_premise(r3$g5, al)$why, note = "terminal rung m = 3200")
          g_none("C4b", role_b, id, a, NA_integer_, "Q90_abs_T2_std", why, region = "negligible", premise = g5_premise(r3$g5, al)$why, note = "terminal rung m = 3200")
          next
        }
        ar <- arm_rows(r, a, c("T1_over_sigma_wor", "T2_over_sigma_wor")); pr <- g5_premise(r3$g5, al)
        grade_stat(gate, "C4a", "scientific", id, cellkey, a, a, ar$rows, "Q90_abs_T1_std", STAT$q90_abs("T1_over_sigma_wor"), "negligible",
                   premise = pr, note = "terminal rung m = 3200", nonfinite = ar$nonfinite_ok, unresolved = ar$unresolved)
        dm <- if (role_b == "scientific") NULL else function(dsp) list(
          disposition = if (role_b == "boundary") "BOUNDARY (alpha = 1/4: O not o; reported descriptively, never PASS/RED-FLAG)" else "STRESS (alpha < 1/4: no prediction)",
          basis = "descriptive (premise not claimed)")
        grade_stat(gate, "C4b", role_b, id, cellkey, a, a, ar$rows, "Q90_abs_T2_std", STAT$q90_abs("T2_over_sigma_wor"), "negligible",
                   premise = if (role_b == "scientific") pr else OKP, note = "terminal rung m = 3200", nonfinite = ar$nonfinite_ok, disp_map = dm,
                   unresolved = ar$unresolved)
      }
    }

    # ---------------- E2 FPC oracle cells (design-conditional): C3b, C3a ----------------
    if (grepl("^P-fpc-or-J", id)) {
      J <- as.integer(sub("^P-fpc-or-J", "", id))
      if (nrow(r)) for (p in sort(unique(r$pop_id))) {
        aa <- arm_rows(r[r$pop_id == p, , drop = FALSE], ARM_ORACLE, c("Z_wor", "psi_hat", "se_hat"))
        d_all <- aa$rows
        d <- if (J == 4L) d_all[d_all$rep_id <= 2000L, , drop = FALSE] else d_all
        if (nrow(d) < 2L) next
        sm <- d$sigma_m[1]; sw <- d$sigma_m_wor[1]
        mH <- as.integer(cellkey$m %||% 100L) - as.integer(cellkey$H %||% 50L)   # design df m - H (FPC cells: 100 - 50)
        disp <- data.frame(cell = id, J = J, f = 2 / J, population = p, draws = nrow(d), sigma_ratio = sm / sw)
        for (lv in GRADE_LEVELS) {
          z <- stats::qnorm(1 - (1 - lv) / 2); tq <- stats::qt(1 - (1 - lv) / 2, mH)
          disp[[sprintf("cov_wor_%02d", 100 * lv)]] <- mean(abs(d$psi_hat - d$Psi_N) <= z * sw)
          disp[[sprintf("cov_det_sigma_m_%02d", 100 * lv)]] <- mean(abs(d$psi_hat - d$Psi_N) <= z * sm)
          disp[[sprintf("cov_se_%02d", 100 * lv)]] <- mean(abs(d$psi_hat - d$Psi_N) <= z * d$se_hat)
          disp[[sprintf("cov_se_t_df%d_%02d", mH, 100 * lv)]] <- mean(abs(d$psi_hat - d$Psi_N) <= tq * d$se_hat)
          disp[[sprintf("overlay_%02d", 100 * lv)]] <- 2 * stats::pnorm(z * sm / sw) - 1
        }
        disp$c3a_statistic <- STAT$c3a()(d_all, rep(1, nrow(d_all)))
        disp$var_ratio_wor <- stats::var(d_all$psi_hat) / sw^2; disp$mean_se2_over_sigma2_m <- mean(d_all$se_hat^2) / sm^2
        add_display("e2_fpc_oracle_by_population", disp)
        zd <- z_margins(gate, "C3b premise: Z_oracle", "premise check", id, cellkey, ARM_ORACLE, d, kind = "draws", pop = p, nonfinite = aa$nonfinite_ok,
                        note = sprintf("decided per population on its own %d redraws (§6 :775-778)", nrow(d)))
        adequate <- all(zd == "PASS")
        for (lv in GRADE_LEVELS)
          grade_stat(gate, "C3b", "scientific", id, cellkey, ARM_ORACLE, ARM_ORACLE, d, sprintf("coverage_minus_%.2f", lv), STAT$cov_det_sigma_m(lv), "coverage",
                     kind = "draws", pop = p, nonfinite = aa$nonfinite_ok,
                     premise = if (adequate) list(ok = TRUE, why = "Z_oracle distributionally adequate (all three C1 margins PASS)")
                               else list(ok = FALSE, why = paste("this population's Z_oracle is not distributionally adequate:", paste(zd, collapse = " / "))),
                     note = "interval (2): oracle +- z_{1-gamma/2} sigma_m (deterministic sigma_m)")
        if (J == 4L) grade_stat(gate, "C3a", "scientific", id, cellkey, ARM_ORACLE, ARM_ORACLE, d_all, "c3a_signed", STAT$c3a(), "c3a", kind = "draws", pop = p,
                                nonfinite = aa$nonfinite_ok,
                                note = sprintf("%d SRS-WOR redraws (R9); G1 binding cleared in %s of %s populations (preparation)", nrow(d_all),
                                               r3$prep$g1_binding$clear %||% "?", r3$prep$g1_binding$n %||% "?"))
      }
    }

    # ---------------- E2 fitted cells (joint, 160 x 5): C1; fitted C3b analogues displayed ----------------
    if (id %in% c("P-fpc-fit-J4", "P-strict-m200")) {
      for (a in intersect(c(ARM_ORACLE, ARM_RA_AIPW, ARM_RB_AIPW), unique(r$arm))) {
        ar <- arm_rows(r, a, c("Z_wor", "R_var", "psi_hat", "se_hat")); d <- ar$rows
        s <- arm_summary(d)
        if (nrow(d)) { f <- pop_freq(d$pop_id); for (lv in GRADE_LEVELS) s[[sprintf("cov_se_%02d", 100 * lv)]] <- f_mean(as.numeric(abs(d$psi_hat - d$Psi_N) <= stats::qnorm(1 - (1 - lv) / 2) * d$se_hat), f) }
        add_display("e2_fitted_by_arm", cbind(data.frame(cell = id, arm = a, failed = fa(a), basis = basis_of(fa(a))), s))
        if (a == ARM_ORACLE) next
        z_margins(gate, "C1", "scientific", id, cellkey, a, d, premise = arm_premise(a, "C1"), nonfinite = ar$nonfinite_ok,
                  note = "E2 fitted cell, f = 0.5 (R7: carries a C1 grade)")
      }
    }

    # ---------------- census cells ----------------
    if (id == "P-vanish-census") {
      tr <- x$trajectory
      if (nrow(tr)) {
        tr$terminal <- as.logical(tr$terminal)
        sp <- split(tr, paste(tr$H, tr$J))
        add_display("e2_vanishing_trajectory", do.call(rbind, lapply(sp, function(t) data.frame(H = t$H[1], J = t$J[1], f = 2 / t$J[1], populations = nrow(t),
          median_bridge = stats::median(t$bridge), q90_bridge = unname(stats::quantile(t$bridge, 0.9, type = 7, na.rm = TRUE)), max_bridge = max(t$bridge),
          mean_Delta_m = mean(t$Delta_m), var_Psi_N_minus_Psi = stats::var(t$Psi_N_minus_Psi),
          ratio_mean_Delta_over_var_Smod = mean(t$Delta_m) / stats::var(t$Psi_N_minus_Psi), terminal = any(t$terminal)))))
      }
      tv <- if (nrow(tr)) tr[tr$terminal %in% TRUE, , drop = FALSE] else tr
      vals <- if (nrow(tv)) as.numeric(tv$bridge) else numeric(0)
      ck <- list(regime = "many", H = 200L, m_h = 2L, J = 800L, m = 400L, f = "1/400", e2c = "vanishing")
      grade_stat(gate, "sec. Psi", "scientific", id, ck, "census", character(0), NULL, "Q90_bridge", function(v) unname(stats::quantile(v, 0.9, type = 7)), "negligible",
                 kind = "census", values = vals[is.finite(vals)], nonfinite = sum(!is.finite(vals)),
                 premise = list(ok = TRUE, why = "vanishing fraction (sufficient condition m/M_pop -> 0)"),
                 note = sprintf("terminal rung (200, 800), %d populations, population-only percentile bootstrap (R6)", sum(is.finite(vals))))
    }
    if (id == "P-strict-census" && nrow(x$trajectory)) {
      tr <- x$trajectory
      add_display("e2_strictness_witness", do.call(rbind, lapply(split(tr, tr$H), function(t) data.frame(H = t$H[1], populations = nrow(t),
        min_Delta_over_sigma2_m = min(t$Delta_over_sigma2_m), median_Delta_over_sigma2_m = stats::median(t$Delta_over_sigma2_m),
        max_Delta_over_sigma2_m = max(t$Delta_over_sigma2_m), median_sigma_ratio = stats::median(t$sigma_ratio)))))
    }
    if (id == "P-unequal-census" && nrow(x$pair)) {
      pr <- x$pair
      add_display("e2_unequal_pair", do.call(rbind, lapply(split(pr, pr$alloc), function(t) data.frame(alloc = t$alloc[1], populations = nrow(t),
        median_Delta_over_sigma2_m = stats::median(t$Delta_over_sigma2_m), median_Delta_m = stats::median(t$Delta_m)))))
    }

    # ---------------- E3 ----------------
    if (id == "P-e3-rc-m400" || grepl("^P-e3-vsens", id)) {
      for (a in unique(r$arm)) {
        ar <- arm_rows(r, a, c("Z_wor", "R_var")); s <- arm_summary(ar$rows)
        add_display(if (id == "P-e3-rc-m400") "e3_rc_foil_by_arm" else "e3_v_sensitivity_by_arm",
                    cbind(data.frame(cell = id, arm = a, failed = fa(a), basis = basis_of(fa(a)),
                                     c2_status = if (startsWith(a, "Rc")) "UNGRADABLE (R-c: premise-failing tier, reported only, §12.4 :1368)" else "diagnostic (E3 carries no C2 grade)"), s))
      }
    }
    if (id == "P-e3-ctrl-m400") {
      HOLD$ae_sha <- tryCatch(x$manifests[[1]]$config$round2$ctrl$ae$content_sha256, error = function(e) NULL) %||% NA_character_
      host <- arm_rows(r, ARM_RA_AIPW, c("Z_shift", "R_var_nomh", "R_var_unit", "R_var_nostrat"))
      ae <- arm_rows(r, ARM_AE_AIPW, "T2_over_sigma_wor")
      ctrls <- list(
        list(lab = "C1 shift: mean Z keyed to Psi_N + 0.30 sigma_m,wor (R8)", pat = "^C1 shift", src = host, arm = ARM_RA_AIPW, st = STAT$mean_z("Z_shift"), reg = "z_location", name = "mean_Z_shift"),
        list(lab = "C2 wrong clustering: no m_h/(m_h-1)", pat = "no m_h", src = host, arm = ARM_RA_AIPW, st = STAT$mean("R_var_nomh"), reg = "wrong_clustering", name = "mean_R_var_nomh"),
        list(lab = "C2 wrong clustering: unit-level clusters", pat = "unit-level", src = host, arm = ARM_RA_AIPW, st = STAT$mean("R_var_unit"), reg = "wrong_clustering", name = "mean_R_var_unit"),
        list(lab = "C2 wrong clustering: no stratification", pat = "no stratification", src = host, arm = ARM_RA_AIPW, st = STAT$mean("R_var_nostrat"), reg = "wrong_clustering", name = "mean_R_var_nostrat"),
        list(lab = "C4 aligned error: Q_0.9(|T2|/sigma_m,wor), AE-AIPW", pat = "^C4 aligned", src = ae, arm = ARM_AE_AIPW, st = STAT$q90_abs("T2_over_sigma_wor"), reg = "negligible", name = "Q90_abs_T2_std"))
      for (cc in ctrls) {
        g6 <- ctrl_of(cc$pat)
        role6 <- if (is.null(g6)) "G6 limb (a) row missing" else if (isTRUE(as.logical(g6$graded))) "graded power control" else g6$status
        dm <- if (identical(role6, "graded power control")) function(dsp) list(disposition = control_disposition(dsp)) else
          function(dsp) list(disposition = paste0("DIAGNOSTIC (G6 limb (a): ", role6, "; never graded as a power control)"), basis = "descriptive (diagnostic)")
        grade_stat(gate, "Power control", "power control", id, cellkey, cc$lab, cc$arm, cc$src$rows, cc$name, cc$st, cc$reg,
                   premise = list(ok = TRUE, why = role6), note = "predicted to fail its row (§8 POWER CONTROLS)", nonfinite = cc$src$nonfinite_ok, disp_map = dm)
      }
      if (nrow(host$rows) && "Z_gf" %in% names(host$rows))
        add_display("e3_controls_reported", data.frame(control = "C1 g-formula-keyed Z (reported, no verdict, R8)", mean = f_mean(host$rows$Z_gf, pop_freq(host$rows$pop_id)),
                                                       rows = nrow(host$rows), basis = basis_of(fa(ARM_RA_AIPW))))
    }

    # ---------------- E4 (C7, implementation check) ----------------
    if (id == "P-e4-array" && nrow(r)) {
      ec <- x$e4_census[x$e4_census$array %in% "regular", , drop = FALSE]
      cfg <- if (length(x$manifests)) x$manifests[[1]]$config else list()
      for (H in sort(unique(r$H))) {
        e <- ec[ec$H == H, , drop = FALSE]; d <- r[r$H == H, , drop = FALSE]
        if (!nrow(e)) next
        d$Z_wor <- (d$psi_hat - e$Psi_N[1]) / sqrt(e$sigma2_m_wor[1])
        nonfin <- sum(!is.finite(d$Z_wor)); d <- d[is.finite(d$Z_wor), , drop = FALSE]
        if (nrow(d) >= 2L) add_display("e4_by_rung", data.frame(H = H, draws = nrow(d), mean_Z = mean(d$Z_wor), sd_Z = stats::sd(d$Z_wor), ks_Z = f_ks(d$Z_wor, rep(1, nrow(d))),
                                             var_ratio_wor = stats::var(d$psi_hat) / e$sigma2_m_wor[1], mean_se2_over_sigma2_m = mean(d$se2_hat) / e$sigma2_m[1],
                                             lindeberg = e$lindeberg[1]))
        mh <- as.integer(cfg$m_h %||% 2L); J <- as.integer(cfg$J %||% 20L)
        ck <- list(regime = "many", H = as.integer(H), m_h = mh, J = J, m = mh * as.integer(H), f = frac_str(mh, J), e2c = "e4")
        z_margins(gate, "C7", "implementation check", id, ck, ARM_ORACLE, d, kind = "draws", pop = H, nonfinite = nonfin, rung = H,
                  note = sprintf("E4 rung H = %d (population column = H)", H))
      }
    }
    NULL }, error = function(e) conditionMessage(e))
    if (!is.null(err)) {
      CMP[[id]]$status <- "GRADING ERROR"; CMP[[id]]$detail <- paste(CMP[[id]]$detail, "| grading error:", err)
      say("[grade] %s: GRADING ERROR -- %s", id, err)
    }
    say("[grade] %s graded in %.1f s", id, proc.time()[["elapsed"]] - t_c)
    rm(x, r); invisible(gc(FALSE))
  }

  # ---------------- C6 (needs the AIPW composites of P-ms-H800 and P-ms-H400) ----------------
  grades <- if (length(G)) do.call(rbind, G) else data.frame()
  cell_comp <- function(claim, cell, arm) {
    g <- grades[grades$claim == claim & grades$cell == cell & grades$arm == arm, , drop = FALSE]
    if (!nrow(g)) "UNGRADABLE (not graded)" else compose(g$disposition)
  }
  aipw_cell <- c(Ra = "P-ms-H800", Rb = "P-ms-H400")
  for (tm in c(ARM_RA_TMLE, ARM_RB_TMLE)) {
    fam <- substr(tm, 1, 2); aipw <- if (fam == "Ra") ARM_RA_AIPW else ARM_RB_AIPW
    stats_c6 <- c("Q90_abs_paired_psi_std", "Q90_abs_se2_diff_std")
    if (is.null(HOLD$gate400)) { for (s in stats_c6) g_none("C6", "scientific", "P-ms-H400", tm, NA, s, "NOT GRADED (no m = 800 output)", region = "negligible"); next }
    c1 <- cell_comp("C1", aipw_cell[[fam]], aipw); c2 <- cell_comp("C2", aipw_cell[[fam]], aipw)
    ar <- if (nrow(HOLD$tmle)) arm_rows(HOLD$tmle, tm, c("paired_tmle_minus_aipw", "se2_diff_over_sigma2_m", "sigma_m_wor")) else list(rows = data.frame(), failed = 0L, nonfinite_ok = 0L)
    d <- ar$rows
    pre <- arm_premise(aipw, "C1")
    gram <- if (nrow(d) && "gram_premise_met" %in% names(d)) all(as.logical(d$gram_premise_met) %in% TRUE) else FALSE
    why <- c(if (c1 != "PASS") sprintf("matched AIPW C1 (%s) is %s", aipw_cell[[fam]], c1),
             if (c2 != "PASS") sprintf("matched AIPW C2 (%s) is %s", aipw_cell[[fam]], c2),
             if (!pre$ok) pre$why,
             if (!gram) sprintf("Gram premise lambda_min >= c_H = %.2f not met in every retained replicate (min %s)", GRADE_TMLE_C_H,
                                if (nrow(d)) format(suppressWarnings(min(as.numeric(d$lambda_min), na.rm = TRUE)), digits = 4) else "NA"))
    prem <- if (length(why)) list(ok = FALSE, why = paste(why, collapse = "; ")) else list(ok = TRUE, why = sprintf("matched AIPW C1 and C2 PASS; Gram met in all %d replicates", nrow(d)))
    if (nrow(d)) d$paired_std <- d$paired_tmle_minus_aipw / d$sigma_m_wor
    grade_stat(HOLD$gate400, "C6", "scientific", "P-ms-H400", HOLD$tmle_key, tm, c(tm, aipw), d, "Q90_abs_paired_psi_std", STAT$q90_abs("paired_std"), "negligible",
               premise = prem, note = "m = 800 (TMLE terminal rung)", nonfinite = ar$nonfinite_ok)
    grade_stat(HOLD$gate400, "C6", "scientific", "P-ms-H400", HOLD$tmle_key, tm, c(tm, aipw), d, "Q90_abs_se2_diff_std", STAT$q90_abs("se2_diff_over_sigma2_m"), "negligible",
               premise = prem, note = "R19: |se2_TMLE - se2_AIPW| / sigma2_m,b", nonfinite = ar$nonfinite_ok)
  }
  grades <- if (length(G)) do.call(rbind, G) else data.frame()

  # ---------------- C5 (from G5 and C4a at a premise-violating rung, alpha < 1/4) ----------------
  # The premise-violation cells are alpha < 1/4; alpha = 1/4 is BOUNDARY, never a C5 witness.
  # The only such rung, alpha = 0.20, is STOPPED after G5.
  low <- RD_ALPHAS[RD_ALPHAS < 0.25 - 1e-9]
  c5_detail <- vapply(low, function(al) if (grade_rd_stopped(al)) sprintf("alpha %.2f: STOPPED after G5 (ruling R27)", al) else
                        sprintf("alpha %.2f: G5 %s; C4a %s", al, if (g5_premise(r3$g5, al)$ok) "PASS" else "not PASS", cell_comp("C4a", "P-fs-m3200", rd_arm(al))), "")
  c5_ok <- any(vapply(low, function(al) !grade_rd_stopped(al) && g5_premise(r3$g5, al)$ok && cell_comp("C4a", "P-fs-m3200", rd_arm(al)) == "PASS", TRUE))
  c5_why <- if (all(vapply(low, grade_rd_stopped, TRUE)))
    sprintf("UNGRADABLE (the only premise-violating rung, alpha < 1/4 (spec :128), is %s, STOPPED after G5 (ruling R27); alpha = 0.25 is BOUNDARY (§1.7a))",
            paste(sprintf("%.2f", low), collapse = ", ")) else "UNGRADABLE (the engineered rates or the numerical floor do not establish the separation cleanly)"
  add_grade("C5", "scientific", "P-fs-m3200", "R-d (alpha < 1/4)", NA, "premise separation",
            ungraded_row(if (c5_ok) "PASS" else c5_why, basis = "derived"),
            paste(c5_detail, collapse = "; "), "never RED-FLAG; never from T2 (§8 C5)", "§8 :949, :128")
  grades <- if (length(G)) do.call(rbind, G) else data.frame()

  # ---------------- composition: one disposition per joint (fitted) claim cell; design-conditional rows stay indexed ----------------
  joint <- grades[grades$role == "scientific" & is.na(grades$population) & grades$claim %in% c("C1", "C2", "C4a", "C6"), , drop = FALSE]
  composed <- if (nrow(joint)) do.call(rbind, lapply(split(joint, paste(joint$claim, joint$cell, joint$arm, sep = "|")), function(g)
    data.frame(claim = g$claim[1], cell = g$cell[1], arm = g$arm[1], metrics = nrow(g), disposition = compose(g$disposition),
               metric_dispositions = paste(sprintf("%s: %s", g$statistic, g$disposition), collapse = "; "), stringsAsFactors = FALSE))) else data.frame()
  single <- grades[(grades$claim %in% c("C4b", "C5", "sec. Psi")) | grades$role == "power control", , drop = FALSE]
  if (nrow(single)) composed <- rbind(composed, data.frame(claim = single$claim, cell = single$cell, arm = single$arm, metrics = 1L,
                                                           disposition = single$disposition, metric_dispositions = single$statistic, stringsAsFactors = FALSE))
  indexed <- grades[!is.na(grades$population), , drop = FALSE]
  idx_summary <- if (nrow(indexed)) do.call(rbind, lapply(split(indexed, paste(indexed$claim, indexed$cell, indexed$statistic, sep = "|")), function(g)
    data.frame(claim = g$claim[1], role = g$role[1], cell = g$cell[1], statistic = g$statistic[1], populations = nrow(g),
               PASS = sum(g$disposition == "PASS"), RED_FLAG = sum(g$disposition == "RED-FLAG"),
               UNGRADABLE = sum(startsWith(g$disposition, "UNGRADABLE")), BLOCKED = sum(startsWith(g$disposition, "BLOCKED")),
               NOT_GRADED = sum(startsWith(g$disposition, "NOT GRADED")), stringsAsFactors = FALSE))) else data.frame()

  # ---------------- outputs ----------------
  w <- function(df, f) if (!is.null(df) && length(df) && nrow(as.data.frame(df))) utils::write.csv(df, file.path(out_dir, f), row.names = FALSE, na = "NA")
  bind <- function(l) if (length(l)) do.call(rbind, l) else data.frame()
  completeness <- bind(CMP); provenance <- bind(PROV); impl <- bind(IMPL); failed <- bind(FAILS); blocks <- bind(BLK)
  w(grades, "grades.csv"); w(composed, "claims.csv"); w(idx_summary, "population_indexed_summary.csv")
  w(completeness, "completeness.csv"); w(provenance, "provenance.csv"); w(impl, "required_checks.csv"); w(failed, "failed_arms.csv")
  w(blocks, "blocked.csv"); w(r3$files, "round3_inputs.csv")
  for (nm in names(DSP)) utils::write.csv(DSP[[nm]], file.path(out_dir, "displays", paste0(nm, ".csv")), row.names = FALSE, na = "NA")
  commits <- unique(stats::na.omit(provenance$package_commit))
  r3_commit <- r3$exec$package_commit %||% NA_character_
  rel <- release_identity_check(commits, unlist(ENG, recursive = FALSE), release_root, r3_commit)
  checks_prov <- data.frame(
    check = c("every production cell agrees with the frozen design (cells.tsv vs frozen_cells())",
              "every production run folder: mode prod, the production root, the frozen v9.8 spec",
              "one package commit across every production run folder",
              if (rel$applies) rel$check else "that commit is the one the round-3 assembly PASSed (production guard)",
              "round-3 execution-gate verdict",
              "kappa_AE: the production panel's record is the frozen round-3 record (content SHA-256)"),
    ok = c(!any(grepl("^design mismatch", completeness$detail)),
           nrow(provenance) > 0 && all(provenance$ok), length(commits) == 1L, if (rel$applies) rel$ok else length(commits) == 1L && identical(commits, r3_commit),
           identical(r3$exec$verdict, "PASS"), identical(HOLD$ae_sha %||% NA_character_, r3$kappa$content_sha256 %||% "missing")),
    detail = c(paste(completeness$cell[grepl("^design mismatch", completeness$detail)], collapse = ", "),
               sprintf("%d folders, %d not ok", nrow(provenance), if (nrow(provenance)) sum(!provenance$ok) else 0L), paste(commits, collapse = ", "),
               if (rel$applies) rel$detail else sprintf("production %s; round 3 %s", paste(commits, collapse = ","), r3_commit),
               r3$exec$verdict %||% "missing",
               sprintf("production %s; round 3 %s", HOLD$ae_sha %||% "missing", r3$kappa$content_sha256 %||% "missing")),
    stringsAsFactors = FALSE)
  w(checks_prov, "provenance_checks.csv")
  list(provenance_checks = checks_prov, grades = grades, claims = composed, indexed = idx_summary, completeness = completeness, provenance = provenance,
       required_checks = impl, failed_arms = failed, blocked = blocks, displays = DSP, round3 = r3, package_commits = commits,
       seconds = proc.time()[["elapsed"]] - t0)
}
