# THE GRADING PASS: post-processing of the saved PRODUCTION
# output into classifications. No population, sample or fit is made; the production output is only read.
#
# It runs through the engine's locked wrapper (scripts/simulation/submit.sh grade does this), so the
# locked library is used read-only and the production folder is never modified:
#   Rscript src/simulation/env/run.R src/simulation/runs/run_grade.R \
#     [--in <production folder>] [--round3 <preparation folder>] [--assembly-roots <execution-check folder>] \
#     [--cells config/simulation/production-cells.tsv] [--out-root <folder>]
# (the defaults are folders under the engine folder). B_boot = 10,000 and the grading-pass label
# "v9.8-grade-1" are frozen: there is no option to change either.
#
# Writes <out-root>/grade-<UTC>/: GRADES.md, grades.csv (every graded quantity: estimate, 95% interval, region,
# disposition, basis), claims.csv (one composed disposition per joint claim cell),
# population_indexed_summary.csv (design-conditional rows), completeness.csv (against the frozen design),
# blocked.csv (the required-check audit), failed_arms.csv (every FAILED fit, superseded runs included),
# provenance.csv, provenance_checks.csv, required_checks.csv, round3_inputs.csv, displays/*.csv,
# grading_manifest.json.
# Any grade using an arm with a FAILED fit is UNGRADABLE (fitting
# failures) with the survivors' statistic labelled descriptive; every failure is preserved; nothing is retried
# with new seeds. Failed or missing required checks, and non-finite required values marked OK, BLOCK the grades
# they affect.
# Exit 0: graded, every cell complete, nothing blocked, every provenance check met. Exit 3: some cell incomplete
# (NOT GRADED, rerun listed) or some grade BLOCKED. Exit 4: a provenance check failed (the output's identities do not
# match). Exit 5: a grading error in some cell. The report is written in every case.

via_wrapper <- endsWith(getOption("sim.engine.script", ""), "/runs/run_grade.R")
if (!via_wrapper) stop("run through a locked wrapper: Rscript <engine>/env/run.R <sim-grade>/simulation/engine/runs/run_grade.R")
t_start <- proc.time()[["elapsed"]]
stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
code_engine <- normalizePath(file.path(dirname(getOption("sim.engine.script")), ".."), winslash = "/", mustWork = TRUE)
code_repo <- normalizePath(file.path(code_engine, "..", ".."), winslash = "/", mustWork = TRUE)
prod_engine <- getOption("sim.engine.dir")
args <- getOption("sim.engine.args")
arg_of <- function(name, default) { i <- which(args == name); if (length(i) && length(args) > i[1]) args[i[1] + 1L] else default }
for (f in c("rng_streams.R", "run_mode.R", "learner_rd.R", "estimators.R", "tmle.R", "learner_rb.R", "controls.R",
            "grading.R", "grade_plan.R", "grade_io.R", "grade_run.R"))
  source(file.path(code_engine, "R", f), local = TRUE)

in_root <- arg_of("--in", file.path(prod_engine, "runs", "out", "production"))
round3_root <- arg_of("--round3", file.path(prod_engine, "runs", "out", "round3"))
# The execution-check summary may live in its own folder (--assembly-roots, an assemble/ subfolder);
# the preparation outputs and kappa_AE are under --round3. The newest execution_gates.json across them is read.
assembly_roots <- strsplit(arg_of("--assembly-roots", file.path(prod_engine, "runs", "out", "round3d")), ",", fixed = TRUE)[[1]]
cells_tsv <- arg_of("--cells", file.path(prod_engine, "arc", "production", "cells.tsv"))
out_dir <- file.path(arg_of("--out-root", file.path(code_engine, "runs", "out", "grading")), paste0("grade-", stamp))
spec <- spec_freeze_check(code_repo)
if (!isTRUE(spec$ok)) stop(sprintf("the design specification hash does not match the recorded one (SHA-256 %s, record %s)", spec$actual, spec$expected), call. = FALSE)
cells <- utils::read.delim(cells_tsv, stringsAsFactors = FALSE, colClasses = "character")
cat(sprintf("[grade] spec v9.8 verified (%s); %d production cells from %s; output %s; B_boot = %d, pass %s\n",
            substr(spec$actual, 1, 12), nrow(cells), in_root, out_dir, GRADE_B_BOOT, GRADE_PASS))

res <- grade_run(in_root, round3_root, cells, out_dir, frozen_sha = spec$actual, assembly_roots = assembly_roots)

# ---- the report ---------------------------------------------------------------------------------------
cl <- res$claims; gr <- res$grades; cmp <- res$completeness; pc <- res$provenance_checks
incomplete <- cmp[cmp$status != "COMPLETE", , drop = FALSE]
rerun <- if (nrow(incomplete)) sprintf("bash simulation/engine/arc/submit.sh --round production --only %s --approve \"<name, date, core-h>\"",
                                       paste(incomplete$cell, collapse = ",")) else ""
blocked_rows <- if (nrow(gr)) gr[startsWith(gr$disposition, "BLOCKED"), , drop = FALSE] else gr
r25_rows <- if (nrow(gr)) gr[startsWith(gr$disposition, "UNGRADABLE (fitting failures"), , drop = FALSE] else gr
tab <- function(d, cols) if (!nrow(d)) "(none)" else c(paste0("| ", paste(cols, collapse = " | "), " |"), paste0("|", strrep("---|", length(cols))),
  apply(d[, cols, drop = FALSE], 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
fmt <- function(x) ifelse(is.na(x), "", formatC(as.numeric(x), digits = 4, format = "fg", flag = "#"))
g_show <- if (nrow(gr)) transform(gr, estimate = fmt(estimate), interval = ifelse(is.na(ci_lo), "", sprintf("[%s, %s]", fmt(ci_lo), fmt(ci_hi))),
                                   region = ifelse(is.na(region_lo), "", sprintf("[%s, %s]", region_lo, region_hi)),
                                   population = ifelse(is.na(population), "", population)) else gr
sci <- if (nrow(g_show)) g_show[g_show$role %in% c("scientific", "boundary", "stress") & g_show$population == "", , drop = FALSE] else g_show
pcs <- if (nrow(g_show)) g_show[g_show$role == "power control", , drop = FALSE] else g_show
md <- c(sprintf("# Grading pass (spec v9.8): %s", stamp), "",
  sprintf("> Post-processing of the saved production output by the frozen §6 plan: §6.1's operator on the 95%% percentile interval of B_boot = %d keyed resamples (level 5, pass `%s`); joint cells resample population blocks only, every block keeping all of its draws, populations equally weighted (R6); design-conditional cells resample their draws and report population-indexed dispositions (R20). Spec SHA-256 %s. Every frozen cell is reported regardless of outcome (§6 publication rule).",
          GRADE_B_BOOT, GRADE_PASS, spec$actual), "",
  "## Provenance", "", tab(transform(pc, ok = ifelse(ok, "yes", "**NO**")), c("check", "ok", "detail")), "",
  "## Completeness", "", tab(cmp, c("cell", "status", "detail", "runs", "superseded")), "",
  if (nrow(incomplete)) c("**Incomplete or refused cells are NOT GRADED** (nothing is graded on partial output; completeness is judged against the frozen design's arms, rungs, populations and replicate IDs). Missing output may be rerun -- the keyed RNG reproduces the same draws and a rerun task supersedes its old folder, while every FAILED fit of any folder stays recorded (R25); FAILED fits are never grounds for a rerun:", "", "```", rerun, "```", "") else "Every production cell is complete against the frozen design.", "",
  "## Blocked grades (failed or missing required checks; non-finite required values marked OK)", "",
  if (nrow(res$blocked)) tab(res$blocked, c("cell", "block")) else "None.", "",
  sprintf("%d graded quantities are BLOCKED.", nrow(blocked_rows)), "",
  "## Fitting failures (author's grading ruling R25)", "",
  "Any grade using an arm with a FAILED fit in that cell is UNGRADABLE (fitting failures); the survivors' statistic is shown, labelled descriptive. Every failure is preserved, superseded runs included; nothing was retried with new seeds.", "",
  tab(if (nrow(res$failed_arms)) res$failed_arms[res$failed_arms$failed > 0, , drop = FALSE] else res$failed_arms, c("cell", "arm", "replicates", "failed", "failed_in_superseded_runs_only")), "",
  sprintf("%d graded quantities are UNGRADABLE (fitting failures).", nrow(r25_rows)), "",
  "## §8 claims: one disposition per graded cell (§6.1a; design-conditional cells are not composed, R20)", "",
  tab(cl, c("claim", "cell", "arm", "metrics", "disposition")), "",
  "## Graded quantities (scientific, boundary and stress rows)", "",
  tab(sci, c("claim", "cell", "arm", "statistic", "estimate", "interval", "region", "disposition", "basis", "failed_replicates")), "",
  "## Power controls (E3 panel; predicted to fail; graded only where G6 limb (a) passed)", "",
  tab(pcs, c("arm", "statistic", "estimate", "interval", "region", "disposition", "premise")), "",
  "## Population-indexed rows (design-conditional: C1 Z_oracle, C2 oracle limb, C3a, C3b and its premise, C7)", "",
  tab(res$indexed, c("claim", "role", "cell", "statistic", "populations", "PASS", "RED_FLAG", "UNGRADABLE", "BLOCKED", "NOT_GRADED")), "",
  "Every row, with its estimate and interval: `grades.csv`.", "",
  "## Implementation checks in production (required rows of every run)", "",
  tab(res$required_checks, c("cell", "required_checks", "declared", "failed", "first_failure")), "",
  "## Displays (no margin, no disposition)", "", paste0("- `displays/", names(res$displays), ".csv`"), "",
  "Ra-TMLE's T1 and T2 appear in no grade and no display (author's ruling R28): they were not independently checked, and the runs' `results.csv` keeps them as archived diagnostics.", "",
  sprintf("Grading wall time %.1f min.", res$seconds / 60))
writeLines(md, file.path(out_dir, "GRADES.md"), useBytes = TRUE)
sha <- function(p) if (file.exists(p)) digest::digest(file = p, algo = "sha256") else NA_character_
code_files <- c(file.path("R", c("rng_streams.R", "run_mode.R", "grading.R", "grade_plan.R", "grade_io.R", "grade_run.R")), "runs/run_grade.R")
jsonlite::write_json(list(created_utc = stamp, spec = spec, grading_pass = GRADE_PASS, B_boot = GRADE_B_BOOT, level = GRADE_LEVEL,
                          production_package_commits = res$package_commits,
                          grading_package_commit = if (file.exists(file.path(code_repo, "PACKAGE-COMMIT.txt"))) readLines(file.path(code_repo, "PACKAGE-COMMIT.txt"), warn = FALSE)[1] else NA,
                          inputs = list(production = in_root, round3 = round3_root, cells = cells_tsv, cells_sha256 = sha(cells_tsv)),
                          round3_files = res$round3$files, code_sha256 = stats::setNames(as.list(vapply(file.path(code_engine, code_files), sha, "")), code_files),
                          runtime = list(R = R.version.string, platform = R.version$platform), seconds = res$seconds),
                     file.path(out_dir, "grading_manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = I(17), null = "null", na = "null")
status <- if (any(cmp$status == "GRADING ERROR")) 5L else if (!all(pc$ok)) 4L else if (nrow(incomplete) || nrow(blocked_rows)) 3L else 0L
cat(sprintf("[grade] wrote %s in %.1f min; exit %d\n", out_dir, (proc.time()[["elapsed"]] - t_start) / 60, status))
quit(save = "no", status = status)
