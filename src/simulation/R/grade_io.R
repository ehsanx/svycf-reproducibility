# Reading saved PRODUCTION output for the grading pass: run folders, provenance, completeness
# against the FROZEN design, the required-check audit, and every FAILED fit.
#   - Completeness is checked against frozen_cells() -- the expected arms, rungs, populations and replicate IDs of
#     every production cell, fixed in this file -- never against what the runs say about themselves. The
#     production cells.tsv must agree with it; a cell that does not is refused.
#   - Nothing is graded on partial output: an incomplete cell is NOT GRADED, with what is missing and the rerun
#     command (a rerun reproduces the same keyed draws; FAILED fits are never grounds for a rerun).
#   - Failed or missing required checks, and non-finite required values in replicates marked OK, BLOCK the grades
#     they affect.
#   - Every FAILED fit is preserved: the failures of superseded folders are kept with those of the used folders.
#   - The R-d trajectories alpha = 0.20 and 0.35 were stopped after their learner-rate check (G5), so
#     production does not run them and the design here does not expect them. A replicate whose T2 is
#     UNRESOLVED (T2_status, R/t2_refine.R) is neither FAILED nor non-finite-OK: arm_rows() counts it
#     apart, for the grades that use T1 or T2 only.
# Requires R/run_mode.R (parse_ids), R/rng_streams.R (SIM_PROD_ROOT) and the arm constants (R/estimators.R,
# R/tmle.R, R/learner_rb.R, R/controls.R, R/learner_rd.R).

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

# The R-d trajectories STOPPED after G5. Their recorded
# G5 failures are in the execution-check summary; production runs the others (production-cells.tsv --rd-alphas).
GRADE_RD_STOPPED <- c(0.20, 0.35)
grade_rd_stopped <- function(alpha) any(abs(alpha - GRADE_RD_STOPPED) < 1e-9)
grade_rd_run <- function() RD_ALPHAS[!vapply(RD_ALPHAS, grade_rd_stopped, TRUE)]

# Ra-TMLE's T1 and T2 were not independently checked (the U9 gate and the T2 validation cover
# Ra-AIPW, the R-d rungs and AE only). They stay archived in the runs' results.csv (the runner is unchanged) and are
# excluded from the paper's results and displays. No grade uses them: C6 uses the paired psi and se differences.
GRADE_REMAINDER_EXCLUDED_ARMS <- ARM_RA_TMLE
GRADE_R28_NOTE <- "T1/T2 not shown: not independently checked; archived in results.csv, excluded from evidence and displays (ruling R28)"
grade_remainder_excluded <- function(arm) arm %in% GRADE_REMAINDER_EXCLUDED_ARMS

# data.table's reader parses the runners' %.17g doubles exactly and is much faster on production-sized files.
grade_read_csv <- function(p) {
  if (!file.exists(p)) return(data.frame())
  if (requireNamespace("data.table", quietly = TRUE))
    return(data.table::fread(p, data.table = FALSE, na.strings = "NA", showProgress = FALSE))
  utils::read.csv(p, stringsAsFactors = FALSE)
}

# ---- the frozen production design ----------------------------------------------------------------
# Arms: (terminal rungs; R-b many-strata only, A1; R-d fixed-strata m >= 200; at m = 3200 only the
# R-d rungs with the internal oracle), (E3 foil with the R-a controls and the coarsened arm),
# (the AE arm exists iff the frozen kappa_AE record is buildable). Counts:. Rungs: and
# A2 (strict H = 25..400 at J = 4; vanishing (25,10)..(200,80) + (200,800); unequal (50,20) x two allocations),
# (E4: H in {25,...,400}, 2000 draws). Replicate IDs are 1..R (production, RNG contract).
frozen_cells <- function(ae_buildable = TRUE) {
  O <- ARM_ORACLE; RA <- ARM_RA_AIPW; RAT <- ARM_RA_TMLE; RB <- ARM_RB_AIPW; RBT <- ARM_RB_TMLE
  RD <- vapply(grade_rd_run(), rd_arm, "")                 # the stopped rungs are not run in production
  L <- list()
  add <- function(id, pops, R, arms = character(0), rungs = NULL) L[[id]] <<- list(pops = as.integer(pops), R = as.integer(R), arms = arms, rungs = rungs)
  for (H in c(25L, 50L, 100L, 200L, 400L)) add(sprintf("P-ms-H%d", H), 1:400, 5, c(O, RA, RAT, RB, RBT))
  add("P-ms-H800", 1:400, 5, c(O, RA))
  for (H in c(25L, 50L, 100L, 200L, 400L, 800L)) add(sprintf("P-ms-or-H%d", H), 100001:100012, 4000, O)
  for (m in c(50L, 100L)) add(sprintf("P-fs-m%d", m), 1:400, 5, c(O, RA, RAT))
  for (m in c(200L, 400L, 800L)) add(sprintf("P-fs-m%d", m), 1:400, 5, c(O, RA, RAT, RD))
  add("P-fs-m1600", 1:400, 5, c(O, RA, RD))
  add("P-fs-m3200", 1:400, 5, c(O, RD))
  for (m in c(50L, 100L, 200L, 400L, 800L, 1600L)) add(sprintf("P-fs-or-m%d", m), 100001:100012, 4000, O)
  for (J in c(100L, 20L, 8L)) add(sprintf("P-fpc-or-J%d", J), 1:12, 2000, O)
  add("P-fpc-or-J4", 1:12, 10000, O)
  add("P-fpc-fit-J4", 1:160, 5, c(O, RA, RB))
  add("P-strict-m200", 1:160, 5, c(O, RA, RB))
  add("P-calib-J4", 1L, 8000, O)
  add("P-wr-calib", 1L, 8000, O)
  add("P-strict-census", 1:120, 0, rungs = data.frame(H = c(25L, 50L, 100L, 200L, 400L), J = 4L))
  add("P-unequal-census", 1:120, 0, rungs = data.frame(H = 50L, J = 20L, alloc = c("high_on_high", "high_on_low"), stringsAsFactors = FALSE))
  add("P-vanish-census", 1:120, 0, rungs = data.frame(H = c(25L, 50L, 100L, 200L, 200L), J = c(10L, 20L, 40L, 80L, 800L)))
  add("P-e3-rc-m400", 1:400, 5, c(O, RA, ARM_RA_NOCF, ARM_RC_PSUCF, ARM_RC_NOCF, ARM_RC_INDCF, ARM_RC_INDCF_COARSE))
  add("P-e3-vsens-V2", 1:400, 5, c(O, RB))
  add("P-e3-vsens-V5", 1:400, 5, c(O, RB))
  add("P-e3-ctrl-m400", 1:160, 5, c(O, RA, if (isTRUE(ae_buildable)) ARM_AE_AIPW))
  add("P-e4-array", 1L, 2000, rungs = data.frame(H = c(25L, 50L, 100L, 200L, 400L)))
  L
}

rung_key <- function(d) {
  k <- paste0("H", d$H)
  if ("J" %in% names(d)) k <- paste0(k, "_J", d$J)
  if ("alloc" %in% names(d)) k <- paste0(k, "_", d$alloc)
  k
}

# The run folders of one cell: the newest folder of every task (names end in -tNNN-<UTC stamp>; a rerun task
# supersedes), or the newest folder for a one-task cell. Returns list(used, superseded).
cell_run_folders <- function(cell_dir, tasks) {
  d <- list.dirs(cell_dir, recursive = FALSE, full.names = TRUE)
  d <- d[file.exists(file.path(d, "manifest.json"))]
  if (!length(d)) return(list(used = character(0), superseded = character(0)))
  d <- d[order(basename(d), decreasing = TRUE)]
  if (tasks > 1L) {
    tag <- ifelse(grepl("-t[0-9]{3}-[0-9TZ]+$", basename(d)), sub("^.*-(t[0-9]{3})-[0-9TZ]+$", "\\1", basename(d)), NA_character_)
    keep <- !is.na(tag) & !duplicated(tag)
    list(used = d[keep], superseded = d[!keep])
  } else list(used = d[1], superseded = d[-1])
}

run_task_of <- function(folder) {
  b <- basename(folder)
  if (grepl("-t[0-9]{3}-[0-9TZ]+$", b)) as.integer(sub("^.*-t([0-9]{3})-[0-9TZ]+$", "\\1", b)) else 0L
}

# What one production cell must contain: the FROZEN expectation, checked against its cells.tsv line.
cell_expectation <- function(row, frozen) {
  a <- strsplit(row$args, "[[:space:]]+")[[1]]
  val <- function(flag) { i <- which(a == flag); if (length(i)) a[i[1] + 1L] else NA_character_ }
  pops_tsv <- tryCatch(parse_ids(val("--pops"), "--pops"), error = function(e) integer(0))
  R_tsv <- suppressWarnings(as.integer(val("--R"))); if (is.na(R_tsv)) R_tsv <- 0L
  fz <- frozen[[row$id]]
  mismatch <- if (is.null(fz)) "cell not in the frozen production design" else
    c(if (!identical(as.integer(pops_tsv), fz$pops)) sprintf("cells.tsv populations (%d) differ from the frozen design (%d)", length(pops_tsv), length(fz$pops)),
      if (!identical(R_tsv, fz$R)) sprintf("cells.tsv R = %d differs from the frozen R = %d", R_tsv, fz$R))
  list(cell = row$id, script = row$script, pops = if (is.null(fz)) pops_tsv else fz$pops, R = if (is.null(fz)) R_tsv else fz$R,
       reps = if (!is.null(fz) && fz$R > 0L) seq_len(fz$R) else integer(0), arms = fz$arms %||% character(0), rungs = fz$rungs,
       tasks = as.integer(row$tasks), block = as.integer(row$block), traj = val("--traj"),
       mismatch = paste(mismatch, collapse = "; "))
}

# Provenance of one run folder against the frozen production run.
run_provenance <- function(folder, frozen_sha, prod_root = SIM_PROD_ROOT) {
  man <- tryCatch(jsonlite::fromJSON(file.path(folder, "manifest.json"), simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(man)) return(data.frame(folder = basename(folder), mode = NA, root = NA, spec_sha256 = NA, package_commit = NA,
                                      ok = FALSE, problem = "manifest.json unreadable", stringsAsFactors = FALSE))
  spec_sha <- man$spec_freeze$sha256 %||% man$spec_freeze$actual %||% NA_character_
  commit <- man$git$package_commit %||% man$package_commit %||% NA_character_
  root <- man$config$root %||% NA_character_; mode <- man$mode %||% NA_character_
  prob <- c(if (!identical(mode, "prod")) sprintf("mode %s, not prod", mode),
            if (!identical(root, prod_root)) sprintf("root %s, not the production root", root),
            if (!identical(spec_sha, frozen_sha)) "spec SHA-256 differs from the frozen v9.8 record")
  data.frame(folder = basename(folder), mode = mode, root = root, spec_sha256 = spec_sha, package_commit = commit,
             ok = !length(prob), problem = paste(prob, collapse = "; "), stringsAsFactors = FALSE)
}

# Read one production cell: every file the grading pass uses (used folders), and every FAILED fit of every
# folder, superseded ones included (every failure is preserved).
read_prod_cell <- function(in_root, row, frozen_sha, prod_root = SIM_PROD_ROOT, frozen = frozen_cells()) {
  ex <- cell_expectation(row, frozen)
  fol <- cell_run_folders(file.path(in_root, row$id), ex$tasks)
  prov <- if (length(fol$used)) do.call(rbind, lapply(fol$used, run_provenance, frozen_sha = frozen_sha, prod_root = prod_root)) else data.frame()
  rd <- function(name, folders = fol$used) {
    l <- lapply(folders, function(f) { x <- grade_read_csv(file.path(f, name)); if (nrow(x)) { x$run_task <- run_task_of(f); x$run_folder <- basename(f) }; x })
    l <- l[vapply(l, nrow, 0L) > 0L]
    if (!length(l)) return(data.frame())
    cols <- unique(unlist(lapply(l, names)))
    do.call(rbind, lapply(l, function(x) { for (cl in setdiff(cols, names(x))) x[[cl]] <- NA; x[, cols] }))
  }
  mans <- lapply(fol$used, function(f) tryCatch(jsonlite::fromJSON(file.path(f, "manifest.json"), simplifyVector = TRUE), error = function(e) NULL))
  res <- rd("results.csv")
  fail_of <- function(r) if (nrow(r) && "arm_status" %in% names(r)) unique(r[r$arm_status %in% "FAILED", c("pop_id", "rep_id", "arm"), drop = FALSE]) else
    data.frame(pop_id = integer(0), rep_id = integer(0), arm = character(0))
  f_used <- fail_of(res)
  f_sup <- if (length(fol$superseded)) fail_of(rd("results.csv", fol$superseded)) else fail_of(data.frame())
  key <- function(d) paste(d$pop_id, d$rep_id, d$arm)
  failed <- unique(rbind(f_used, f_sup))
  list(expect = ex, folders = fol, provenance = prov, manifests = mans,
       results = res, census = rd("census.csv"), checks = rd("checks.csv"),
       trajectory = rd("census_trajectory.csv"), pair = rd("census_pair.csv"), e4_census = rd("e4_census.csv"),
       failed = failed, failed_superseded_only = sum(!key(f_sup) %in% key(f_used)),
       required_declared = vapply(mans, function(m) as.numeric(m$checks$n_required %||% NA), 0))
}

# Completeness of one read cell against the FROZEN design. Sampled cells: exactly the frozen populations and
# arms, and for every (population, arm) exactly the replicate IDs 1..R, once each (OK or FAILED rows). Census
# cells: exactly the frozen populations, each at exactly the frozen rungs. E4: exactly the frozen rungs, each with
# replicate IDs 1..R. Anything missing, duplicated or unexpected is INCOMPLETE; a design or provenance mismatch is
# REFUSED. Returns list(status, detail, missing_pops, missing_tasks, duplicates).
cell_completeness <- function(x) {
  ex <- x$expect
  out <- function(status, detail, mp = integer(0), mt = integer(0), dup = 0L)
    list(status = status, detail = detail, missing_pops = mp, missing_tasks = mt, duplicates = dup)
  if (nzchar(ex$mismatch)) return(out("REFUSED", paste("design mismatch:", ex$mismatch)))
  if (!length(x$folders$used)) return(out("INCOMPLETE", "no run folder", ex$pops, seq_len(ex$tasks)))
  if (nrow(x$provenance) && !all(x$provenance$ok))
    return(out("REFUSED", paste(unique(x$provenance$problem[!x$provenance$ok]), collapse = "; ")))
  task_of_pop <- function(p) if (ex$tasks > 1L) ((match(p, ex$pops) - 1L) %/% ex$block) + 1L else 1L
  tasks_of <- function(pops) sort(unique(vapply(pops, task_of_pop, 1L)))
  if (ex$script == "e4") {
    r <- x$results
    if (!nrow(r)) return(out("INCOMPLETE", "no E4 draws", mt = 1L))
    want <- ex$rungs$H; have <- sort(unique(r$H))
    miss_H <- setdiff(want, have); extra_H <- setdiff(have, want)
    bad_reps <- want[vapply(want, function(H) { rr <- r$rep_id[r$H == H]; !(length(rr) == ex$R && setequal(rr, seq_len(ex$R))) }, TRUE)]
    dup <- sum(duplicated(r[, c("H", "rep_id")]))
    if (length(miss_H) || length(extra_H) || length(bad_reps) || dup)
      return(out("INCOMPLETE", sprintf("E4: missing rungs %s; unexpected rungs %s; rungs without exactly reps 1..%d: %s; duplicated draws %d",
                                       paste(miss_H, collapse = ","), paste(extra_H, collapse = ","), ex$R, paste(bad_reps, collapse = ","), dup), mt = 1L, dup = dup))
    return(out("COMPLETE", sprintf("%d rungs x reps 1..%d", length(want), ex$R)))
  }
  if (ex$script == "census") {
    tab <- if (identical(ex$traj, "unequal")) x$pair else x$trajectory
    if (!nrow(tab)) return(out("INCOMPLETE", "no census rows", ex$pops, seq_len(ex$tasks)))
    want <- sort(rung_key(ex$rungs))
    tab$rk <- rung_key(tab[, intersect(names(ex$rungs), names(tab)), drop = FALSE])
    mp <- setdiff(ex$pops, tab$pop_id); extra <- setdiff(unique(tab$pop_id), ex$pops)
    bad <- vapply(split(tab$rk, tab$pop_id), function(k) !identical(sort(k), want), TRUE)
    dup <- sum(duplicated(tab[, c("pop_id", "rk")]))
    badp <- as.integer(names(bad)[bad])
    if (length(mp) || length(extra) || length(badp) || dup)
      return(out("INCOMPLETE", sprintf("%d of %d populations missing; %d populations without exactly the frozen rungs (%s); %d duplicated rows; %d unexpected populations",
                                       length(mp), length(ex$pops), length(badp), paste(want, collapse = " "), dup, length(extra)),
                 sort(c(mp, badp)), tasks_of(sort(intersect(c(mp, badp), ex$pops))), dup))
    return(out("COMPLETE", sprintf("%d populations x the %d frozen rungs", length(ex$pops), length(want))))
  }
  r <- x$results
  if (!nrow(r)) return(out("INCOMPLETE", "no results rows", ex$pops, seq_len(ex$tasks)))
  arms_have <- unique(r$arm)
  miss_arm <- setdiff(ex$arms, arms_have); extra_arm <- setdiff(arms_have, ex$arms)
  mp <- setdiff(ex$pops, unique(r$pop_id)); extra_pop <- setdiff(unique(r$pop_id), ex$pops)
  dup <- sum(duplicated(r[, c("pop_id", "rep_id", "arm")]))
  short <- integer(0)
  for (a in intersect(ex$arms, arms_have)) {
    ra <- r[r$arm == a, c("pop_id", "rep_id")]
    sp <- split(ra$rep_id, factor(ra$pop_id, levels = ex$pops))
    ok <- vapply(sp, function(v) length(v) == length(ex$reps) && setequal(v, ex$reps), TRUE)
    short <- c(short, ex$pops[!ok])
  }
  short <- setdiff(unique(short), mp)
  if (length(miss_arm) || length(extra_arm) || length(mp) || length(extra_pop) || length(short) || dup) {
    allp <- sort(unique(c(mp, short)))
    return(out("INCOMPLETE", sprintf("missing arms [%s]; unexpected arms [%s]; %d of %d populations missing; %d populations without exactly reps 1..%d for every arm; %d duplicated rows; %d unexpected populations",
                                     paste(miss_arm, collapse = ", "), paste(extra_arm, collapse = ", "), length(mp), length(ex$pops), length(short), ex$R, dup, length(extra_pop)),
               allp, if (length(miss_arm)) seq_len(ex$tasks) else tasks_of(allp), dup))
  }
  out("COMPLETE", sprintf("%d populations x reps 1..%d x %d frozen arms", length(ex$pops), ex$R, length(ex$arms)))
}

# The required-check audit of one cell. Returns the grades it blocks:
#   all    reasons that block every grade of the cell;
#   arms   named list: reasons that block the grades using that arm;
#   rungs  named list (E4, by H): reasons that block that rung.
# Blocks: (1) a run whose checks.csv holds fewer required rows than its manifest declares (missing checks);
# (2) any required row that did not pass (a row of arm "all", or with no arm, blocks the whole cell);
# (3) missing required families: U1 for every replicate and U8 for every (replicate, arm) marked OK (sampled cells);
# the census FPC algebra at every population and rung (census cells); U1, U8 and finite outputs at every E4 rung.
required_check_audit <- function(x) {
  ex <- x$expect; ck <- x$checks
  A <- list(all = character(0), arms = list(), rungs = list())
  add_arm <- function(a, why) A$arms[[a]] <<- unique(c(A$arms[[a]], why))
  if (!nrow(ck) || !all(c("required", "pass", "check") %in% names(ck))) { A$all <- "no required checks recorded"; return(A) }
  req <- as.logical(ck$required) %in% TRUE
  n_req <- tapply(req, factor(ck$run_folder, levels = basename(x$folders$used)), sum)
  n_req[is.na(n_req)] <- 0                               # a folder with no checks.csv holds no required rows
  decl <- x$required_declared
  short <- which(is.finite(decl) & as.numeric(n_req) < decl)
  if (length(short)) A$all <- c(A$all, sprintf("missing required checks: %s holds %d of the %d required rows its manifest declares",
                                             basename(x$folders$used)[short], as.integer(n_req[short]), as.integer(decl[short])))
  rq <- ck[req, , drop = FALSE]
  bad <- rq[!(as.logical(rq$pass) %in% TRUE), , drop = FALSE]
  if (nrow(bad)) {
    lab <- function(b) sprintf("failed required check: %s (%d row%s)", b$check[1], nrow(b), if (nrow(b) > 1) "s" else "")
    if (ex$script == "e4") {
      for (H in unique(bad$H)) { b <- bad[bad$H %in% H, ]; if (is.na(H)) A$all <- c(A$all, lab(b)) else A$rungs[[as.character(H)]] <- c(A$rungs[[as.character(H)]], lab(b)) }
    } else if (ex$script == "census" || !"arm" %in% names(bad)) {
      A$all <- c(A$all, vapply(split(bad, bad$check), lab, ""))
    } else {
      for (a in unique(bad$arm)) { b <- bad[bad$arm %in% a, ]
        if (is.na(a) || a %in% c("all", "")) A$all <- c(A$all, vapply(split(b, b$check), lab, "")) else for (bb in split(b, b$check)) add_arm(a, lab(bb)) }
    }
  }
  if (ex$script == "slice" && nrow(x$results)) {
    r <- x$results; okr <- if ("arm_status" %in% names(r)) r[r$arm_status %in% "OK", ] else r
    u1 <- rq[startsWith(rq$check, "U1"), ]; u8 <- rq[startsWith(rq$check, "U8"), ]
    rep_keys <- unique(paste(r$pop_id, r$rep_id))
    miss_u1 <- setdiff(rep_keys, unique(paste(u1$pop_id, u1$rep_id)))
    if (length(miss_u1)) A$all <- c(A$all, sprintf("missing required checks: no U1 rows for %d replicate(s)", length(miss_u1)))
    have8 <- unique(paste(u8$pop_id, u8$rep_id, u8$arm))
    mk <- paste(okr$pop_id, okr$rep_id, okr$arm)
    for (a in unique(okr$arm)) { m <- sum(!mk[okr$arm == a] %in% have8)
      if (m) add_arm(a, sprintf("missing required checks: no U8 rows for %d replicate(s) marked OK", m)) }
  }
  if (ex$script == "census") {
    fa <- rq[startsWith(rq$check, "census FPC algebra"), ]
    per <- table(factor(fa$pop_id, levels = ex$pops))
    n_r <- nrow(ex$rungs)
    if (any(per != n_r)) A$all <- c(A$all, sprintf("missing required checks: %d population(s) without the census FPC algebra at all %d frozen rungs", sum(per != n_r), n_r))
  }
  if (ex$script == "e4") for (H in ex$rungs$H) {
    fam <- c("U1", "U8", "finite outputs")
    have <- vapply(fam, function(f) any(rq$H %in% H & startsWith(rq$check, f)), TRUE)
    if (!all(have)) A$rungs[[as.character(H)]] <- c(A$rungs[[as.character(H)]], sprintf("missing required checks at H = %d: %s", H, paste(fam[!have], collapse = ", ")))
  }
  A
}

# Rows of one arm usable by a statistic: arm_status OK and finite in the needed columns. FAILED rows are counted;
# OK rows with a non-finite needed value are counted separately (they BLOCK the grade, never silently dropped).
# When a needed column is T1 or T2, an OK row whose T2 is UNRESOLVED (T2_status) is counted as unresolved,
# not as non-finite: its grade is UNGRADABLE with the survivors' statistic descriptive (R/grade_run.R).
T2_COLS <- c("T1", "T2", "T1_over_sigma_wor", "T2_over_sigma_wor")
arm_rows <- function(r, arm, cols) {
  a <- r[r$arm == arm, , drop = FALSE]
  st <- if ("arm_status" %in% names(a)) a$arm_status else rep("OK", nrow(a))
  ok <- st %in% "OK"
  unres <- if (any(cols %in% T2_COLS) && "T2_status" %in% names(a)) ok & a$T2_status %in% "unresolved" else rep(FALSE, nrow(a))
  fin <- ok & !unres
  for (cl in cols) fin <- fin & is.finite(suppressWarnings(as.numeric(a[[cl]])))
  list(rows = a[fin, , drop = FALSE], failed = sum(st %in% "FAILED"), nonfinite_ok = sum(ok & !unres & !fin), unresolved = sum(unres),
       total = nrow(a))
}
