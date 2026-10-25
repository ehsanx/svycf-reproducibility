# ---------------------------------------------------------------------------
# 01-count-audit.R -- the SPA per-facility client-count audit, behind a frozen
# feasibility gate
# (config/application/gate.yaml).
#
# Exactly one mode flag is required. There is no default mode.
#
#   --selftest
#       Synthetic data in a temporary directory only. Exercises the counting,
#       binning, zero-facility handling, gate parsing and evaluation, the code
#       preflight (packaged copy; git checkout too where git is installed), the
#       no-data-before-preflight guarantee, the structure-mode output filter and
#       non-overwrite. Opens no real data file.
#
#   --structure --survey <S>
#       Schema and linkage checks on the real files. Prints variable presence,
#       class and label, value-label checks, and TRUE/FALSE linkage results
#       ONLY: no row counts, frequency tables, unmatched-ID counts, identifiers
#       or importer diagnostics. Allowed before the gate is frozen.
#
#   --count --survey <S>
#       The audit. Opens NO data file unless (1) the gate parses, is status:
#       FROZEN, has every field set and lists <S>, and (2) the gate and this
#       script are unchanged: in a git checkout, tracked, in
#       HEAD and identical to HEAD in index and working tree; in a packaged copy
#       without git, listed in its code manifest with a matching
#       md5. Then writes aggregate metrics and the
#       PASS/FAIL result to the results folder, never overwriting and never writing a
#       facility identifier, and appends one line to RUNS.md there.
#
# Surveys: MW_SPA_2013-14, HT_SPA_2013, HT_SPA_2017-18.
#
# DATA POLICY (DHS Terms of Use): record-level data never leave the private data
# folder. Outputs are aggregates.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(haven)
  library(yaml)
})

# --- configuration -----------------------------------------------------------
# File stems (zip names without .zip) per survey and data structure. Recode = SR,
# raw = SP (dhsprogram.com/data/File-Types-and-Names.cfm).
FILES <- list(
  "MW_SPA_2013-14" = list(
    recode = c(FC = "MWFC6JDTSR", AN = "MWAN6IDTSR", FP = "MWFP6IDTSR", SC = "MWSC6IDTSR"),
    raw    = c(FC = "MWFC6KDTSP", AN = "MWAN6KDTSP", FP = "MWFP6KDTSP", SC = "MWSC6KDTSP")),
  "HT_SPA_2013" = list(
    recode = c(FC = "HTFC6ADTSR", AN = "HTAN6ADTSR", FP = "HTFP6ADTSR", SC = "HTSC6ADTSR"),
    raw    = c(FC = "HTFC6BDTSP", AN = "HTAN6BDTSP", FP = "HTFP6BDTSP", SC = "HTSC6BDTSP")),
  "HT_SPA_2017-18" = list(
    recode = c(FC = "HTFC7ADTSR", AN = "HTAN7ADTSR", FP = "HTFP7ADTSR", SC = "HTSC7ADTSR"),
    raw    = c(FC = "HTFC7ADTSP", AN = "HTAN7ADTSP", FP = "HTFP7ADTSP", SC = "HTSC7ADTSP"))
)
SURVEYS  <- names(FILES)
SERVICES <- c("AN", "FP", "SC")

# Variable names (matched case-insensitively), CONFIRMED on all three surveys with
# --structure and metadata-only reads. Value labels: code 1 = facility
# completed (v010a/result), observation complete (c020a), exit complete (c020b);
# raw c101 = client agreed to observation, raw x?100 = client agreed to the exit
# interview (agreement, not completion). This block is part of the script
# that count mode verifies.
#   fc_counts: facility-file per-service counts of clients who completed observation
#   and exit (recode v071 ANC, v072 sick child, v073 FP -- labels identical in all
#   three surveys). NA = no cross-check.
#   region/type may be a vector named by country code where names differ by survey.
VARS <- list(
  recode = list(
    fac_id = "inv_id", fac_done = "v010a", region = "v001", type = "v007",
    cli_id = "inv_id",
    cli_obs  = c(AN = "c020a", FP = "c020a", SC = "c020a"),
    cli_exit = c(AN = "c020b", FP = "c020b", SC = "c020b"),
    fc_counts = c(AN = "v071", FP = "v073", SC = "v072"),
    fc_count_candidates = c("v071", "v072", "v073")),
  raw = list(
    fac_id = "facil", fac_done = "result", region = c(MW = "region", HT = "depart"), type = "factype",
    cli_id = "facil",
    cli_obs  = c(AN = "c101", FP = "c101", SC = "c101"),
    cli_exit = c(AN = "xa100", FP = "xf100", SC = "xc100"),
    fc_counts = c(AN = NA_character_, FP = NA_character_, SC = NA_character_),
    fc_count_candidates = character(0))
)
DONE_CODE <- 1   # DHS: 1 = completed (facility, observation) / agrees (exit consent)

# Facility-type groups for the typegroup strata, keyed by country code, mapping the
# recode v007 codes (from the value labels). Malawi: 1-4 central / district /
# rural-community / other hospital; 5 health centre, 6 maternity; 7 dispensary,
# 8 clinic, 9 health post. Haiti: 1-4 university / regional-department / community
# referral / other hospital; 5-6 health centre with / without beds; 7 dispensary.
TYPE_GROUPS <- list(
  MW = c("1" = "hospital", "2" = "hospital", "3" = "hospital", "4" = "hospital",
         "5" = "health centre", "6" = "health centre",
         "7" = "lower-level", "8" = "lower-level", "9" = "lower-level"),
  HT = c("1" = "hospital", "2" = "hospital", "3" = "hospital", "4" = "hospital",
         "5" = "health centre", "6" = "health centre", "7" = "lower-level"))
STRATA <- c("region_x_type", "region", "type", "typegroup", "region_x_typegroup", "none")

# Unmapped or missing codes give NA, which stratum_of() refuses for analysed facilities.
type_group_of <- function(codes, survey) {
  map <- TYPE_GROUPS[[substr(survey, 1, 2)]]
  if (is.null(map)) fail("no facility-type groups configured for ", survey)
  unname(map[as.character(as.vector(codes))])
}

PLACEHOLDERS <- c("", "tbd", "todo", "null", "na", "?", "xx", "xxx", "fixme", "none yet")
THRESHOLD_NAMES <- c("min_mean_nk", "max_share_nk_0", "max_share_nk_le_1", "max_share_nk_le_2",
                     "min_facilities_per_stratum", "stage2_c", "min_share_nk_gt_c")

# --- small utilities -----------------------------------------------------------
fail <- function(...) stop(paste0(...), call. = FALSE)

# A frame-variable name for this survey: a plain string, or a vector named by
# country code (the first two letters of the survey key).
var_for <- function(x, survey) {
  if (is.null(names(x))) return(x)
  cc <- substr(survey, 1, 2)
  if (!cc %in% names(x)) fail("no variable configured for country ", cc)
  x[[cc]]
}

# A YAML list of distinct allowed strings, with no null, empty or nested entries.
valid_set <- function(x, allowed) {
  v <- unlist(x)
  !is.null(x) && length(x) >= 1 && length(v) == length(x) && is.character(v) &&
    !anyNA(v) && all(v %in% allowed) && !anyDuplicated(v)
}

is_placeholder <- function(x) {
  is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x)) ||
    (is.character(x) && length(x) == 1 && tolower(trimws(x)) %in% PLACEHOLDERS)
}

# Read a .dta with lower-cased names, no importer chatter, optionally only some
# columns / zero rows (metadata only). Errors are re-raised WITHOUT the
# importer's message, which can quote row positions.
read_dta_lower <- function(path, cols = NULL, n_max = Inf) {
  tryCatch({
    meta <- suppressWarnings(suppressMessages(haven::read_dta(path, n_max = 0)))
    actual <- names(meta)
    if (!is.null(cols)) {
      hit <- actual[match(tolower(cols), tolower(actual))]
      hit <- hit[!is.na(hit)]
      d <- suppressWarnings(suppressMessages(
        haven::read_dta(path, col_select = tidyselect::all_of(hit), n_max = n_max)))
    } else {
      d <- if (n_max == 0) meta else
        suppressWarnings(suppressMessages(haven::read_dta(path, n_max = n_max)))
    }
    names(d) <- tolower(names(d))
    d
  }, error = function(e) fail("read error in ", basename(path), " (importer details suppressed)"))
}

find_dta <- function(data_root, survey, structure, module) {
  stem <- FILES[[survey]][[structure]][[module]]
  dir  <- file.path(data_root, survey, structure, stem)
  hits <- list.files(dir, pattern = "\\.dta$", ignore.case = TRUE, recursive = TRUE, full.names = TRUE)
  if (length(hits) != 1) return(NA_character_)
  hits
}

# --- gate --------------------------------------------------------------------------
read_gate_yaml <- function(path) {
  if (!file.exists(path)) fail("gate file not found: ", path)
  txt <- readLines(path, warn = FALSE, encoding = "UTF-8")
  fm  <- which(trimws(txt) == "---")
  if (length(fm) < 2 || fm[1] != 1) fail("gate has no YAML front matter delimited by --- lines")
  yaml::yaml.load(paste(txt[(fm[1] + 1):(fm[2] - 1)], collapse = "\n"))
}

# Returns a validated gate list, or stops naming every problem. Never partial.
validate_gate <- function(g) {
  bad <- character(0)
  add <- function(...) bad <<- c(bad, paste0(...))
  if (!identical(g$status, "FROZEN")) add("status must be FROZEN (is ", format(g$status %||% "missing"), ")")
  if (is_placeholder(g$frozen_by)) add("frozen_by unset")
  if (is_placeholder(g$frozen_on) || is.na(as.Date(as.character(g$frozen_on), optional = TRUE)))
    add("frozen_on unset or not YYYY-MM-DD")
  if (!valid_set(g$surveys, SURVEYS)) add("surveys must be a non-empty list of known surveys")
  if (is_placeholder(g$count_source) || !g$count_source %in% c("recode", "raw")) add("count_source must be recode|raw")
  if (is_placeholder(g$client_record) || !g$client_record %in% c("obs_and_exit", "obs_only"))
    add("client_record must be obs_and_exit|obs_only")
  if (!valid_set(g$services, SERVICES)) add("services must be a non-empty subset of AN, FP, SC")
  if (is_placeholder(g$zero_record_facilities) || !g$zero_record_facilities %in% c("keep", "drop"))
    add("zero_record_facilities must be keep|drop")
  if (is_placeholder(g$strata) || !g$strata %in% STRATA)
    add("strata must be one of ", paste(STRATA, collapse = "|"))
  th <- g$thresholds
  if (!is.list(th)) { add("thresholds block missing"); th <- list() }
  is_off <- function(x) identical(x, "off")
  for (nm in THRESHOLD_NAMES) {
    x <- th[[nm]]
    if (is_off(x)) next
    if (isFALSE(x)) { add(nm, ': write "off" in quotes (unquoted off is YAML for false)'); next }
    if (is_placeholder(x) || !is.numeric(x) || length(x) != 1 || !is.finite(x)) { add(nm, " unset (number or \"off\")"); next }
    if (grepl("share", nm) && (x < 0 || x > 1)) add(nm, " must be in [0, 1]")
    if (nm %in% c("min_facilities_per_stratum", "stage2_c") && (x < 1 || x != round(x))) add(nm, " must be an integer >= 1")
    if (nm == "min_mean_nk" && x < 0) add(nm, " must be >= 0")
  }
  if (!is_off(th$min_share_nk_gt_c) && is_off(th$stage2_c)) add("min_share_nk_gt_c needs stage2_c")
  graded <- setdiff(THRESHOLD_NAMES, "stage2_c")
  if (all(vapply(graded, function(nm) is_off(th[[nm]]), logical(1)))) add("every threshold is off: not a gate")
  b <- unlist(g$reporting_bins)
  if (is_placeholder(g$reporting_bins) || !is.numeric(b) || length(b) < 2 || b[1] != 0 || any(diff(b) <= 0) ||
      any(b != round(b)))
    add("reporting_bins must be increasing integer lower edges starting at 0 (at least two)")
  if (length(bad)) fail("GATE NOT USABLE FOR COUNT MODE:\n  - ", paste(bad, collapse = "\n  - "))
  g$surveys <- unlist(g$surveys); g$services <- unlist(g$services); g$reporting_bins <- b
  g
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# --- git preflight -------------------------------------------------------------------
git_run <- function(repo, args) {
  out <- suppressWarnings(system2("git", c("-C", shQuote(repo), args), stdout = TRUE, stderr = TRUE))
  list(status = attr(out, "status") %||% 0L, out = out)
}

# Every path must be tracked, present in HEAD, and identical to HEAD in the index
# and the working tree. Returns provenance hashes; stops otherwise.
preflight_git <- function(repo, rel_paths) {
  prob <- character(0)
  for (p in rel_paths) {
    if (git_run(repo, c("ls-files", "--error-unmatch", "--", p))$status != 0) { prob <- c(prob, paste0(p, ": not tracked")); next }
    if (git_run(repo, c("cat-file", "-e", paste0("HEAD:", p)))$status != 0) { prob <- c(prob, paste0(p, ": not in HEAD")); next }
    if (git_run(repo, c("diff", "--quiet", "--cached", "HEAD", "--", p))$status != 0) prob <- c(prob, paste0(p, ": index differs from HEAD"))
    if (git_run(repo, c("diff", "--quiet", "--", p))$status != 0) prob <- c(prob, paste0(p, ": working tree differs from index"))
  }
  if (length(prob)) fail("PREFLIGHT FAILED (no data file was opened):\n  - ", paste(prob, collapse = "\n  - "))
  one <- function(args) git_run(repo, args)$out[1]
  c(head = one(c("rev-parse", "HEAD")),
    setNames(vapply(rel_paths, function(p) one(c("rev-parse", paste0("HEAD:", p))), ""), paste0("blob:", rel_paths)),
    setNames(vapply(rel_paths, function(p) one(c("log", "-1", "--format=%H", "--", p)), ""), paste0("last_commit:", rel_paths)))
}

# A packaged copy without git can hold CODE-MANIFEST.tsv at its root,
# written from a git checkout: line 1 "# head<TAB><commit>", then
# path, md5 (of the shipped file), blob and last_commit. Every path must be listed and match its md5;
# the provenance returned is the committed one, exactly as preflight_git returns it. A git checkout
# always uses preflight_git, even if a manifest is lying around.
CODE_MANIFEST <- "CODE-MANIFEST.tsv"
# The reproduction release (no git) carries config/application/code-manifest.tsv instead: line 1 "# release<TAB><id>",
# then path and md5 of every file the preflight names. The provenance it returns is the release identity and the
# md5 of each file -- not a commit, which the release does not have.
RELEASE_MANIFEST <- file.path("config", "application", "code-manifest.tsv")
preflight_release <- function(repo, rel_paths) {
  l <- readLines(file.path(repo, RELEASE_MANIFEST), warn = FALSE)
  if (length(l) < 2 || !grepl("^# release\t[A-Za-z0-9._:+-]+$", l[1]) || l[2] != "path\tmd5")
    fail("PREFLIGHT FAILED (no data file was opened): ", RELEASE_MANIFEST, " is malformed")
  m <- utils::read.delim(text = l[-1], colClasses = "character", quote = "")
  prob <- character(0)
  for (p in rel_paths) {
    i <- match(p, m$path)
    if (is.na(i)) { prob <- c(prob, paste0(p, ": not in ", RELEASE_MANIFEST)); next }
    f <- file.path(repo, p)
    if (!file.exists(f)) { prob <- c(prob, paste0(p, ": missing")); next }
    if (!identical(unname(tools::md5sum(f)), m$md5[i])) prob <- c(prob, paste0(p, ": differs from ", RELEASE_MANIFEST, " (md5)"))
  }
  if (length(prob)) fail("PREFLIGHT FAILED (no data file was opened):\n  - ", paste(prob, collapse = "\n  - "))
  i <- match(rel_paths, m$path)
  c(head = paste0("release:", sub("^# release\t", "", l[1])),
    setNames(paste0("md5:", m$md5[i]), paste0("blob:", rel_paths)), setNames(rep("release", length(rel_paths)), paste0("last_commit:", rel_paths)))
}
preflight_code <- function(repo, rel_paths) {
  if (file.exists(file.path(repo, ".git"))) return(preflight_git(repo, rel_paths))
  if (file.exists(file.path(repo, RELEASE_MANIFEST))) return(preflight_release(repo, rel_paths))
  mf <- file.path(repo, CODE_MANIFEST)
  if (!file.exists(mf))
    fail("PREFLIGHT FAILED (no data file was opened): ", repo, " is neither a git checkout nor a packaged copy with ", CODE_MANIFEST)
  l <- readLines(mf, warn = FALSE)
  if (length(l) < 2 || !grepl("^# head\t[0-9a-f]{40}$", l[1]) || l[2] != "path\tmd5\tblob\tlast_commit")
    fail("PREFLIGHT FAILED (no data file was opened): ", CODE_MANIFEST, " is malformed")
  m <- utils::read.delim(text = l[-1], colClasses = "character", quote = "")
  prob <- character(0)
  for (p in rel_paths) {
    i <- match(p, m$path)
    if (is.na(i)) { prob <- c(prob, paste0(p, ": not in ", CODE_MANIFEST)); next }
    f <- file.path(repo, p)
    if (!file.exists(f)) { prob <- c(prob, paste0(p, ": missing")); next }
    if (!identical(unname(tools::md5sum(f)), m$md5[i])) prob <- c(prob, paste0(p, ": differs from ", CODE_MANIFEST, " (md5)"))
  }
  if (length(prob)) fail("PREFLIGHT FAILED (no data file was opened):\n  - ", paste(prob, collapse = "\n  - "))
  i <- match(rel_paths, m$path)
  c(head = sub("^# head\t", "", l[1]),
    setNames(m$blob[i], paste0("blob:", rel_paths)), setNames(m$last_commit[i], paste0("last_commit:", rel_paths)))
}

# --- counting (pure functions; no I/O) ---------------------------------------------
# fac: data.frame(id, done, region, type); cli: named list (by service) of
# data.frame(id, obs, exit). Returns list(nk = data.frame(id, region, type, nk),
# n_unlinked = client records whose facility is not a completed facility).
compute_nk <- function(fac, cli, services, client_record, zero_rule) {
  fac <- fac[!is.na(fac$done) & fac$done == DONE_CODE, intersect(c("id", "region", "type", "typegroup"), names(fac))]
  keep_rec <- function(d) {
    ok <- !is.na(d$obs) & d$obs == DONE_CODE
    if (client_record == "obs_and_exit") ok <- ok & !is.na(d$exit) & d$exit == DONE_CODE
    d$id[ok]
  }
  ids <- unlist(lapply(services, function(s) keep_rec(cli[[s]])), use.names = FALSE)
  linked <- ids %in% fac$id
  tab <- table(factor(ids[linked], levels = fac$id))
  fac$nk <- as.integer(tab[as.character(fac$id)])
  census <- fac                                          # completed facilities, before the zero rule
  if (zero_rule == "drop") fac <- fac[fac$nk > 0, , drop = FALSE]
  list(nk = fac, census = census, n_unlinked = sum(!linked))
}

# Analysed facilities must all have a stratum (fail closed, fixed message). The
# census table used for display maps a missing value to "(missing)" instead.
stratum_of <- function(nk, strata, display = FALSE) {
  vars <- switch(strata, region_x_type = c("region", "type"), region = "region", type = "type",
                 typegroup = "typegroup", region_x_typegroup = c("region", "typegroup"), none = character(0))
  if (length(vars) && !all(vars %in% names(nk))) fail("stratum variable not loaded for strata = ", strata)
  if (display) {
    for (x in vars) nk[[x]][is.na(nk[[x]])] <- "(missing)"
  } else if (length(vars) && anyNA(nk[vars])) {
    fail("stratum variable missing for some analysed facilities; the frozen gate does not define this (nothing written)")
  }
  switch(strata, region_x_type = paste(nk$region, nk$type, sep = " x "),
         region = nk$region, type = nk$type, typegroup = nk$typegroup,
         region_x_typegroup = paste(nk$region, nk$typegroup, sep = " x "), none = rep("all", nrow(nk)))
}

compute_metrics <- function(res, gate) {
  nk <- res$nk$nk; n <- length(nk)
  if (n == 0) fail("no facilities left after applying the gate's facility rules")
  c_take <- gate$thresholds$stage2_c
  st <- table(stratum_of(res$nk, gate$strata))
  if (length(st) < 1 || sum(st) != n) fail("stratum table does not cover every analysed facility (nothing written)")
  list(
    scalars = c(n_facilities = n, n_records = sum(nk), mean_nk = mean(nk), median_nk = median(nk),
                max_nk = max(nk), share_nk_0 = mean(nk == 0), share_nk_le_1 = mean(nk <= 1),
                share_nk_le_2 = mean(nk <= 2),
                share_nk_gt_c = if (identical(c_take, "off")) NA_real_ else mean(nk > c_take),
                n_strata = length(st), min_facilities_per_stratum = min(st),
                n_unlinked_client_records = res$n_unlinked),
    strata_census = table(stratum_of(res$census, gate$strata, display = TRUE)),
    bins = bin_counts(nk, gate$reporting_bins))
}

bin_counts <- function(nk, edges) {
  lab <- c(ifelse(diff(edges) == 1, as.character(head(edges, -1)),
                  paste0(head(edges, -1), "-", tail(edges, -1) - 1)),
           paste0(tail(edges, 1), "+"))
  idx <- findInterval(nk, edges)
  setNames(as.integer(tabulate(idx, nbins = length(edges))), lab)
}

evaluate_gate <- function(m, gate) {
  th <- gate$thresholds; s <- m$scalars
  rules <- list(
    min_mean_nk                = list(val = s[["mean_nk"]],                    ok = function(v, t) v >= t),
    max_share_nk_0             = list(val = s[["share_nk_0"]],                 ok = function(v, t) v <= t),
    max_share_nk_le_1          = list(val = s[["share_nk_le_1"]],              ok = function(v, t) v <= t),
    max_share_nk_le_2          = list(val = s[["share_nk_le_2"]],              ok = function(v, t) v <= t),
    min_facilities_per_stratum = list(val = s[["min_facilities_per_stratum"]], ok = function(v, t) v >= t),
    min_share_nk_gt_c          = list(val = s[["share_nk_gt_c"]],              ok = function(v, t) v >= t))
  rows <- lapply(names(rules), function(nm) {
    t <- th[[nm]]
    if (identical(t, "off")) return(data.frame(criterion = nm, threshold = "off", value = NA_real_, result = "NOT EVALUATED"))
    r <- rules[[nm]]
    data.frame(criterion = nm, threshold = as.character(t), value = r$val,
               result = if (isTRUE(r$ok(r$val, t))) "PASS" else "FAIL")
  })
  tab <- do.call(rbind, rows)
  list(table = tab, overall = if (all(tab$result %in% c("PASS", "NOT EVALUATED"))) "PASS" else "FAIL")
}

# --- data loading (count mode only, after preflight) ------------------------------
load_survey <- function(data_root, survey, gate, reader) {
  v <- VARS[[gate$count_source]]
  fc_path <- find_dta(data_root, survey, gate$count_source, "FC")
  if (is.na(fc_path)) fail("facility file not found (expected exactly one .dta for ", survey, "/", gate$count_source, "/FC)")
  cross <- v$fc_counts[gate$services]
  use_cross <- all(!is.na(cross))
  rv <- var_for(v$region, survey); tv <- var_for(v$type, survey)
  fc <- reader(fc_path, cols = c(v$fac_id, v$fac_done, rv, tv, if (use_cross) cross))
  need <- c(v$fac_id, v$fac_done, rv, tv, if (use_cross) cross)
  if (!all(need %in% names(fc))) fail("facility file lacks configured variables; run --structure")
  fac <- data.frame(id = as.vector(fc[[v$fac_id]]), done = as.vector(fc[[v$fac_done]]),
                    region = as.character(haven::as_factor(fc[[rv]], levels = "both")),
                    type = as.character(haven::as_factor(fc[[tv]], levels = "both")), stringsAsFactors = FALSE)
  if (grepl("typegroup", gate$strata)) {
    if (gate$count_source != "recode") fail("typegroup strata are defined for recode (v007) codes only")
    fac$typegroup <- type_group_of(fc[[tv]], survey)
  }
  cli <- list()
  for (s in gate$services) {
    p <- find_dta(data_root, survey, gate$count_source, s)
    if (is.na(p)) fail("client file not found for ", s)
    d <- reader(p, cols = c(v$cli_id, v$cli_obs[[s]], v$cli_exit[[s]]))
    if (!all(c(v$cli_id, v$cli_obs[[s]], v$cli_exit[[s]]) %in% names(d))) fail(s, " file lacks configured variables; run --structure")
    cli[[s]] <- data.frame(id = as.vector(d[[v$cli_id]]), obs = as.vector(d[[v$cli_obs[[s]]]]),
                           exit = as.vector(d[[v$cli_exit[[s]]]]))
  }
  xc <- NULL
  if (use_cross) {
    done <- !is.na(fac$done) & fac$done == DONE_CODE
    xc <- data.frame(id = fac$id[done],
                     fc_total = Reduce(`+`, lapply(cross, function(x) as.numeric(fc[[x]][done]))))  # works for one service
  }
  list(fac = fac, cli = cli, cross = xc)
}

# The same four conditions --structure reports. Fixed messages: no counts, IDs or positions.
check_linkage <- function(fac, cli) {
  if (anyNA(fac$id)) fail("facility ID missing in the facility file; run --structure (nothing written)")
  if (anyDuplicated(fac$id)) fail("facility ID not unique in the facility file; run --structure (nothing written)")
  for (s in names(cli)) {
    ci <- cli[[s]]$id
    if (is.numeric(ci) != is.numeric(fac$id)) fail(s, " facility ID type differs from the facility file; run --structure (nothing written)")
    if (anyNA(ci)) fail(s, " facility ID missing on client records; run --structure (nothing written)")
    if (!all(ci %in% fac$id)) fail(s, " client facility IDs absent from the facility file; run --structure (nothing written)")
  }
}

# --- count mode --------------------------------------------------------------------
# `reader` is injectable so the self-test can prove no data file is read before
# the gate and preflight pass.
count_mode <- function(survey, gate_path, repo, rel_paths, data_root, manifest_path, out_dir,
                       reader = read_dta_lower, now_utc = Sys.time()) {
  gate <- validate_gate(read_gate_yaml(gate_path))
  if (!survey %in% gate$surveys) fail("survey ", survey, " is not listed in the frozen gate")
  prov <- preflight_code(repo, rel_paths)
  if (!file.exists(manifest_path)) fail("raw/MANIFEST.txt missing: run 00-unpack.R first")
  prov <- c(prov, manifest_md5 = unname(tools::md5sum(manifest_path)),
            run_utc = format(now_utc, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  # ---- nothing above this line opens a data file ----
  dat <- load_survey(data_root, survey, gate, reader)
  check_linkage(dat$fac, dat$cli)
  res <- compute_nk(dat$fac, dat$cli, gate$services, gate$client_record, gate$zero_record_facilities)
  m   <- compute_metrics(res, gate)
  if (!is.null(dat$cross)) {
    j <- merge(res$census[, c("id", "nk")], dat$cross, by = "id")
    m$scalars <- c(m$scalars, crosscheck_share_equal = mean(!is.na(j$fc_total) & j$nk == j$fc_total))
  }
  ev <- evaluate_gate(m, gate)
  write_outputs(out_dir, survey, prov, gate, m, ev)
  invisible(list(metrics = m, eval = ev, provenance = prov))
}

write_outputs <- function(out_dir, survey, prov, gate, m, ev) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- gsub("[-:]", "", prov[["run_utc"]])
  stem  <- sprintf("count-audit_%s_%s_%s", survey, stamp, substr(prov[["head"]], 1, 7))
  csv <- file.path(out_dir, paste0(stem, ".csv")); md <- file.path(out_dir, paste0(stem, ".md"))
  if (file.exists(csv) || file.exists(md)) fail("refusing to overwrite an earlier run: ", stem)
  metrics <- rbind(
    data.frame(section = "metric", name = names(m$scalars), value = unname(m$scalars)),
    data.frame(section = "nk_bin", name = names(m$bins), value = unname(m$bins)),
    data.frame(section = "gate", name = ev$table$criterion, value = ev$table$value))
  write.csv(metrics, csv, row.names = FALSE)
  writeLines(c(
    sprintf("# SPA count audit -- %s -- %s", survey, ev$overall), "",
    "## Provenance", "", paste0("- `", names(prov), "`: `", prov, "`"), "",
    sprintf("Gate: frozen by %s on %s; count_source %s; client_record %s; services %s; zero-record facilities %s; strata %s.",
            gate$frozen_by, gate$frozen_on, gate$count_source, gate$client_record,
            paste(gate$services, collapse = "+"), gate$zero_record_facilities, gate$strata), "",
    "## Gate evaluation", "", "| criterion | threshold | value | result |", "|---|---|---|---|",
    sprintf("| %s | %s | %s | %s |", ev$table$criterion, ev$table$threshold,
            ifelse(is.na(ev$table$value), "", formatC(ev$table$value, digits = 4, format = "g")), ev$table$result), "",
    sprintf("**Overall: %s**", ev$overall), "",
    "## n_k distribution (facilities per bin)", "", "| n_k | facilities |", "|---|---|",
    sprintf("| %s | %d |", names(m$bins), m$bins), "",
    "## Facilities per stratum (census of completed facilities, before the zero-record rule)", "",
    "| stratum | facilities |", "|---|---|",
    sprintf("| %s | %d |", names(m$strata_census), as.integer(m$strata_census))), md)
  runs <- file.path(out_dir, "RUNS.md")
  if (!file.exists(runs)) writeLines(c("# Count-audit runs (append-only)", "",
                                       "| run (UTC) | survey | HEAD | gate blob | script blob | result | files |",
                                       "|---|---|---|---|---|---|---|"), runs)
  blob <- function(k) substr(prov[grep(paste0("^blob:.*", k), names(prov))][1], 1, 12)
  if (anyNA(c(blob("GATE.md"), blob("count-audit.R")))) fail("provenance incomplete: gate or script blob missing")
  cat(sprintf("| %s | %s | %s | %s | %s | %s | %s |\n", prov[["run_utc"]], survey, substr(prov[["head"]], 1, 12),
              blob("GATE.md"), blob("count-audit.R"), ev$overall, stem), file = runs, append = TRUE)
  message(sprintf("%s: %s  (wrote results/%s.{csv,md}; appended RUNS.md)", survey, ev$overall, stem))
  print(ev$table, row.names = FALSE)
}

# --- structure mode ------------------------------------------------------------------
# Returns output lines. Only presence/class/label, value-label checks and
# TRUE/FALSE linkage results; no counts, IDs or importer text.
structure_lines <- function(data_root, survey, reader = read_dta_lower) {
  out <- character(0); say <- function(...) out <<- c(out, paste0(...))
  lab_of <- function(x) {                     # exact = TRUE: "label" must not partially match "labels"
    l <- attr(x, "label", exact = TRUE)
    if (!is.character(l) || length(l) != 1 || is.na(l)) "(no label)" else substr(gsub("[[:cntrl:]]", " ", l), 1, 80)
  }
  has_code <- function(x, code) { vl <- attr(x, "labels"); !is.null(vl) && code %in% unname(vl) }
  for (structure in c("recode", "raw")) {
    v <- VARS[[structure]]
    say("== ", survey, " / ", structure)
    files <- vapply(c("FC", SERVICES), function(mod) find_dta(data_root, survey, structure, mod), "")
    for (mod in names(files)) if (is.na(files[[mod]])) say("[file] ", mod, ": NOT FOUND (need exactly one .dta)")
    fc_ids <- NULL
    if (!is.na(files[["FC"]])) {
      meta <- tryCatch(reader(files[["FC"]], n_max = 0), error = function(e) NULL)
      if (is.null(meta)) say("[file] FC: read error (details suppressed)") else {
        for (nm in unique(c(v$fac_id, v$fac_done, var_for(v$region, survey), var_for(v$type, survey),
                            v$fc_count_candidates, stats::na.omit(v$fc_counts)))) {
          if (nm %in% names(meta)) say("[schema] FC ", nm, ": present, class=", class(meta[[nm]])[1], ", label=", lab_of(meta[[nm]]))
          else say("[schema] FC ", nm, ": MISSING")
        }
        if (v$fac_done %in% names(meta)) say("[codes] FC ", v$fac_done, " has value label for code ", DONE_CODE, ": ", has_code(meta[[v$fac_done]], DONE_CODE))
        idd <- tryCatch(reader(files[["FC"]], cols = v$fac_id), error = function(e) NULL)
        if (!is.null(idd) && v$fac_id %in% names(idd)) {
          fc_ids <- as.vector(idd[[v$fac_id]])
          say("[link] FC ", v$fac_id, " never missing: ", !anyNA(fc_ids))
          say("[link] FC ", v$fac_id, " unique: ", !anyDuplicated(fc_ids))
        }
      }
    }
    for (s in SERVICES) {
      if (is.na(files[[s]])) next
      meta <- tryCatch(reader(files[[s]], n_max = 0), error = function(e) NULL)
      if (is.null(meta)) { say("[file] ", s, ": read error (details suppressed)"); next }
      for (nm in unique(c(v$cli_id, v$cli_obs[[s]], v$cli_exit[[s]]))) {
        if (nm %in% names(meta)) say("[schema] ", s, " ", nm, ": present, class=", class(meta[[nm]])[1], ", label=", lab_of(meta[[nm]]))
        else say("[schema] ", s, " ", nm, ": MISSING")
      }
      for (nm in unique(c(v$cli_obs[[s]], v$cli_exit[[s]])))
        if (nm %in% names(meta)) say("[codes] ", s, " ", nm, " has value label for code ", DONE_CODE, ": ", has_code(meta[[nm]], DONE_CODE))
      if (!is.null(fc_ids) && v$cli_id %in% names(meta)) {
        cd <- tryCatch(reader(files[[s]], cols = v$cli_id), error = function(e) NULL)
        if (!is.null(cd)) {
          ci <- as.vector(cd[[v$cli_id]])
          say("[link] ", s, " ", v$cli_id, " never missing: ", !anyNA(ci))
          say("[link] ", s, " ", v$cli_id, " same type as FC id: ", identical(typeof(ci), typeof(fc_ids)))
          say("[link] ", s, " every facility id found in FC: ", all(ci[!is.na(ci)] %in% fc_ids))
        }
      }
    }
  }
  out
}

# Belt and braces: structure output may contain only the fixed line shapes above.
STRUCTURE_LINE_OK <- paste0(
  "^(== .+ / (recode|raw)",
  "|\\[file\\] [A-Z]{2}: (NOT FOUND \\(need exactly one \\.dta\\)|read error \\(details suppressed\\))",
  "|\\[schema\\] [A-Z]{2} [a-z0-9_]+: (MISSING|present, class=[A-Za-z_.]+, label=.*)",
  "|\\[codes\\] [A-Z]{2} [a-z0-9_]+ has value label for code 1: (TRUE|FALSE)",
  "|\\[link\\] [A-Z]{2} .+: (TRUE|FALSE))$")
check_structure_lines <- function(lines) {
  bad <- !grepl(STRUCTURE_LINE_OK, lines)
  if (any(bad)) fail("structure output contained an unexpected line shape; suppressed")
  lines
}

# --- self-test -----------------------------------------------------------------------
selftest <- function() {
  n_ok <- 0
  check <- function(cond, what) { if (!isTRUE(cond)) fail("SELFTEST FAILED: ", what); n_ok <<- n_ok + 1 }
  expect_error <- function(expr, pattern, what) {
    msg <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
    check(!is.null(msg) && grepl(pattern, msg), paste0(what, " (got: ", msg %||% "no error", ")"))
  }
  tmp <- tempfile("count-audit-selftest-"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  # Synthetic census: 12 completed facilities (ids 101-112) + 2 not completed (113-114).
  # Intended n_k (obs+exit, AN+FP+SC): 101:0 102:1 103:1 104:2 105:3 106:4 107:5 108:5 109:6 110:8 111:10 112:12
  want <- c(0, 1, 1, 2, 3, 4, 5, 5, 6, 8, 10, 12); ids <- 101:112
  fac <- data.frame(id = c(ids, 113L, 114L), done = c(rep(1, 12), 2, 3),
                    region = rep(c("N", "C", "S"), length.out = 14), type = rep(c("hosp", "hc"), each = 7))
  per <- rep(ids, want)                                  # one record per client, round-robin over services
  svc <- rep_len(SERVICES, length(per))
  mk  <- function(s) {
    d <- data.frame(id = per[svc == s], obs = 1, exit = 1)
    rbind(d,
      data.frame(id = 105L, obs = 2, exit = 1),          # observation incomplete -> excluded
      data.frame(id = 106L, obs = 1, exit = 2),          # exit refused -> excluded unless obs_only
      data.frame(id = 999L, obs = 1, exit = 1),          # facility absent from FC -> unlinked
      data.frame(id = 113L, obs = 1, exit = 1))          # facility not completed -> unlinked
  }
  cli <- setNames(lapply(SERVICES, mk), SERVICES)

  r <- compute_nk(fac, cli, SERVICES, "obs_and_exit", "keep")
  check(identical(r$nk$nk, as.integer(want)), "n_k by facility (keep, obs_and_exit)")
  check(r$n_unlinked == 6, "unlinked client records counted (2 per service)")
  r_drop <- compute_nk(fac, cli, SERVICES, "obs_and_exit", "drop")
  check(nrow(r_drop$nk) == 11 && all(r_drop$nk$nk > 0), "zero-record facility dropped")
  r_obs <- compute_nk(fac, cli, SERVICES, "obs_only", "keep")
  check(r_obs$nk$nk[r_obs$nk$id == 106] == 4 + 3, "obs_only counts exit refusals")
  r_an <- compute_nk(fac, cli, "AN", "obs_and_exit", "keep")
  check(sum(r_an$nk$nk) == sum(svc == "AN"), "single-service subset")

  gate_ok <- list(status = "FROZEN", frozen_by = "Test", frozen_on = "2026-01-01", surveys = list("MW_SPA_2013-14"),
                  count_source = "raw", client_record = "obs_and_exit", services = as.list(SERVICES),
                  zero_record_facilities = "keep", strata = "region",
                  thresholds = list(min_mean_nk = 4, max_share_nk_0 = 0.1, max_share_nk_le_1 = 0.25,
                                    max_share_nk_le_2 = "off", min_facilities_per_stratum = 4,
                                    stage2_c = 3, min_share_nk_gt_c = 0.5),
                  reporting_bins = list(0, 1, 2, 3, 5, 10))
  g <- validate_gate(gate_ok)
  m <- compute_metrics(r, g)
  check(abs(m$scalars[["mean_nk"]] - mean(want)) < 1e-12, "mean n_k")
  check(abs(m$scalars[["share_nk_0"]] - 1/12) < 1e-12 && abs(m$scalars[["share_nk_le_1"]] - 3/12) < 1e-12, "left-tail shares")
  check(abs(m$scalars[["share_nk_gt_c"]] - mean(want > 3)) < 1e-12, "share n_k > c")
  check(m$scalars[["min_facilities_per_stratum"]] == 4 && m$scalars[["n_strata"]] == 3, "strata (region)")
  check(identical(unname(m$bins), c(1L, 2L, 1L, 2L, 4L, 2L)) &&
        identical(names(m$bins), c("0", "1", "2", "3-4", "5-9", "10+")), "binning")
  ev <- evaluate_gate(m, g)
  # mean 4.75 >= 4 PASS; share0 .083 <= .1 PASS; share<=1 .25 <= .25 PASS; <=2 off; min/stratum 4 >= 4 PASS; >c 7/12 >= .5 PASS
  check(ev$overall == "PASS" && sum(ev$table$result == "NOT EVALUATED") == 1, "gate evaluation PASS path")
  g2 <- g; g2$thresholds$min_mean_nk <- 5
  check(evaluate_gate(m, g2)$overall == "FAIL", "gate evaluation FAIL path")

  # Gate validation rejects drafts, nulls, placeholders and all-off gates.
  expect_error(validate_gate(modifyList(gate_ok, list(status = "DRAFT"))), "status must be FROZEN", "DRAFT rejected")
  bad <- gate_ok; bad$thresholds$min_mean_nk <- NULL
  expect_error(validate_gate(bad), "min_mean_nk unset", "missing threshold rejected")
  expect_error(validate_gate(modifyList(gate_ok, list(frozen_by = "TBD"))), "frozen_by unset", "placeholder rejected")
  alloff <- gate_ok; alloff$thresholds <- lapply(alloff$thresholds, function(x) "off")
  expect_error(validate_gate(alloff), "every threshold is off", "all-off gate rejected")
  nc <- gate_ok; nc$thresholds$stage2_c <- "off"
  expect_error(validate_gate(nc), "needs stage2_c", "c-share without c rejected")
  with_field <- function(nm, val) { x <- gate_ok; x[nm] <- list(val); x }   # replace, not merge
  expect_error(validate_gate(with_field("services", list(NULL))), "services must be", "services: [null] rejected")
  expect_error(validate_gate(with_field("services", list(list()))), "services must be", "services: [{}] rejected")
  expect_error(validate_gate(with_field("surveys", list())), "surveys must be", "empty surveys rejected")
  expect_error(validate_gate(with_field("reporting_bins", list(0, 1.5, 3))), "integer lower edges",
               "non-integer bins rejected")
  check(identical(validate_gate(with_field("services", list("AN")))$services, "AN"), "single-service gate accepted")

  # Strata: an analysed facility without a stratum fails closed; the census display
  # table keeps it as "(missing)"; under drop, the display table is the pre-drop census.
  fac_na <- fac; fac_na$region[fac_na$id == 102] <- NA
  r_na <- compute_nk(fac_na, cli, SERVICES, "obs_and_exit", "keep")
  expect_error(compute_metrics(r_na, g), "stratum variable missing", "missing stratum fails closed")
  check("(missing)" %in% names(table(stratum_of(r_na$census, "region", display = TRUE))), "display maps missing stratum")
  check(identical(type_group_of(c(1, 4, 5, 6, 7, 9, 10, NA), "MW_SPA_2013-14"),
                  c("hospital", "hospital", "health centre", "health centre", "lower-level", "lower-level", NA, NA)),
        "Malawi type groups (unmapped code -> NA)")
  check(identical(type_group_of(c(4, 6, 7, 8), "HT_SPA_2013"), c("hospital", "health centre", "lower-level", NA)),
        "Haiti type groups")
  fac_tg <- fac; fac_tg$typegroup <- type_group_of(ifelse(fac$type == "hosp", 2, 8), "MW_SPA_2013-14")
  r_tg <- compute_nk(fac_tg, cli, SERVICES, "obs_and_exit", "keep")
  check(identical(as.vector(table(stratum_of(r_tg$nk, "typegroup"))), c(7L, 5L)), "typegroup strata carried through compute_nk")
  fac_tg$typegroup[fac_tg$id == 103] <- NA
  expect_error(stratum_of(compute_nk(fac_tg, cli, SERVICES, "obs_and_exit", "keep")$nk, "typegroup"),
               "stratum variable missing", "unmapped type code fails closed")
  g_drop <- g; g_drop$zero_record_facilities <- "drop"
  m_drop <- compute_metrics(r_drop, g_drop)
  check(m_drop$scalars[["n_facilities"]] == 11 && sum(m_drop$strata_census) == 12,
        "strata table written from the pre-drop census (no n_k = 0 disclosure)")

  # The real DRAFT gate in this repo must be rejected too.
  real_gate <- GATE_FILE
  if (file.exists(real_gate)) {
    rg <- read_gate_yaml(real_gate)
    if (!identical(rg$status, "FROZEN"))
      expect_error(validate_gate(rg), "GATE NOT USABLE", "repository GATE.md (not frozen) rejected")
  }

  # Code preflight and the no-data-before-preflight guarantee. First a packaged copy (no git,
  # a CODE-MANIFEST.tsv), then -- where git is installed -- a throwaway git repo.
  gate_yaml <- c("---", yaml::as.yaml(gate_ok), "---", "", "body")
  draft <- sub("status: FROZEN", "status: DRAFT", gate_yaml, fixed = TRUE)
  pf <- c("GATE.md", "count-audit.R")
  touched <- 0
  spy <- function(...) { touched <<- touched + 1; fail("spy reader should not be reached") }
  run_in <- function(repo) count_mode("MW_SPA_2013-14", file.path(repo, "GATE.md"), repo, pf,
                                      data_root = tmp, manifest_path = file.path(tmp, "none"), out_dir = tmp, reader = spy)
  pkg <- file.path(tmp, "pkg"); dir.create(pkg)
  writeLines(gate_yaml, file.path(pkg, "GATE.md")); writeLines("# script", file.path(pkg, "count-audit.R"))
  freeze <- function(files = pf) writeLines(c(paste0("# head\t", strrep("a", 40)), "path\tmd5\tblob\tlast_commit",
    sprintf("%s\t%s\t%s\t%s", files, unname(tools::md5sum(file.path(pkg, files))), strrep("b", 40), strrep("c", 40))),
    file.path(pkg, CODE_MANIFEST))
  expect_error(run_in(pkg), "neither a git checkout nor", "no git and no manifest fails preflight")
  freeze("GATE.md")
  expect_error(run_in(pkg), "count-audit.R: not in CODE-MANIFEST", "file missing from the manifest fails preflight")
  freeze()
  prov <- preflight_code(pkg, pf)
  check(prov[["head"]] == strrep("a", 40) && prov[["blob:GATE.md"]] == strrep("b", 40) &&
        prov[["last_commit:count-audit.R"]] == strrep("c", 40), "packaged preflight returns the manifest's hashes")
  cat("\n# edit\n", file = file.path(pkg, "GATE.md"), append = TRUE)
  expect_error(run_in(pkg), "GATE.md: differs from CODE-MANIFEST", "modified gate fails packaged preflight")
  writeLines(gate_yaml, file.path(pkg, "GATE.md"))
  expect_error(run_in(pkg), "MANIFEST.txt missing", "unchanged packaged copy passes preflight, stops at manifest")
  writeLines(draft, file.path(pkg, "GATE.md")); freeze()
  expect_error(run_in(pkg), "status must be FROZEN", "frozen-in-manifest DRAFT gate never unlocks count mode")
  writeLines(gate_yaml, file.path(pkg, "GATE.md")); freeze()

  if (nzchar(Sys.which("git"))) {
    repo <- file.path(tmp, "repo"); dir.create(repo)
    g_ <- function(...) git_run(repo, c("-c", "user.name=selftest", "-c", "user.email=selftest@example.invalid", ...))
    g_("init", "-q"); g_("config", "core.autocrlf", "false")
    writeLines(gate_yaml, file.path(repo, "GATE.md")); writeLines("# script", file.path(repo, "count-audit.R"))
    file.copy(file.path(pkg, CODE_MANIFEST), repo)           # a stray manifest never overrides git
    run <- function() run_in(repo)
    expect_error(run(), "not tracked", "untracked files fail preflight")
    g_("add", "GATE.md", "count-audit.R"); g_("commit", "-q", "-m", "freeze")
    prov <- preflight_git(repo, pf)
    check(nchar(prov[["head"]]) == 40 && nchar(prov[["blob:GATE.md"]]) == 40, "preflight returns hashes")
    cat("\n# edit\n", file = file.path(repo, "GATE.md"), append = TRUE)
    expect_error(run(), "working tree differs", "modified gate fails preflight")
    g_("add", "GATE.md")
    expect_error(run(), "index differs from HEAD", "staged change fails preflight")
    g_("reset", "-q", "--hard")
    expect_error(run(), "MANIFEST.txt missing", "clean repo passes preflight, stops at manifest")
    writeLines(draft, file.path(repo, "GATE.md")); g_("commit", "-q", "-am", "draft")
    expect_error(run(), "status must be FROZEN", "committed DRAFT gate never unlocks count mode")
  } else message("selftest: git is not installed -- git-checkout preflight checks skipped ",
                 "(a packaged copy uses the manifest checks above)")
  check(touched == 0, "no data read before gate + preflight passed")

  # End-to-end on synthetic .dta files, plus non-overwrite, in the packaged copy (no git needed).
  repo <- pkg
  root <- file.path(tmp, "data"); sv <- "MW_SPA_2013-14"; v <- VARS$raw
  wdta <- function(mod, df) {
    d <- file.path(root, sv, "raw", FILES[[sv]]$raw[[mod]]); dir.create(d, recursive = TRUE)
    haven::write_dta(df, file.path(d, "X.dta"))
  }
  lab <- function(x) haven::labelled(x, c(completed = 1, other = 2))
  # Facility IDs carry per-facility value labels (names) and no variable label: the
  # structure output must still never list them.
  fac_id_lab <- haven::labelled(fac$id, setNames(fac$id, paste0("Facname", fac$id)))
  wdta("FC", setNames(data.frame(fac_id_lab, lab(fac$done), fac$region, fac$type), c(v$fac_id, v$fac_done, var_for(v$region, sv), var_for(v$type, sv))))
  wcli <- function(cl) for (s in SERVICES) {
    d <- file.path(root, sv, "raw", FILES[[sv]]$raw[[s]]); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    haven::write_dta(setNames(data.frame(cl[[s]]$id, lab(cl[[s]]$obs), lab(cl[[s]]$exit)),
                              c(v$cli_id, v$cli_obs[[s]], v$cli_exit[[s]])), file.path(d, "X.dta"))
  }
  wcli(cli)                                              # includes client records of facility 999 (absent from FC)
  man <- file.path(tmp, "MANIFEST.txt"); writeLines("synthetic", man)
  outd <- file.path(tmp, "results"); t0 <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")

  # Structure mode on these files: fixed shapes only, no counts, no IDs or value labels.
  sl <- structure_lines(root, sv)
  check(length(check_structure_lines(sl)) > 0, "structure output passes the line-shape filter")
  numeric_hit <- !grepl("^== ", sl) & grepl("\\b(14|12|6)\\b", sub("label=.*$", "", sl))
  check(!any(numeric_hit), paste0("structure output has no row counts; offending: ",
                                  paste(sl[numeric_hit], collapse = " | ")))
  check(!any(grepl("Facname|\\b(10[1-9]|11[0-4]|999)\\b", sl)), "structure output lists no facility IDs or value labels")
  sch <- sub(":.*$", "", grep("^\\[schema\\]", sl, value = TRUE))
  check(!anyDuplicated(sch), "exactly one [schema] line per configured variable")
  check(any(sl == paste0("[schema] FC ", v$fac_done, ": present, class=haven_labelled, label=(no label)")),
        "value labels are not echoed as the variable label")
  check(any(grepl("^\\[link\\] FP every facility id found in FC: FALSE$", sl)), "structure linkage is boolean only")
  expect_error(check_structure_lines(c(sl, "rows: 123")), "unexpected line shape", "line-shape filter rejects counts")

  # Broken linkage fails closed after preflight, before anything is computed or written.
  expect_error(suppressMessages(count_mode(sv, file.path(repo, "GATE.md"), repo, c("GATE.md", "count-audit.R"),
                                           root, man, outd, now_utc = t0)),
               "absent from the facility file", "unlinked client IDs fail closed")
  check(!dir.exists(outd) || !length(list.files(outd)), "linkage failure wrote nothing")

  cli_ok <- lapply(cli, function(d) d[d$id != 999L, ])
  wcli(cli_ok)
  invisible(capture.output(res <- suppressMessages(count_mode(sv, file.path(repo, "GATE.md"), repo,
    c("GATE.md", "count-audit.R"), root, man, outd, now_utc = t0))))
  check(res$eval$overall == "PASS" && abs(res$metrics$scalars[["mean_nk"]] - mean(want)) < 1e-12, "end-to-end on .dta")
  expect_error(suppressMessages(count_mode(sv, file.path(repo, "GATE.md"), repo, c("GATE.md", "count-audit.R"),
                                           root, man, outd, now_utc = t0)), "refusing to overwrite", "non-overwrite")
  runs <- readLines(file.path(outd, "RUNS.md"))
  check(sum(grepl("^\\| 2026-01-01", runs)) == 1, "RUNS.md appended once (failed rerun wrote nothing)")
  outtxt <- paste(c(readLines(list.files(outd, "\\.md$", full.names = TRUE)[1]),
                    readLines(list.files(outd, "\\.csv$", full.names = TRUE)[1])), collapse = "\n")
  check(!grepl("\\b(10[1-9]|11[0-4]|999)\\b|Facname", outtxt), "no facility identifier in outputs")

  # Recode path with a single-service gate and the facility-file cross-check.
  vr <- VARS$recode
  rdta <- function(mod, df) {
    d <- file.path(root, sv, "recode", FILES[[sv]]$recode[[mod]]); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    haven::write_dta(df, file.path(d, "X.dta"))
  }
  an_true <- as.integer(table(factor(cli_ok$AN$id[cli_ok$AN$obs == 1 & cli_ok$AN$exit == 1], levels = fac$id)))
  an_fc <- an_true; an_fc[fac$id == 110] <- an_fc[fac$id == 110] + 1L    # one deliberate disagreement
  rdta("FC", setNames(data.frame(fac$id, lab(fac$done), fac$region, fac$type, an_fc),
                      c(vr$fac_id, vr$fac_done, vr$region, vr$type, vr$fc_counts[["AN"]])))
  rdta("AN", setNames(data.frame(cli_ok$AN$id, lab(cli_ok$AN$obs), lab(cli_ok$AN$exit)),
                      c(vr$cli_id, vr$cli_obs[["AN"]], vr$cli_exit[["AN"]])))
  g_rec <- gate_ok; g_rec$count_source <- "recode"; g_rec$services <- list("AN")
  g_rec$thresholds$min_mean_nk <- 0; g_rec$thresholds$min_share_nk_gt_c <- "off"; g_rec$thresholds$max_share_nk_0 <- "off"
  g_rec$thresholds$max_share_nk_le_1 <- "off"
  writeLines(c("---", yaml::as.yaml(g_rec), "---"), file.path(repo, "GATE.md"))
  freeze()                                               # the recode gate is now the frozen one
  invisible(capture.output(rr <- suppressMessages(count_mode(sv, file.path(repo, "GATE.md"), repo,
    c("GATE.md", "count-audit.R"), root, man, outd, now_utc = t0 + 60))))
  check(rr$metrics$scalars[["n_records"]] == sum(an_true[fac$done == 1]), "recode single-service count (completed facilities only)")
  check(abs(rr$metrics$scalars[["crosscheck_share_equal"]] - 11/12) < 1e-12, "cross-check detects the one disagreement")

  message(sprintf("selftest: %d checks passed (synthetic data only; no real file opened)", n_ok))
}

# --- entry point ---------------------------------------------------------------------
DHS <- local({
  for (cand in c(file.path("src", "application"), ".", "..")) {
    if (file.exists(file.path(cand, "R", "01-count-audit.R"))) return(normalizePath(cand, winslash = "/"))
  }
  fail("cannot locate src/application/ (run from the repository root)")
})
REPO <- normalizePath(file.path(DHS, "..", ".."), winslash = "/")
GATE_FILE <- file.path(REPO, "config", "application", "gate.yaml")
DATA_ROOT <- function() Sys.getenv("SPA_RAW_DIR", file.path(REPO, "data", "spa"))
OUT_ROOT  <- function() Sys.getenv("SPA_RESULTS_DIR", file.path(REPO, "outputs", "application"))
rel  <- function(p) {
  a <- normalizePath(p, winslash = "/"); pre <- paste0(REPO, "/")
  if (!startsWith(a, pre)) fail("path outside the repository: ", p)
  substring(a, nchar(pre) + 1)
}

main <- function(args) {
  modes <- intersect(args, c("--selftest", "--structure", "--count"))
  usage <- "usage: Rscript src/application/R/01-count-audit.R (--selftest | --structure --survey <S> | --count --survey <S>)"
  if (length(modes) != 1) fail("exactly one mode flag is required; there is no default.\n", usage)
  survey <- if ("--survey" %in% args) args[which(args == "--survey") + 1] else NA
  if (modes != "--selftest" && (is.na(survey) || !survey %in% SURVEYS))
    fail("--survey must be one of: ", paste(SURVEYS, collapse = ", "), "\n", usage)
  switch(modes,
    "--selftest"  = selftest(),
    "--structure" = writeLines(check_structure_lines(structure_lines(DATA_ROOT(), survey))),
    "--count"     = {
      ca     <- commandArgs(FALSE)
      ca     <- ca[seq_len(match("--args", ca, nomatch = length(ca) + 1L) - 1L)]   # R's own options only
      self   <- sub("^--file=", "", grep("^--file=", ca, value = TRUE))
      script <- normalizePath(file.path(DHS, "R", "01-count-audit.R"), winslash = "/")
      if (length(self) != 1 || !identical(normalizePath(self, winslash = "/", mustWork = FALSE), script))
        fail("count mode must be run directly: Rscript src/application/R/01-count-audit.R --count --survey <S>")
      count_mode(survey, GATE_FILE, REPO, c(rel(GATE_FILE), rel(script)),
                 DATA_ROOT(), file.path(DATA_ROOT(), "MANIFEST.txt"), OUT_ROOT())
    })
}

if (sys.nframe() == 0) main(commandArgs(trailingOnly = TRUE))
