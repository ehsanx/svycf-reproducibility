# Run modes, identities and the frozen-input check.
#
#   dev   SYNTHETIC root; populations 1..B, reps 1..R; DEV output. Kept for unit tests and code
#         checks only: no recorded run uses it.
#   dry   a dry run: the production root; production population IDs, shared
#         with production; rep IDs on the disjoint range 1000001, 1000002, ..., so no dry-run
#         replicate is a production one. Census-only cells (R = 0) have
#         no rep axis and use the reserved population IDs 1000001-1000005 instead.
#   prod  production: the production root, production population IDs, rep IDs 1..R.
#
# dry and prod refuse to run unless every frozen input listed in config/simulation/frozen-design.json has its
# recorded SHA-256, and record the simulation design specification's SHA-256 from that file (the
# specification document itself is not distributed).
# Requires R/rng_streams.R (SIM_PROD_ROOT) sourced first.

SIM_SPEC_VERSION <- "v9.8"
SIM_FROZEN_DESIGN <- file.path("config", "simulation", "frozen-design.json")
DRY_ID_BASE <- 1000000L                 # dry-run rep IDs (and census-only dry populations) start above it
SIM_DEV_ROOT <- "SYNTHETIC-DEV-E2E"

# Production population IDs per master: fitted 1..B (MS, FS 1..400; FPC 1..160) and the E1
# oracle anchors 100001..100012 on MS and FS.
PROD_POPS <- list(MS = c(1:400, 100001:100012), FS = c(1:400, 100001:100012), FPC = 1:160)

# "1:10", "100001:100012", "1,3,5" or "1:5,9" -> sorted unique positive integers.
parse_ids <- function(x, what = "--pops") {
  if (is.null(x) || !nzchar(x)) stop(what, " is empty", call. = FALSE)
  parts <- strsplit(x, ",", fixed = TRUE)[[1]]
  ids <- unlist(lapply(parts, function(p) {
    ab <- strsplit(trimws(p), ":", fixed = TRUE)[[1]]
    v <- suppressWarnings(as.numeric(ab))
    if (!length(v) || length(v) > 2L || anyNA(v) || any(v != floor(v)) || any(v < 1))
      stop(what, ": cannot read '", p, "' (use a:b or a comma list of positive integers)", call. = FALSE)
    if (length(v) == 2L) { if (v[2] < v[1]) stop(what, ": empty range '", p, "'", call. = FALSE); seq(v[1], v[2]) } else v
  }))
  ids <- as.integer(ids)
  if (anyDuplicated(ids)) stop(what, ": duplicated ids", call. = FALSE)
  sort(ids)
}

# The frozen-design record: the specification's recorded SHA-256 and the SHA-256 of every frozen input (paths
# relative to the repository root).
spec_frozen_record <- function(repo_root) {
  f <- file.path(repo_root, SIM_FROZEN_DESIGN)
  if (!file.exists(f)) return(NULL)
  jsonlite::fromJSON(f, simplifyVector = FALSE)
}

spec_frozen_sha <- function(repo_root) {
  rec <- spec_frozen_record(repo_root)
  sha <- if (is.null(rec) || is.null(rec$specification_sha256)) "" else rec$specification_sha256
  if (grepl("^[0-9a-f]{64}$", sha)) sha else NA_character_
}

# The recorded specification identifier when every frozen input has its recorded SHA-256; otherwise a message that
# names the inputs that differ (and so can never equal the identifier).
spec_sha256 <- function(repo_root) {
  rec <- spec_frozen_record(repo_root)
  if (is.null(rec)) return(NA_character_)
  files <- unlist(rec$frozen_files)
  got <- vapply(names(files), function(p) {
    f <- file.path(repo_root, p)
    if (file.exists(f)) digest::digest(file = f, algo = "sha256") else "missing"
  }, "")
  bad <- names(files)[got != files]
  if (length(bad)) paste("frozen inputs differ:", paste(bad, collapse = ", ")) else spec_frozen_sha(repo_root)
}

spec_freeze_check <- function(repo_root) {
  exp <- spec_frozen_sha(repo_root); act <- spec_sha256(repo_root)
  list(ok = !is.na(exp) && identical(exp, act), expected = exp, actual = act, version = SIM_SPEC_VERSION,
       basis = "config/simulation/frozen-design.json: frozen inputs verified, specification SHA-256 recorded")
}

run_mode_label <- function(mode) switch(mode,
  dev = "DEV: SYNTHETIC root; a code check, not a result",
  dry = "v9.8 DRY RUN (round 3, §12.5): implementation, gates and runtime only; grades nothing; not a confirmatory result",
  prod = "v9.8 PRODUCTION: confirmatory output, graded only by the frozen §6 plan after the run",
  stop("unknown mode ", mode, call. = FALSE))

# Resolves the mode, root, population and rep IDs from the runner's arguments.
#   arg_of(name, default)   the runner's argument reader
#   census_only             TRUE for an R = 0 cell (no rep axis)
#   master_id               the cell's master, to validate production population IDs
run_mode_setup <- function(arg_of, repo_root, B_default = 2L, R_default = 10L, census_only = FALSE, master_id = NULL) {
  mode <- arg_of("--mode", "dev")
  if (!mode %in% c("dev", "dry", "prod")) stop("--mode must be dev, dry or prod", call. = FALSE)
  B <- suppressWarnings(as.integer(arg_of("--B", as.character(B_default))))
  R <- suppressWarnings(as.integer(arg_of("--R", as.character(R_default))))
  if (is.na(B) || B < 1L) stop("--B must be a positive integer", call. = FALSE)
  if (is.na(R) || R < 0L) stop("--R must be a non-negative integer", call. = FALSE)
  pops_arg <- arg_of("--pops", "")
  reps_arg <- arg_of("--reps", "")
  if (mode == "dev") {
    pops <- if (nzchar(pops_arg)) parse_ids(pops_arg) else seq_len(B)
    reps <- if (nzchar(reps_arg)) parse_ids(reps_arg, "--reps") else seq_len(R)
    root <- SIM_DEV_ROOT
    spec <- list(ok = NA, expected = spec_frozen_sha(repo_root), actual = spec_sha256(repo_root), version = SIM_SPEC_VERSION)
  } else {
    spec <- spec_freeze_check(repo_root)
    if (!isTRUE(spec$ok))
      stop(sprintf("--mode %s needs the frozen spec %s: %s (frozen-design record %s: %s)",
                   mode, SIM_SPEC_VERSION, spec$actual, SIM_FROZEN_DESIGN, spec$expected), call. = FALSE)
    if (!nzchar(pops_arg)) stop("--mode ", mode, " needs explicit --pops (production population IDs)", call. = FALSE)
    pops <- parse_ids(pops_arg)
    root <- SIM_PROD_ROOT
    if (census_only) {
      if (mode == "dry" && any(pops <= DRY_ID_BASE))
        stop("a census-only dry-run cell uses the reserved population IDs 1000001-1000005", call. = FALSE)
      if (mode == "prod" && any(pops > DRY_ID_BASE)) stop("production census cells use production population IDs", call. = FALSE)
      reps <- integer(0)
    } else {
      if (any(pops > DRY_ID_BASE)) stop("sampled cells use production population IDs (<= ", DRY_ID_BASE, ")", call. = FALSE)
      reps <- if (nzchar(reps_arg)) parse_ids(reps_arg, "--reps") else if (mode == "dry") DRY_ID_BASE + seq_len(R) else seq_len(R)
      if (mode == "dry" && any(reps <= DRY_ID_BASE))
        stop("dry-run replicates use the disjoint rep-ID range > ", DRY_ID_BASE, " (dry runs and production use disjoint replicate IDs)", call. = FALSE)
      if (mode == "prod" && any(reps > DRY_ID_BASE)) stop("production replicates use rep IDs <= ", DRY_ID_BASE, call. = FALSE)
    }
    if (!is.null(master_id) && !(census_only && mode == "dry")) {
      ok <- PROD_POPS[[master_id]]
      if (is.null(ok)) stop("unknown master ", master_id, call. = FALSE)
      if (!all(pops %in% ok)) stop(sprintf("population IDs outside master %s's production range: %s", master_id,
                                           paste(head(setdiff(pops, ok), 5), collapse = ",")), call. = FALSE)
    }
  }
  # Array tasks: --task i --block k takes the i-th block of k populations.
  task <- suppressWarnings(as.integer(arg_of("--task", "0"))); block <- suppressWarnings(as.integer(arg_of("--block", "0")))
  if (is.na(task) || is.na(block) || task < 0L || block < 0L) stop("--task and --block must be non-negative integers", call. = FALSE)
  if (task > 0L) {
    if (block < 1L) stop("--task needs --block", call. = FALSE)
    first <- (task - 1L) * block + 1L
    if (first > length(pops)) stop(sprintf("--task %d is beyond the %d populations", task, length(pops)), call. = FALSE)
    pops <- pops[first:min(task * block, length(pops))]
  }
  list(mode = mode, root = root, pops = pops, reps = reps, B = length(pops), R = length(reps),
       task = task, block = block, label = run_mode_label(mode), spec = spec, census_only = census_only)
}
