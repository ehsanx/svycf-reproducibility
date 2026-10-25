# install_packages.R -- the SPA analysis's R packages, pinned and isolated. Installs every package listed in
# config/application/packages.csv at its pinned version -- from the CRAN snapshot of 2026-09-24 (which serves
# exactly those versions; they include every version the recorded runs logged), and svytmle from the vendored
# build vendor/svytmle_0.1.2.tar.gz -- into one library that holds nothing else: SPA_R_LIBRARY, by default
# src/application/library. The user and site libraries are neither read nor written.
#   Rscript scripts/application/install_packages.R             install what is missing or at another version
#   Rscript scripts/application/install_packages.R --dry-run   print the plan; install nothing
#   Rscript scripts/application/install_packages.R --verify    fresh session: load every package from that library,
#                                                              check every version and the svytmle build stamp, and
#                                                              write outputs/application/SETUP-OK.txt
# Run from the repository root (scripts/application/setup.sh does).
args <- commandArgs(trailingOnly = TRUE)
mode <- if ("--dry-run" %in% args) "dry" else if ("--verify" %in% args) "verify" else "install"
SNAP <- Sys.getenv("CRAN_SNAPSHOT", "https://packagemanager.posit.co/cran/2026-09-24")
TARBALL <- file.path("vendor", "svytmle_0.1.2.tar.gz")
LIB <- Sys.getenv("SPA_R_LIBRARY", file.path("src", "application", "library"))
PINS <- utils::read.csv(file.path("config", "application", "packages.csv"), colClasses = "character")
pin <- grep("^\\s*svytmle_packaged:", readLines(file.path("config", "application", "protocol.yaml"), warn = FALSE), value = TRUE)
if (length(pin) != 1) stop("the protocol must have exactly one svytmle_packaged line", call. = FALSE)
PINNED <- sub('^\\s*svytmle_packaged:\\s*"([^"]+)".*$', "\\1", pin)
CRAN <- PINS[PINS$source != "vendor", , drop = FALSE]

dir.create(LIB, recursive = TRUE, showWarnings = FALSE)
LIB <- normalizePath(LIB, winslash = "/")
.libPaths(c(LIB, .Library))                                    # this library first, then R's own; nothing else
same_version <- function(a, b) !is.na(a) && package_version(a) == package_version(b)   # "2.0-42" == "2.0.42"
installed <- function() installed.packages(lib.loc = LIB, fields = "Packaged")
version_in_lib <- function(ip, p) { i <- match(p, ip[, "Package"]); if (is.na(i)) NA_character_ else unname(ip[i, "Version"]) }
svytmle_ok <- function(ip) { i <- match("svytmle", ip[, "Package"]); !is.na(i) && identical(unname(ip[i, "Packaged"]), PINNED) }

if (mode == "verify") {
  rows <- vapply(PINS$package, function(p) {
    ns <- tryCatch(loadNamespace(p, lib.loc = LIB), error = function(e) stop("cannot load ", p, ": ", conditionMessage(e), call. = FALSE))
    v <- as.character(packageVersion(p, lib.loc = LIB)); w <- PINS$version[PINS$package == p]
    if (!same_version(v, w)) stop(p, " is ", v, ", pinned ", w, call. = FALSE)
    where <- normalizePath(dirname(getNamespaceInfo(ns, "path")), winslash = "/")
    if (!identical(where, LIB)) stop(p, " loads from ", where, ", not from ", LIB, call. = FALSE)
    sprintf("%-14s %-12s", p, v)
  }, "")
  pk <- utils::packageDescription("svytmle", lib.loc = LIB)$Packaged
  if (!identical(pk, PINNED)) stop("svytmle is not the vendored build (Packaged '", pk, "', pinned '", PINNED, "')", call. = FALSE)
  out <- c(sprintf("# SPA setup OK %s; %s; library %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                   R.version.string, LIB),
           sprintf("# %d pinned packages at their pinned versions; svytmle Packaged: %s", nrow(PINS), pk), rows)
  dir.create(file.path("outputs", "application"), recursive = TRUE, showWarnings = FALSE)
  writeLines(out, file.path("outputs", "application", "SETUP-OK.txt"))
  writeLines(out[1:2])
  quit(save = "no")
}

ip <- installed()
todo <- CRAN$package[!vapply(seq_len(nrow(CRAN)), function(i) same_version(version_in_lib(ip, CRAN$package[i]), CRAN$version[i]), TRUE)]
av <- available.packages(repos = SNAP, type = "source")
missing_in_snapshot <- setdiff(todo, rownames(av))
if (length(missing_in_snapshot)) stop("not in the snapshot: ", paste(missing_in_snapshot, collapse = ", "), call. = FALSE)
wrong <- todo[!vapply(todo, function(p) same_version(av[p, "Version"], CRAN$version[CRAN$package == p]), TRUE)]
if (length(wrong)) stop("the snapshot serves other versions of: ", paste(wrong, collapse = ", "), call. = FALSE)
cat("library:", LIB, "\n")
cat(sprintf("pinned packages already in the library: %d of %d\n", nrow(CRAN) - length(todo), nrow(CRAN)))
cat(sprintf("to install from %s (%d):\n", SNAP, length(todo)))
for (p in todo) cat(sprintf("  %-14s %s\n", p, CRAN$version[CRAN$package == p]))
cat("svytmle (vendored build, Packaged ", PINNED, "): ", if (svytmle_ok(ip)) "already in the library" else "to install", "\n", sep = "")
if (mode == "dry") quit(save = "no")

if (length(todo)) install.packages(todo, lib = LIB, repos = SNAP, type = "source", dependencies = FALSE, Ncpus = 2)
if (!svytmle_ok(installed())) install.packages(TARBALL, lib = LIB, repos = NULL, type = "source")
ip <- installed()
bad <- PINS$package[!vapply(seq_len(nrow(PINS)), function(i) same_version(version_in_lib(ip, PINS$package[i]), PINS$version[i]), TRUE)]
if (length(bad)) stop("after installing, these are missing or at another version: ", paste(bad, collapse = ", "), call. = FALSE)
if (!svytmle_ok(ip)) stop("svytmle in the library is not the vendored build", call. = FALSE)
cat("installs done\n")
