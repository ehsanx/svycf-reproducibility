# ---------------------------------------------------------------------------
# 00-unpack.R -- file the DHS Program zips into a fixed layout and unzip them.
#
# Usage (from the repo root):
#   Rscript src/application/R/00-unpack.R                 # zips already in data/spa/zips/
#   Rscript src/application/R/00-unpack.R --from <dir>    # also MOVE DHS zips from <dir>
#                                                         # (e.g. your Downloads folder)
#
# Layout it produces (data/ is ignored by git):
#   data/spa/zips/<FILE>.zip                  the downloads, untouched
#   data/spa/<survey>/<recode|raw>/<FILE>/    unzipped contents
#   data/spa/MANIFEST.txt                     name, bytes, md5 of every zip
#
# DATA POLICY. DHS Terms of Use: the datasets "must not be passed on to other
# researchers without the written consent of The DHS Program." The script ends
# (in a git checkout) by asking git whether any data file under data/ could be staged,
# and stops with an error if one could.
# ---------------------------------------------------------------------------

DHS <- local({
  for (cand in c(file.path("src", "application"), ".", "..")) {
    if (file.exists(file.path(cand, "R", "00-unpack.R"))) return(normalizePath(cand))
  }
  stop("cannot locate src/application/ (run from the repository root)")
})
RAW  <- Sys.getenv("SPA_RAW_DIR", file.path(DHS, "..", "..", "data", "spa"))   # a private (chmod 700) directory
ZIPS <- file.path(RAW, "zips")
dir.create(ZIPS, recursive = TRUE, showWarnings = FALSE)

# The SPA surveys of the study. Key = country code + recode version
# as it appears in the file name. DHS naming is CC DD VV FF DS .zip: country, file
# type, release version, format (DT Stata, SD SAS, SV SPSS, FL flat) and, for SPA
# only, data structure (SR = SPA recode, SP = SPA raw) -- dhsprogram.com/data/
# File-Types-and-Names.cfm. Recode and raw carry different version codes.
SPA_SURVEYS <- c(
  MW6I = "MW_SPA_2013-14", MW6J = "MW_SPA_2013-14", MW6K = "MW_SPA_2013-14",
  HT6A = "HT_SPA_2013",    HT6B = "HT_SPA_2013",
  HT7A = "HT_SPA_2017-18"
)

# --- optionally pull zips in from another folder ---------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2 && args[1] == "--from") {
  src <- normalizePath(args[2], mustWork = TRUE)
  found <- list.files(src, pattern = "^(HT|MW)[A-Z]{2}[0-9][0-9A-Z](DT|SD|SV|FL)(SP|SR)?\\.zip$",
                      ignore.case = TRUE, full.names = TRUE)
  for (f in found) {
    dest <- file.path(ZIPS, basename(f))
    if (file.exists(dest)) { message("already in raw/zips, leaving source alone: ", basename(f)); next }
    if (!file.rename(f, dest)) {                       # cross-drive: copy then remove
      stopifnot(file.copy(f, dest)); unlink(f)
    }
    message("moved ", basename(f), " -> raw/zips/")
  }
}

# --- classify and unzip ------------------------------------------------------
zips <- list.files(ZIPS, pattern = "\\.zip$", ignore.case = TRUE, full.names = TRUE)
if (!length(zips)) stop("no zips in ", ZIPS, " -- download them there, or pass --from <dir>")

rx <- "^([A-Z]{2})([A-Z]{2})([0-9][0-9A-Z])(DT|SD|SV|FL)(SP|SR)?\\.ZIP$"
for (z in zips) {
  nm <- toupper(basename(z))
  m  <- regmatches(nm, regexec(rx, nm))[[1]]
  if (!length(m)) { warning("unrecognised file name, skipped: ", nm); next }
  cc <- m[2]; ver <- m[4]; kind <- m[6]
  survey <- if (nzchar(kind)) {
    SPA_SURVEYS[paste0(cc, ver)] %||% paste0(cc, "_SPA_v", ver)
  } else paste0(cc, "_DHS_v", ver)                      # household-survey recodes
  sub  <- if (kind == "SP") "raw" else "recode"              # SR = recode, SP = raw
  dest <- file.path(RAW, survey, sub, sub("\\.ZIP$", "", nm))
  if (dir.exists(dest) && length(list.files(dest))) next
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(z, exdir = dest)
  message(sprintf("%-16s -> raw/%s/%s/", nm, survey, sub))
}

# --- manifest ------------------------------------------------------------------
info <- file.info(zips)
man  <- data.frame(file = basename(zips), bytes = info$size,
                   md5 = unname(tools::md5sum(zips)), stringsAsFactors = FALSE)
writeLines(c(sprintf("# DHS zips as downloaded; %d files", nrow(man)),
             sprintf("%-20s %12s  %s", man$file, format(man$bytes, big.mark = ","), man$md5)),
           file.path(RAW, "MANIFEST.txt"))

# --- privacy check: no data file may be stageable by git ------------------------
# A copy without git (with a code manifest) is not a repository, so nothing in it
# can be committed; the check applies to git checkouts only.
repo <- normalizePath(file.path(DHS, "..", ".."))
if (!file.exists(file.path(repo, ".git")) && (file.exists(file.path(repo, "CODE-MANIFEST.tsv")) ||
                                               file.exists(file.path(repo, "config", "application", "code-manifest.tsv")))) {
  message("privacy check not applicable: packaged copy without git (nothing here can be committed)")
} else {
  st <- suppressWarnings(system2("git", c("-C", shQuote(repo), "status", "--porcelain", "--untracked-files=all",
                                          "--", "data"), stdout = TRUE, stderr = FALSE))
  if (!is.null(attr(st, "status")) && attr(st, "status") != 0)
    stop("PRIVACY CHECK COULD NOT RUN (git failed); treat every file under application/dhs/ as unverified.")
  paths <- sub("^.. ", "", st)
  ok <- rep(FALSE, length(paths))                         # nothing under data/ may be visible to git
  if (any(!ok)) {
    stop("PRIVACY CHECK FAILED -- these files under application/dhs/ are visible to git:\n  ",
         paste(paths[!ok], collapse = "\n  "),
         "\nExtend the DHS block in .gitignore before doing anything else.")
  }
  message("privacy check passed: no data file under application/dhs/ is visible to git")
}
