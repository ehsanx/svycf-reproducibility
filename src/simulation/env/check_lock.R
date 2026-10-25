# Validate the restored engine library against renv.lock PACKAGE BY PACKAGE, and check
# the R version SEPARATELY against the R versions used for the recorded runs.
#
#   Rscript --no-init-file src/simulation/env/check_lock.R
#
# The lock's package versions are binding on every
# platform. The R version is not one value: the lock was taken under R 4.5.1 (Windows),
# and every recorded cluster run used R 4.4.0 (Linux); both are listed in APPROVED_R.
# renv::status() reports any R difference as "out of sync", so it cannot be
# this check. Each run records its R version.
# Exit status 0 only if every locked package resolves to its locked version AND the R
# version is one of those two.

APPROVED_R <- c("4.5.1" = "laptop; the R version recorded in renv.lock",
                "4.4.0" = "UBC ARC Sockeye (gcc/9.4.0 r/4.4.0); approved by the author 2026-10-06")

engine_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) != 1L) stop("run with Rscript so the script path is known")
  normalizePath(file.path(dirname(f), ".."), winslash = "/", mustWork = TRUE)
})

# renv from the project library itself, as env/run.R does.
plat <- R.version$platform
rmin <- paste0("R-", R.version$major, ".", sub("\\..*$", "", R.version$minor))
cand <- Sys.glob(file.path(engine_dir, "renv", "library", "*", rmin, plat, "renv"))
if (length(cand) != 1L)
  stop("no unique renv in the project library for ", rmin, "/", plat, "; restore the lock first")
invisible(loadNamespace("renv", lib.loc = dirname(cand)))
suppressMessages(renv::load(engine_dir, quiet = TRUE))

lock <- jsonlite::fromJSON(file.path(engine_dir, "renv.lock"), simplifyVector = FALSE)
bad <- character()
for (p in lock$Packages) {
  where <- find.package(p$Package, lib.loc = .libPaths(), quiet = TRUE)
  d <- if (length(where))
    suppressWarnings(utils::packageDescription(p$Package, lib.loc = dirname(where[1L]), fields = "Version"))
  else NA_character_
  if (!identical(d, p$Version)) bad <- c(bad, sprintf("%s (lock %s, found %s)", p$Package, p$Version, d))
}
rv <- paste(R.version$major, R.version$minor, sep = ".")
r_ok <- rv %in% names(APPROVED_R)
cat(sprintf("packages: %d of %d locked packages at their locked version%s\n",
            length(lock$Packages) - length(bad), length(lock$Packages),
            if (length(bad)) paste0(" -- MISMATCH: ", paste(bad, collapse = "; ")) else ""))
cat(sprintf("R: %s (lock records %s) -- %s\n", rv, lock$R$Version,
            if (r_ok) paste("APPROVED:", APPROVED_R[[rv]]) else
              paste("NOT APPROVED; approved runtimes:", paste(names(APPROVED_R), collapse = ", "))))
d <- utils::packageDescription("svytmle")
cat("svytmle", d$Version, "Packaged:", d$Packaged, "\n")
if (!identical(d$Packaged, "2026-09-24 19:40:00 UTC; wilds")) bad <- c(bad, "svytmle Packaged stamp")
quit(save = "no", status = if (length(bad) || !r_ok) 1L else 0L)
