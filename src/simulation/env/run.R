# Explicit entry point into the simulation engine's locked environment.
#
#   Rscript src/simulation/env/run.R <script.R> [args...]
#
# There is deliberately no .Rprofile in src/simulation and none at the repository
# root: the SPA analysis (src/application) uses its own package library
# and the simulation is entered only through this
# wrapper. The script runs with the renv project library from
# src/simulation/renv.lock. Its arguments are in getOption("sim.engine.args");
# commandArgs() still shows this wrapper's own command line.
#
# Thread pinning: determinism over
# speed. BLAS/OpenMP libraries read their thread counts when they load, and R's
# BLAS loads at start-up, so setting the variables inside this process would be
# too late. The wrapper therefore sets them and re-executes itself once in a
# fresh Rscript. The child also skips the user's ~/.Rprofile (--no-init-file).
# The
# values are recorded by env/record_runtime.R and in options(sim.engine.threads).

pin <- c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
         BLIS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
         RCPP_PARALLEL_NUM_THREADS = "1", R_DATATABLE_NUM_THREADS = "1",
         # ranger's num.threads = NULL resolves through this variable first
         # (ranger 0.18.0: R_RANGER_NUM_THREADS, then option, then Ncpus, then 2).
         R_RANGER_NUM_THREADS = "1")

self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(self) != 1L) stop("run with Rscript: Rscript simulation/engine/env/run.R <script.R> [args...]")
self <- normalizePath(self, winslash = "/", mustWork = TRUE)
engine_dir <- normalizePath(file.path(dirname(self), ".."), winslash = "/", mustWork = TRUE)
args <- commandArgs(TRUE)
if (!length(args)) stop("usage: Rscript simulation/engine/env/run.R <script.R> [args...]")

if (!identical(Sys.getenv("SIM_ENGINE_PINNED"), "1")) {
  inherited <- Sys.getenv(names(pin), unset = NA)
  do.call(Sys.setenv, c(as.list(pin), SIM_ENGINE_PINNED = "1",
                        SIM_ENGINE_INHERITED_THREADS = paste(names(pin), inherited, sep = "=", collapse = ";")))
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  # system2 quotes the command itself; the arguments are quoted here.
  status <- system2(rscript, c("--no-init-file", shQuote(self), shQuote(args)))
  quit(save = "no", status = status)
}

# --- child: thread variables are now in force from process start ---
options(ranger.num.threads = 1L, sim.engine.threads = as.list(Sys.getenv(names(pin))))

# svytmle is restored from the md5-pinned tarball in env/pkg (see env/snapshot.R).
Sys.setenv(RENV_PATHS_CELLAR = file.path(engine_dir, "env", "pkg"))

# Load renv from the project library itself, so no user-library renv is needed
# and the renv version is the locked one.
plat <- R.version$platform
rmin <- paste0("R-", R.version$major, ".", sub("\\..*$", "", R.version$minor))
cand <- Sys.glob(file.path(engine_dir, "renv", "library", "*", rmin, plat, "renv"))
if (length(cand) != 1L)
  stop("no unique renv in the project library for ", rmin, "/", plat,
       " under ", engine_dir, "/renv/library; restore the lock first (env/README.md)")
invisible(loadNamespace("renv", lib.loc = dirname(cand)))
renv::load(engine_dir, quiet = TRUE)

plib <- normalizePath(renv::paths$library(project = engine_dir), winslash = "/", mustWork = TRUE)
if (!identical(normalizePath(.libPaths()[1], winslash = "/"), plib))
  stop("renv::load did not put the project library first: ", .libPaths()[1])

# Refuse to run outside the lock: every locked package must resolve, on the
# library path renv::load() set, to its locked version. That path is the project
# library and then renv's sandbox of R's own library: a clean renv::restore()
# leaves a recommended package (Matrix, lattice, KernSmooth, codetools) in the
# sandbox when R already ships the locked version, so checking the project
# library alone refuses a correct restore (finding 1).
# The user library is never on this path. Cheap (DESCRIPTION reads only).
lock <- jsonlite::fromJSON(file.path(engine_dir, "renv.lock"), simplifyVector = FALSE)
bad <- character()
for (p in lock$Packages) {
  where <- find.package(p$Package, lib.loc = .libPaths(), quiet = TRUE)
  d <- if (length(where))
    suppressWarnings(utils::packageDescription(p$Package, lib.loc = dirname(where[1L]), fields = "Version"))
  else NA_character_
  if (!identical(d, p$Version)) bad <- c(bad, sprintf("%s (lock %s, found %s)", p$Package, p$Version, d))
}
if (length(bad)) stop("library path does not match renv.lock: ", paste(bad, collapse = "; "))

# The script is sys.source()d, so inside it sys.nframe() > 0 and --file names
# this wrapper. Entry scripts therefore identify themselves through
# getOption("sim.engine.script") and read their arguments from sim.engine.args.
script <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
options(sim.engine.dir = engine_dir, sim.engine.args = args[-1], sim.engine.script = script)
sys.source(script, envir = new.env(parent = globalenv()), chdir = FALSE, toplevel.env = globalenv())
