# Build the saved quadrature-rule artifact and its manifest.
#
#   Rscript src/simulation/env/run.R src/simulation/grids/build_quad_rules.R
#
# Writes grids/quad-rules-dev-v1.tsv (every base rule the runs use) and
# grids/MANIFEST.json (its sha256 plus the P0 grid hashes it must reproduce at H = 25).
# The recorded runs use the saved artifact; rebuilding it is not needed to reproduce them.
# After any rebuild, run grids/record_grid_hashes.R to record the hashes at every H
# (25, 100, 400); this script records H = 25 only.
# Quadrature only: no population, sample or estimator is touched.

via_wrapper <- endsWith(getOption("sim.engine.script", ""), "/grids/build_quad_rules.R")
if (!via_wrapper) stop("run through the locked wrapper: Rscript simulation/engine/env/run.R simulation/engine/grids/build_quad_rules.R")
engine_dir <- getOption("sim.engine.dir")
for (f in c("quadrature.R", "quad_rules_io.R", "nuisance_truth.R", "p0_grid.R"))
  source(file.path(engine_dir, "R", f), local = TRUE)

FILE <- "quad-rules-dev-v1.tsv"
# Every base rule the runs draw: primary P0 grid GL(6) x GL(2); independent grid
# CC(9) x CC(5); Gauss-Hermite 64 (eta0, R-a marginalisation, the census) and 128 (U9).
SPEC <- data.frame(rule = c("gauss-legendre", "gauss-legendre", "clenshaw-curtis", "clenshaw-curtis",
                            "gauss-hermite", "gauss-hermite"),
                   n = c(6L, 2L, 9L, 5L, 64L, 128L), stringsAsFactors = FALSE)
H <- 25L
PRIMARY <- list(rule = "gauss-legendre", n_c = 6L, n_psu = 2L)    # p0_context() defaults
SECOND <- list(rule = "clenshaw-curtis", n_c = 9L, n_psu = 5L)

path <- file.path(engine_dir, "grids", FILE)
quad_rules_deactivate()
h_computed <- c(primary = p0_grid_hash(do.call(p0_grid, c(list(H = H), PRIMARY))),
                second = p0_grid_hash(do.call(p0_grid, c(list(H = H), SECOND))))
sha <- quad_rules_save(quad_rules_build(SPEC), path)
ld <- quad_rules_load(path, sha)
quad_rules_activate(ld, required = TRUE)
h_loaded <- c(primary = p0_grid_hash(do.call(p0_grid, c(list(H = H), PRIMARY))),
              second = p0_grid_hash(do.call(p0_grid, c(list(H = H), SECOND))))
quad_rules_deactivate()
if (!identical(h_computed, h_loaded))
  stop("on the saving platform the loaded rules must reproduce the computed grids bit for bit")

man <- list(
  file = FILE, sha256 = sha, status = "DEV -- not a production freeze (sign-off sequence pending)",
  created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  built_on = list(R = R.version.string, platform = R.version$platform, os = utils::sessionInfo()$running,
                  La_version = La_version()),
  rules = SPEC,
  expected_grid_hash = list(H = H, primary = list(definition = PRIMARY, sha256 = unname(h_loaded["primary"])),
                            second = list(definition = SECOND, sha256 = unname(h_loaded["second"]))),
  note = paste("p0_grid_hash() of the H = 25 grids built from this artifact. Any platform that loads the",
               "artifact must reproduce these; tests/testthat/test-quad-rules.R asserts it."))
jsonlite::write_json(man, file.path(engine_dir, "grids", "MANIFEST.json"), auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("wrote %s (sha256 %s)\nprimary grid %s\nsecond grid  %s\n", FILE, sha,
            substr(h_loaded["primary"], 1, 12), substr(h_loaded["second"], 1, 12)))
