# Record, in grids/MANIFEST.json, the P0 grid hashes that the saved quadrature
# artifact must reproduce at every H the runs use: fixed-strata H = 5, the E2
# FPC sweep H = 50, the E2 vanishing-fraction terminal rung H = 200,
# and many-strata H = 25, 100, 400.
#
#   Rscript src/simulation/env/run.R src/simulation/grids/record_grid_hashes.R
#
# Adds "expected_grid_hashes" (one entry per H) and leaves the artifact, its sha256 and
# the H = 25 "expected_grid_hash" record untouched; the H = 25 hashes computed here
# must equal that record, or nothing is written. The P0 grid depends on H only through
# the K_B class shares (R/p0_grid.R), so H = 100 and H = 400 have the same nodes and
# weights; their hashes differ because the grid definition records H.
# Quadrature only: no population, sample or estimator is touched.

via_wrapper <- endsWith(getOption("sim.engine.script", ""), "/grids/record_grid_hashes.R")
if (!via_wrapper) stop("run through the locked wrapper: Rscript simulation/engine/env/run.R simulation/engine/grids/record_grid_hashes.R")
engine_dir <- getOption("sim.engine.dir")
for (f in c("quadrature.R", "quad_rules_io.R", "nuisance_truth.R", "p0_grid.R"))
  source(file.path(engine_dir, "R", f), local = TRUE)

H_LIST <- c(5L, 25L, 50L, 100L, 200L, 400L, 800L)   # fixed-strata H = 5; E2 FPC sweep H = 50; many-strata m = 50 ... 1600; E2 vanishing-fraction terminal rung H = 200
path <- file.path(engine_dir, "grids", "MANIFEST.json")
man <- jsonlite::fromJSON(path, simplifyVector = FALSE)
dev <- quad_rules_dev(engine_dir)                               # sha256-verified load
quad_rules_activate(dev, required = TRUE)
e <- man$expected_grid_hash
def <- function(d) list(rule = d$rule, n_c = as.integer(d$n_c), n_psu = as.integer(d$n_psu))
hashes <- lapply(H_LIST, function(H) list(
  H = H,
  primary_sha256 = p0_grid_hash(do.call(p0_grid, c(list(H = H), def(e$primary$definition)))),
  second_sha256 = p0_grid_hash(do.call(p0_grid, c(list(H = H), def(e$second$definition))))))
quad_rules_deactivate()
h25 <- hashes[[which(H_LIST == e$H)]]
if (!identical(h25$primary_sha256, e$primary$sha256) || !identical(h25$second_sha256, e$second$sha256))
  stop("the H = ", e$H, " hashes computed now differ from the recorded expected_grid_hash -- nothing written")

man$expected_grid_hashes <- hashes
man$status <- "FROZEN with spec v9.8 (2026-10-11, the R3 design freeze): the base quadrature rules every v9.8 run loads (sha256 unchanged since 2026-10-07)"
man$expected_grid_hashes_recorded <- list(
  utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  R = R.version.string, platform = R.version$platform,
  definitions = "primary and second as in expected_grid_hash (p0_context() defaults)",
  note = paste("p0_grid_hash() at each H, from the saved artifact. Every platform must reproduce them;",
               "tests/testthat/test-quad-rules.R asserts all of them and runs/run_e2e_dev.R checks its own H",
               "as a required check."))
jsonlite::write_json(man, path, auto_unbox = TRUE, pretty = TRUE)
for (h in hashes) cat(sprintf("H = %3d  primary %s  second %s\n", h$H, substr(h$primary_sha256, 1, 12),
                              substr(h$second_sha256, 1, 12)))
