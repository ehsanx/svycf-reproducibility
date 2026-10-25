# Build the saved base rules of the T2 refinement ladder (R/t2_refine.R) and their manifest.
#
#   Rscript src/simulation/env/run.R src/simulation/grids/build_refine_rules.R
#
# Writes grids/quad-rules-refine-v1.tsv and grids/REFINE-MANIFEST.json. The frozen artifact (quad-rules-dev-v1.tsv,
# grids/MANIFEST.json) is NOT touched: R/t2_refine.R merges these rules into the active store only for keys the
# frozen artifact lacks, and stops if a shared key disagrees. The rules: Gauss-Legendre 8, 3, 10, 4, 12 (the ladder
# GL(8,3), GL(10,4), GL(12,6); 6 is already frozen) and Clenshaw-Curtis 13, 7 (U9's brute force on a refined fold,
# CC(13,7)). Saved as exact hex floats, so every platform that loads them builds the same nodes and weights.
# Quadrature only: no population, sample or estimator is touched.

via_wrapper <- endsWith(getOption("sim.engine.script", ""), "/grids/build_refine_rules.R")
if (!via_wrapper) stop("run through the locked wrapper: Rscript simulation/engine/env/run.R simulation/engine/grids/build_refine_rules.R")
engine_dir <- getOption("sim.engine.dir")
for (f in c("quadrature.R", "quad_rules_io.R")) source(file.path(engine_dir, "R", f), local = TRUE)

FILE <- "quad-rules-refine-v1.tsv"
SPEC <- data.frame(rule = c(rep("gauss-legendre", 5), rep("clenshaw-curtis", 2)),
                   n = c(8L, 3L, 10L, 4L, 12L, 13L, 7L), stringsAsFactors = FALSE)
path <- file.path(engine_dir, "grids", FILE)
quad_rules_deactivate()
built <- quad_rules_build(SPEC)
sha <- quad_rules_save(built, path)
ld <- quad_rules_load(path, sha)
same <- vapply(names(built), function(k) identical(built[[k]]$x, ld$rules[[k]]$x) && identical(built[[k]]$w, ld$rules[[k]]$w), TRUE)
if (!all(same)) stop("the saved rules do not reload bit for bit: ", paste(names(same)[!same], collapse = ", "))
wsum <- vapply(ld$rules, function(r) sum(r$w), 0)
man <- list(
  file = FILE, sha256 = sha,
  status = "PROPOSED (round 3d, 2026-10-15): the base rules of the T2 refinement ladder; not part of the frozen v9.8 artifact",
  created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  built_on = list(R = R.version.string, platform = R.version$platform, os = utils::sessionInfo()$running, La_version = La_version()),
  rules = SPEC,
  weight_sums = as.list(wsum),
  note = "Base rules on [-1, 1] (weights sum to 2). R/t2_refine.R t2_refine_rules_extend() adds them to the active store after the frozen artifact.")
jsonlite::write_json(man, file.path(engine_dir, "grids", "REFINE-MANIFEST.json"), auto_unbox = TRUE, pretty = TRUE, digits = I(17))
cat(sprintf("wrote %s (sha256 %s): %s\n", FILE, sha, paste(names(ld$rules), collapse = ", ")))
