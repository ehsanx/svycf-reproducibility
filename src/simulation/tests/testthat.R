# Runs the engine's testthat suite against whatever library is already active.
# It never loads an environment itself: the simulation's locked library is
# brought in explicitly by the caller (renv::load(<engine dir>) or the engine's
# wrapper), and there is no .Rprofile anywhere in this tree.
#
#   Rscript src/simulation/tests/testthat.R [filter]
#   Rscript src/simulation/env/run.R src/simulation/tests/testthat.R [filter]
# Under env/run.R this file is sys.source()d: --file names the wrapper, so the
# path and arguments come from options sim.engine.script / sim.engine.args.
via_wrapper <- endsWith(getOption("sim.engine.script", ""), "/tests/testthat.R")
args <- commandArgs(trailingOnly = FALSE)
self <- if (via_wrapper) getOption("sim.engine.script") else sub("^--file=", "", args[grep("^--file=", args)])
tests_dir <- if (length(self)) {
  file.path(dirname(normalizePath(self[1L], winslash = "/")), "testthat")
} else {
  "tests/testthat"
}
filter <- if (via_wrapper) getOption("sim.engine.args") else commandArgs(trailingOnly = TRUE)
testthat::test_dir(tests_dir, filter = if (length(filter)) filter[1L] else NULL,
                   reporter = "summary", stop_on_failure = TRUE)
