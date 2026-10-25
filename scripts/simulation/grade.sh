#!/bin/bash
# The summaries of the production output (src/simulation/runs/run_grade.R), on one core. Submitted by
# `scripts/simulation/submit.sh grade`, which passes run_grade.R's arguments: --in, --round3, --assembly-roots, --cells
# and --out-root. Exit status 0: every cell complete, nothing withheld, every provenance check met; 3: a cell
# incomplete or a statistic withheld; 4: a provenance check not met; 5: an error in some cell.
set -euo pipefail
: "${SIM_ROOT:?SIM_ROOT is not set -- submit with scripts/simulation/submit.sh grade}"
ROOT="$SIM_ROOT"
cd "$ROOT"
source scripts/simulation/env.sh
Rscript "$ENGINE/env/run.R" "$ENGINE/runs/run_grade.R" "$@"
