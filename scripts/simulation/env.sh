# Sourced by setup.sh and job.sh (not run on its own). Sets ROOT, ENGINE and the environment of the recorded runs:
# the simulation's own locked library only (no user or site library), svytmle from vendor/, one BLAS/OpenMP thread.
ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
ENGINE="$ROOT/src/simulation"
if [ -f "$ROOT/config/slurm.env" ]; then
  # shellcheck disable=SC1091
  source "$ROOT/config/slurm.env"
fi
if [ -n "${SIM_MODULES:-}" ] && command -v module >/dev/null 2>&1; then
  # shellcheck disable=SC2086
  module load $SIM_MODULES
fi
unset R_LIBS
export R_LIBS_USER="$ENGINE/renv/no-user-library"      # a path that does not exist is ignored by R
export R_LIBS_SITE="$ENGINE/renv/no-site-library"
export R_PROFILE_USER="$ENGINE/renv/no-user-profile"
export RENV_CONFIG_CACHE_ENABLED=FALSE
export RENV_PATHS_SANDBOX="$ENGINE/renv/sandbox"
export RENV_PATHS_CELLAR="$ROOT/vendor"
export RENV_PROJECT="$ENGINE"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 R_RANGER_NUM_THREADS=1 \
       R_DATATABLE_NUM_THREADS=1
