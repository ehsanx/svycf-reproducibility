# Sourced by setup.sh, submit.sh and job.sh (not run on its own). Sets SPA_ROOT, loads SPA_MODULES from
# config/slurm.env, and makes every Rscript use only the pinned package library (SPA_R_LIBRARY, by default
# src/application/library) and R's own library: no user or site library, no user profile.
SPA_ROOT="${SPA_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
if [ -f "$SPA_ROOT/config/slurm.env" ]; then
  # shellcheck disable=SC1091
  source "$SPA_ROOT/config/slurm.env"
fi
if [ -n "${SPA_MODULES:-}" ] && command -v module >/dev/null 2>&1; then
  # shellcheck disable=SC2086
  module load $SPA_MODULES
fi
export SPA_R_LIBRARY="${SPA_R_LIBRARY:-$SPA_ROOT/src/application/library}"
export R_LIBS="$SPA_R_LIBRARY"
export R_LIBS_USER="$SPA_ROOT/src/application/no-user-library"     # a path that does not exist is ignored by R
export R_LIBS_SITE="$SPA_ROOT/src/application/no-site-library"
export R_PROFILE_USER="$SPA_ROOT/src/application/no-user-profile"
