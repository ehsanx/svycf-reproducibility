#!/bin/bash
# Install the R packages of the SPA analysis into their own library (src/application/library, or SPA_R_LIBRARY):
# the 62 CRAN packages of config/application/packages.csv at their pinned versions, from the CRAN snapshot of
# 2026-09-24, and the vendored svytmle build. Nothing is read from or written to your own R library. Needs internet
# and a C/C++/Fortran toolchain (packages are built from source); on a cluster, run it on a login node. Run once:
#   bash scripts/application/setup.sh             # install, then load everything as the jobs will
#   bash scripts/application/setup.sh --dry-run   # show what would be installed; install nothing
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
source scripts/application/env.sh
echo "c976da847d9d28069123835b4bb71281  vendor/svytmle_0.1.2.tar.gz" | md5sum -c --quiet - \
  || { echo "vendor/svytmle_0.1.2.tar.gz does not match its recorded md5 -- stop"; exit 1; }
export CRAN_SNAPSHOT="${SPA_CRAN_SNAPSHOT:-https://packagemanager.posit.co/cran/2026-09-24}"
mkdir -p outputs/application
if [ "${1:-}" = "--dry-run" ]; then Rscript scripts/application/install_packages.R --dry-run; exit 0; fi
rm -f outputs/application/SETUP-OK.txt
Rscript scripts/application/install_packages.R
Rscript scripts/application/install_packages.R --verify
echo "setup done (outputs/application/SETUP-OK.txt)"
