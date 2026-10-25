#!/bin/bash
# Restore the simulation's locked R library (src/simulation/renv.lock) into src/simulation/renv/library, with
# svytmle from vendor/, then check every package against the lock. Needs internet access and a C, C++ and Fortran
# toolchain (packages are built from source); on a cluster, run it on a login node. Run once:
#   bash scripts/simulation/setup.sh
# A rerun installs nothing and only repeats the check. It also writes PACKAGE-COMMIT.txt, the package identity that
# every simulation run records (see config/execution-identity.json).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
source scripts/simulation/env.sh
mkdir -p outputs/simulation
echo "c976da847d9d28069123835b4bb71281  vendor/svytmle_0.1.2.tar.gz" | md5sum -c --quiet - \
  || { echo "vendor/svytmle_0.1.2.tar.gz does not match its recorded md5 -- stop"; exit 1; }
command -v cmake >/dev/null 2>&1 || echo "note: cmake not found; nloptr (needed by lme4) may fail to build"
echo "R: $(Rscript --version 2>&1 | head -n 1)"
Rscript --no-init-file -e '
source(file.path(Sys.getenv("RENV_PROJECT"), "renv", "activate.R"))
if (packageVersion("renv") != "1.2.3") stop("renv is ", packageVersion("renv"), ", the lock records 1.2.3")
options(renv.config.install.jobs = 2L)
renv::restore(project = Sys.getenv("RENV_PROJECT"), prompt = FALSE)
d <- packageDescription("svytmle")
if (!identical(d$Packaged, "2026-09-24 19:40:00 UTC; wilds")) stop("svytmle is not the vendored build: Packaged ", d$Packaged)
cat("svytmle", d$Version, "Packaged:", d$Packaged, "\n")'
set +e
Rscript --no-init-file "$ENGINE/env/check_lock.R" | tee outputs/simulation/setup-check.txt
rc=${PIPESTATUS[0]}
set -e
if grep -q "MISMATCH" outputs/simulation/setup-check.txt \
   || ! grep -q "^svytmle 0.1.2 Packaged: 2026-09-24 19:40:00 UTC; wilds" outputs/simulation/setup-check.txt; then
  echo "setup FAILED: the library does not match the lock (outputs/simulation/setup-check.txt)"; exit 1
fi
if [ "$rc" -ne 0 ]; then
  echo "note: this R version is not one of the two used for the recorded runs (4.4.0 on the cluster, 4.5.1 locally)."
  echo "      The locked package versions are installed; numerical results may differ in the last digits."
fi
id=$(sed -n 's/.*"release_id": *"\([^"]*\)".*/\1/p' config/execution-identity.json)
[ -n "$id" ] || { echo "config/execution-identity.json has no release_id"; exit 1; }
printf 'release:%s\n' "$id" > PACKAGE-COMMIT.txt
echo "setup done; runs will record the package identity release:$id"
