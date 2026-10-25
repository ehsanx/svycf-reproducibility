#!/bin/bash
# One simulation cell, or one array task of a cell. Submitted by scripts/simulation/submit.sh (which passes every
# Slurm setting on the sbatch command line), or run directly by `submit.sh local`. Reads the cell's line from the
# stage's cell table:   id  minutes  tests  script  args  [mem  tasks  block  after ...]
# Output: outputs/simulation/<stage>/<cell>/<run folder>/. One core; the engine pins every thread count to 1.
set -euo pipefail
: "${SIM_ROOT:?SIM_ROOT is not set -- submit with scripts/simulation/submit.sh}"
: "${SIM_STAGE:?SIM_STAGE is not set}"
: "${SIM_CELL:?SIM_CELL is not set}"
: "${SIM_CELLS:?SIM_CELLS is not set}"
ROOT="$SIM_ROOT"
cd "$ROOT"
source scripts/simulation/env.sh
LINE=$(awk -F'\t' -v id="$SIM_CELL" 'NR > 1 && $1 == id' "$SIM_CELLS")
[ -n "$LINE" ] || { echo "[sim] no cell '$SIM_CELL' in $SIM_CELLS"; exit 2; }
IFS=$'\t' read -r ID MIN TESTS SCRIPT ARGS MEM TASKS BLOCK AFTER REST <<< "$LINE"
TASK_ID="${SIM_TASK:-${SLURM_ARRAY_TASK_ID:-}}"
if [ -n "$TASK_ID" ] && [ "${BLOCK:-0}" -gt 0 ] 2>/dev/null; then ARGS="$ARGS --task $TASK_ID --block $BLOCK"; fi
case "$SCRIPT" in
  slice)  RSCRIPT="$ENGINE/runs/run_e2e_dev.R" ;;
  census) RSCRIPT="$ENGINE/runs/run_e2_census.R" ;;
  e4)     RSCRIPT="$ENGINE/runs/run_e4.R" ;;
  prep)   RSCRIPT="$ENGINE/runs/run_prep.R" ;;
  *) echo "[sim] unknown script '$SCRIPT'"; exit 2 ;;
esac
OUT="$ROOT/outputs/simulation/$SIM_STAGE/$ID"
mkdir -p "$OUT"
echo "[sim] $SIM_STAGE cell $ID${TASK_ID:+ task $TASK_ID}: $SCRIPT $ARGS ($(date -u +%FT%TZ))"
T0=$(date +%s)
# shellcheck disable=SC2086
Rscript "$ENGINE/env/run.R" "$RSCRIPT" $ARGS --out-root "$OUT"
echo "[sim] cell $ID${TASK_ID:+ task $TASK_ID} done in $(( $(date +%s) - T0 )) s"
