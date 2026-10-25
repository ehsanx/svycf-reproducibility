#!/bin/bash
# One job of the SPA analysis; submitted by scripts/application/submit.sh, which passes every Slurm setting.
#   SPA_JOB=prep      check the downloaded zips, unpack them into a private folder, count audit, synthetic self-tests
#   SPA_JOB=task      one line of a task table (SPA_TASKS; line = array index): overlap, reference splits or resamples
#   SPA_JOB=finalize  collect and summarise one survey (SPA_SURVEY)
# Record-level data stay in SPA_RAW_DIR (default data/spa, mode 700) and in memory; outputs are aggregates.
set -euo pipefail
: "${SPA_ROOT:?SPA_ROOT is not set -- submit with scripts/application/submit.sh}"
cd "$SPA_ROOT"
source scripts/application/env.sh                                  # modules and the pinned package library
grep -q "^# SPA setup OK" outputs/application/SETUP-OK.txt 2>/dev/null \
  || { echo "run scripts/application/setup.sh first (no outputs/application/SETUP-OK.txt)"; exit 1; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1      # one BLAS thread per forked worker
export SPA_RAW_DIR="${SPA_RAW_DIR:-$SPA_ROOT/data/spa}"
export SPA_RESULTS_DIR="${SPA_RESULTS_DIR:-$SPA_ROOT/outputs/application}"
private() { local m; m=$(stat -c %a "$SPA_RAW_DIR"); [ $(( 8#$m & 8#077 )) -eq 0 ] \
            || { echo "$SPA_RAW_DIR (mode $m) allows group or other access -- refusing"; exit 1; }; }
echo "[spa] $SPA_JOB ${SPA_SURVEY:-} task ${SLURM_ARRAY_TASK_ID:--} on $(hostname) $(date -u +%FT%TZ)"
case "$SPA_JOB" in
  prep)
    mkdir -p "$SPA_RAW_DIR/zips"; chmod 700 "$SPA_RAW_DIR"; private
    bad=0
    while IFS=, read -r file _survey _type _structure _bytes md5; do
      [ "$file" = file ] && continue
      if [ ! -f "$SPA_RAW_DIR/zips/$file" ]; then echo "missing $SPA_RAW_DIR/zips/$file"; bad=1
      elif [ "$(md5sum < "$SPA_RAW_DIR/zips/$file" | cut -d' ' -f1)" != "$md5" ]; then echo "md5 differs: $file"; bad=1; fi
    done < config/application/data-files.csv
    [ "$bad" -eq 0 ] || { echo "the downloaded zips do not match config/application/data-files.csv -- stop"; exit 1; }
    Rscript src/application/R/00-unpack.R
    chmod -R go-rwx "$SPA_RAW_DIR"
    Rscript src/application/R/01-count-audit.R --selftest
    for S in MW_SPA_2013-14 HT_SPA_2017-18 HT_SPA_2013; do Rscript src/application/R/01-count-audit.R --count --survey "$S"; done
    SPA_WORKERS=1 Rscript src/application/R/02-resample.R --selftest ;;
  task)
    private
    : "${SPA_TASKS:?SPA_TASKS is not set}"
    read -r MODE SURVEY F RANGE <<< "$(sed -n "${SLURM_ARRAY_TASK_ID:?}p" "$SPA_TASKS")"
    case "$MODE" in
      overlap)   Rscript src/application/R/02-resample.R --overlap   --survey "$SURVEY" ;;
      reference) Rscript src/application/R/02-resample.R --reference --survey "$SURVEY" --splits "$RANGE" ;;
      run)       Rscript src/application/R/02-resample.R --run       --survey "$SURVEY" --f "$F" --reps "$RANGE" ;;
      *) echo "unknown mode $MODE"; exit 2 ;;
    esac ;;
  finalize)
    private
    : "${SPA_SURVEY:?SPA_SURVEY is not set}"
    Rscript src/application/R/02-resample.R --collect   --survey "$SPA_SURVEY"
    Rscript src/application/R/02-resample.R --summarise --survey "$SPA_SURVEY" ;;
  *) echo "SPA_JOB must be prep, task or finalize"; exit 2 ;;
esac
echo "[spa] done $(date -u +%FT%TZ)"
