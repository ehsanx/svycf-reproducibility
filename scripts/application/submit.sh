#!/bin/bash
# Submit the SPA analysis to Slurm. Run from anywhere, after scripts/application/setup.sh and after placing the 24
# approved DHS zips in data/spa/zips/ (or SPA_RAW_DIR/zips).
#   bash scripts/application/submit.sh prep [--print]                    # check, unpack, count audit, self-tests
#   bash scripts/application/submit.sh run <survey>... [--after JOBID] [--print]
# Surveys: MW_SPA_2013-14, HT_SPA_2017-18, HT_SPA_2013 (the recorded order). Per survey: the census overlap, the
# resample chunks (32 cores each) and the reference splits (SPA_REF_CHUNK cores), then collect and summarise.
# Surveys named together run one after another; a failure cancels the later ones.
# Settings from config/slurm.env: SLURM_ACCOUNT, SLURM_PARTITION, SPA_MODULES, SLURM_MAIL, SLURM_EXTRA, and the
# walltimes SPA_RUN_TIME (default 01:00:00) and SPA_REF_TIME (default 05:30:00), about three times the recorded
# run times. Outputs: outputs/application/.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
if [ -f config/slurm.env ]; then
  # shellcheck disable=SC1091
  source config/slurm.env
fi
CMD="${1:-}"; [ $# -gt 0 ] && shift
PRINT=0; AFTER=""; SURVEYS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --print) PRINT=1; shift ;;
    --after) AFTER="${2:?--after needs a job id}"; shift 2 ;;
    MW_SPA_2013-14|HT_SPA_2017-18|HT_SPA_2013) SURVEYS+=("$1"); shift ;;
    *) echo "unknown argument: $1"; exit 2 ;;
  esac
done
SB=(sbatch --parsable)
if [ -n "${SLURM_ACCOUNT:-}" ]; then SB+=(--account="$SLURM_ACCOUNT"); fi
if [ -n "${SLURM_PARTITION:-}" ]; then SB+=(--partition="$SLURM_PARTITION"); fi
if [ -n "${SLURM_MAIL:-}" ]; then SB+=(--mail-user="$SLURM_MAIL" --mail-type=FAIL); fi
read -r -a EXTRA <<< "${SLURM_EXTRA:-}"
SB+=(${EXTRA[@]+"${EXTRA[@]}"})
RAW="${SPA_RAW_DIR:-$ROOT/data/spa}"
EXP="--export=ALL,SPA_ROOT=$ROOT,SPA_RAW_DIR=$RAW"
LOGS=outputs/application/logs
mkdir -p "$LOGS"
submit() {               # prints the command with --print; otherwise submits it and prints the job id, or stops
  if [ "$PRINT" = 1 ]; then printf '%q ' "$@" >&2; echo >&2; echo "PRINTED"; return 0; fi
  local x
  if ! x=$("$@"); then echo "sbatch failed -- nothing more submitted" >&2; return 1; fi
  x="${x%%;*}"
  [[ "$x" =~ ^[0-9]+$ ]] || { echo "sbatch returned no job id ('$x') -- nothing more submitted" >&2; return 1; }
  echo "$x"
}
common=(--nodes=1 --ntasks=1 --chdir="$ROOT")

case "$CMD" in
  prep)
    id=$(submit "${SB[@]}" "${common[@]}" --job-name=spa_prep --cpus-per-task=1 --mem=8G --time=01:00:00 \
         --output="$LOGS/%x_%j.out" --error="$LOGS/%x_%j.err" "$EXP,SPA_JOB=prep" scripts/application/job.sh)
    [ "$PRINT" = 1 ] || echo "submitted the preparation job $id" ;;

  run)
    [ "${#SURVEYS[@]}" -ge 1 ] || { echo "name the survey(s): MW_SPA_2013-14 HT_SPA_2017-18 HT_SPA_2013"; exit 2; }
    grep -q "^# SPA setup OK" outputs/application/SETUP-OK.txt 2>/dev/null || [ "$PRINT" = 1 ] \
      || { echo "run scripts/application/setup.sh first"; exit 1; }
    RUN_TIME="${SPA_RUN_TIME:-01:00:00}"; REF_TIME="${SPA_REF_TIME:-05:30:00}"
    export SPA_REP_CHUNK="${SPA_REP_CHUNK:-64}" SPA_REF_CHUNK="${SPA_REF_CHUNK:-20}"
    source scripts/application/env.sh                              # modules and the pinned package library
    Rscript scripts/application/make_tasks.R > outputs/application/tasks.tsv      # reads the protocol only
    prev="$AFTER"
    for S in "${SURVEYS[@]}"; do
      for m in overlap run reference; do
        awk -F'\t' -v s="$S" -v m="$m" '$2 == s && $1 == m' outputs/application/tasks.tsv > "outputs/application/tasks_${S}_$m.tsv"
      done
      n_o=$(wc -l < "outputs/application/tasks_${S}_overlap.tsv"); n_r=$(wc -l < "outputs/application/tasks_${S}_run.tsv")
      n_f=$(wc -l < "outputs/application/tasks_${S}_reference.tsv")
      dep=(); if [ -n "$prev" ]; then dep=(--dependency=afterok:"$prev" --kill-on-invalid-dep=yes); fi
      out=(--output="$LOGS/%x_%A_%a.out" --error="$LOGS/%x_%A_%a.err")
      jo=$(submit "${SB[@]}" "${common[@]}" ${dep[@]+"${dep[@]}"} "${out[@]}" --job-name=spa_overlap --array=1-"$n_o" \
           --cpus-per-task=1 --mem=8G --time=01:00:00 "$EXP,SPA_JOB=task,SPA_TASKS=$ROOT/outputs/application/tasks_${S}_overlap.tsv" \
           scripts/application/job.sh)
      jr=$(submit "${SB[@]}" "${common[@]}" ${dep[@]+"${dep[@]}"} "${out[@]}" --job-name=spa_resample --array=1-"$n_r" \
           --cpus-per-task=32 --mem=64G --time="$RUN_TIME" "$EXP,SPA_JOB=task,SPA_TASKS=$ROOT/outputs/application/tasks_${S}_run.tsv" \
           scripts/application/job.sh)
      jf=$(submit "${SB[@]}" "${common[@]}" ${dep[@]+"${dep[@]}"} "${out[@]}" --job-name=spa_reference --array=1-"$n_f" \
           --cpus-per-task="$SPA_REF_CHUNK" --mem=32G --time="$REF_TIME" \
           "$EXP,SPA_JOB=task,SPA_TASKS=$ROOT/outputs/application/tasks_${S}_reference.tsv" scripts/application/job.sh)
      prev=$(submit "${SB[@]}" "${common[@]}" --dependency=afterok:"$jo":"$jr":"$jf" --kill-on-invalid-dep=yes \
             --job-name=spa_finalize --cpus-per-task=1 --mem=16G --time=00:30:00 --output="$LOGS/%x_%j.out" \
             --error="$LOGS/%x_%j.err" "$EXP,SPA_JOB=finalize,SPA_SURVEY=$S" scripts/application/job.sh)
      [ "$PRINT" = 1 ] || echo "$S: overlap $jo, resamples $jr ($n_r x 32 cores, $RUN_TIME), reference $jf ($n_f x $SPA_REF_CHUNK cores, $REF_TIME), finalize $prev"
    done ;;

  *) sed -n '2,12p' "$0"; exit 2 ;;
esac
