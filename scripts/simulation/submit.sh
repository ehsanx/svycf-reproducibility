#!/bin/bash
# Submit the simulation study to Slurm, or run one cell without Slurm. Run from anywhere, after setup.sh.
#   bash scripts/simulation/submit.sh list prep|production          # cells, limits, tasks and core-hours; submits nothing
#   bash scripts/simulation/submit.sh prep [--only a,b] [--print] [--fresh-kappa]
#   bash scripts/simulation/submit.sh production [--only a,b] [--then-grade] [--print]
#   bash scripts/simulation/submit.sh grade [--after-production | --after ID[:ID...]] [--in FOLDER]
#                                           [--preparation recorded|outputs] [--print]
#   bash scripts/simulation/submit.sh local prep|production <cell> [<task>]
# Slurm settings come from config/slurm.env (copy config/slurm.env.example): SLURM_ACCOUNT, SLURM_PARTITION,
# SIM_MODULES, SLURM_MAIL and SLURM_EXTRA. --print shows every sbatch command without submitting.
# Cell tables: config/simulation/preparation-cells.tsv and config/simulation/production-cells.tsv (id, minutes, tests,
# script, args, mem, tasks, block, after). Production time limits come from config/simulation/production-limits.csv,
# the limits of the recorded run. Every submission is recorded in outputs/simulation/logs/<stage>/jobs.tsv, one batch
# per call. The summaries start only after every job of a production batch has ended: --then-grade queues them with
# the batch, --after-production waits for every job of the latest batch. Outputs: outputs/simulation/<stage>/<cell>/.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
if [ -f config/slurm.env ]; then
  # shellcheck disable=SC1091
  source config/slurm.env
fi
STAGE="${1:-}"; [ $# -gt 0 ] && shift
KAPPA=src/simulation/runs/frozen/kappa_ae.json
KAPPA_REC=results/simulation/preparation/prep-panel/prep-panel-20261012T081105Z/kappa_ae.json
KAPPA_SHA=ca5748d394915b72ef911482ccb3f44e03456d8d5246f29f20846cd2f263c359
LIMITS=config/simulation/production-limits.csv
cells_of() { case "$1" in prep) echo config/simulation/preparation-cells.tsv ;;
                          production) echo config/simulation/production-cells.tsv ;;
                          *) echo "unknown stage '$1' (prep or production)" >&2; return 2 ;; esac; }
col() { awk -F'\t' -v id="$2" -v c="$3" -v d="$4" 'NR > 1 && $1 == id {v = $c; print (v == "" ? d : v)}' "$1"; }
limit_of() {             # stage, cell: production uses the recorded limits; preparation the cell table's minutes
  local l=""
  if [ "$1" = production ]; then l=$(awk -F, -v id="$2" 'NR > 1 && $1 == id {print $2}' "$LIMITS"); fi
  if [ -n "$l" ]; then echo "$l"; else col "$(cells_of "$1")" "$2" 2 10; fi
}
sha() { if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }
need_setup() {
  grep -q '^release:' PACKAGE-COMMIT.txt 2>/dev/null || { echo "run scripts/simulation/setup.sh first (no PACKAGE-COMMIT.txt)"; exit 1; }
}
install_kappa() {        # the recorded kappa_AE record, the frozen input of the E3 control cell
  [ "$(sha "$KAPPA_REC")" = "$KAPPA_SHA" ] || { echo "$KAPPA_REC does not match its recorded SHA-256 -- stop"; exit 1; }
  if [ ! -f "$KAPPA" ]; then
    mkdir -p "$(dirname "$KAPPA")"; cp "$KAPPA_REC" "$KAPPA"; echo "installed the recorded kappa_AE record at $KAPPA"
  elif [ "$(sha "$KAPPA")" != "$KAPPA_SHA" ]; then
    echo "$KAPPA differs from the recorded record (a re-derived one?). Production uses the recorded record:"
    echo "move $KAPPA aside and rerun, or compare it with $KAPPA_REC first."; exit 1
  fi
}
list_cells() {
  local T; T=$(cells_of "$1"); local total=0
  printf '  %-18s %6s %5s %5s  %s\n' cell min mem tasks args
  while IFS=$'\t' read -r id min _tests script args mem tasks _block _after _rest; do
    [ "$id" = id ] && continue
    min=$(limit_of "$1" "$id")
    printf '  %-18s %6s %5s %5s  %s %s\n' "$id" "$min" "${mem:-8G}" "${tasks:-1}" "$script" "$args"
    total=$(( total + min * ${tasks:-1} ))
  done < "$T"
  awk -v t="$total" 'BEGIN {printf "  at most %.1f core-hours if every task runs to its limit (one core per task)\n", t / 60}'
}

SB=(sbatch --parsable)
if [ -n "${SLURM_ACCOUNT:-}" ]; then SB+=(--account="$SLURM_ACCOUNT"); fi
if [ -n "${SLURM_PARTITION:-}" ]; then SB+=(--partition="$SLURM_PARTITION"); fi
if [ -n "${SLURM_MAIL:-}" ]; then SB+=(--mail-user="$SLURM_MAIL" --mail-type=FAIL); fi
read -r -a EXTRA <<< "${SLURM_EXTRA:-}"
SB+=(${EXTRA[@]+"${EXTRA[@]}"})
PRINT=0; ONLY=""; AFTER=""; AFTER_PROD=0; THEN_GRADE=0; PREP=recorded; FRESH=0; IN=outputs/simulation/production
parse_opts() {
  while [ $# -gt 0 ]; do
    case "$1" in
      --print) PRINT=1; shift ;;
      --only) ONLY="${2:?--only needs a comma-separated list of cells}"; shift 2 ;;
      --after) AFTER="${2:?--after needs one or more job ids, separated by colons}"; AFTER="${AFTER//,/:}"; shift 2 ;;
      --after-production) AFTER_PROD=1; shift ;;
      --then-grade) THEN_GRADE=1; shift ;;
      --preparation) PREP="${2:?--preparation needs recorded or outputs}"; shift 2 ;;
      --fresh-kappa) FRESH=1; shift ;;
      --in) IN="${2:?--in needs a folder}"; shift 2 ;;
      *) echo "unknown option: $1"; exit 2 ;;
    esac
  done
}
submit() {               # prints the command with --print; otherwise submits it and prints the job id, or stops
  if [ "$PRINT" = 1 ]; then printf '%q ' "$@" >&2; echo >&2; echo "PRINTED"; return 0; fi
  local x
  if ! x=$("$@"); then echo "sbatch failed -- nothing more submitted" >&2; return 1; fi
  x="${x%%;*}"
  [[ "$x" =~ ^[0-9]+$ ]] || { echo "sbatch returned no job id ('$x') -- nothing more submitted" >&2; return 1; }
  echo "$x"
}
submit_grade() {         # $1: colon-separated job ids the summaries wait for (all of them must have ended), or empty
  case "$PREP" in recorded) R3=results/simulation/preparation ;; outputs) R3=outputs/simulation/prep ;;
                  *) echo "--preparation recorded|outputs"; exit 2 ;; esac
  case "$IN" in /*) ;; *) IN="$ROOT/$IN" ;; esac
  local args=(--in "$IN" --round3 "$ROOT/$R3" --assembly-roots "$ROOT/results/simulation/execution-check"
              --cells "$ROOT/config/simulation/production-cells.tsv" --out-root "$ROOT/outputs/simulation/grading")
  mkdir -p outputs/simulation/logs/grading
  local opts=(--job-name=sim_grade --time=300 --mem=8G --cpus-per-task=1 --ntasks=1 --nodes=1 --chdir="$ROOT"
              --export="ALL,SIM_ROOT=$ROOT"
              --output=outputs/simulation/logs/grading/%x_%j.out --error=outputs/simulation/logs/grading/%x_%j.err)
  if [ -n "$1" ]; then opts+=(--dependency="afterany:$1"); fi
  local id; id=$(submit "${SB[@]}" "${opts[@]}" scripts/simulation/grade.sh "${args[@]}")
  if [ "$PRINT" != 1 ]; then
    echo "submitted the summaries as job $id${1:+ (after jobs $1)}; output in outputs/simulation/grading/"
  fi
}

case "$STAGE" in
  list)
    T="${1:?list prep|production}"; echo "$T cells ($(cells_of "$T")):"; list_cells "$T" ;;

  prep|production)
    parse_opts "$@"
    [ "$PRINT" = 1 ] || need_setup
    T=$(cells_of "$STAGE")
    if [ "$STAGE" = production ]; then install_kappa
    elif [ "$FRESH" = 1 ]; then
      [ ! -f "$KAPPA" ] || { echo "--fresh-kappa: $KAPPA exists; move it aside to derive kappa_AE afresh"; exit 1; }
      echo "kappa_AE will be derived afresh by prep-panel (written to $KAPPA)"
    else install_kappa; fi
    mapfile -t ALL < <(awk -F'\t' 'NR > 1 && NF >= 4 {print $1}' "$T")
    if [ -n "$ONLY" ]; then IFS=',' read -r -a SEL <<< "$ONLY"; else SEL=("${ALL[@]}"); fi
    for c in "${SEL[@]}"; do printf '%s\n' "${ALL[@]}" | grep -qx -- "$c" || { echo "unknown cell '$c' in $T"; exit 2; }; done
    mkdir -p "outputs/simulation/logs/$STAGE"
    LEDGER="outputs/simulation/logs/$STAGE/jobs.tsv"
    BATCH=$(date -u +%Y%m%dT%H%M%SZ)
    if [ "$PRINT" != 1 ] && [ ! -f "$LEDGER" ]; then printf 'batch\tsubmitted_utc\tcell\tjob\tlimit_min\ttasks\n' > "$LEDGER"; fi
    declare -A ID_OF=()
    # cells without a dependency first; prep-summary ("after" = prep) waits for every other preparation job
    ORDERED=(); LAST=()
    for c in "${SEL[@]}"; do
      if [ "$(col "$T" "$c" 9 -)" = prep ]; then LAST+=("$c"); else ORDERED+=("$c"); fi
    done
    ORDERED+=(${LAST[@]+"${LAST[@]}"})
    for c in "${ORDERED[@]}"; do
      m=$(limit_of "$STAGE" "$c"); mem=$(col "$T" "$c" 6 8G); t=$(col "$T" "$c" 7 1); after=$(col "$T" "$c" 9 -)
      opts=(--job-name="sim_${STAGE}_$c" --time="$m" --mem="$mem" --cpus-per-task=1 --ntasks=1 --nodes=1
            --chdir="$ROOT" --export="ALL,SIM_ROOT=$ROOT,SIM_STAGE=$STAGE,SIM_CELL=$c,SIM_CELLS=$ROOT/$T")
      if [ "$t" -gt 1 ]; then
        opts+=(--array="1-$t" --output="outputs/simulation/logs/$STAGE/%x_%A_%a.out" --error="outputs/simulation/logs/$STAGE/%x_%A_%a.err")
      else
        opts+=(--output="outputs/simulation/logs/$STAGE/%x_%j.out" --error="outputs/simulation/logs/$STAGE/%x_%j.err")
      fi
      if [ "$after" = prep ]; then
        d=""; for k in "${!ID_OF[@]}"; do d="$d:${ID_OF[$k]}"; done
        if [ -n "$d" ]; then opts+=(--dependency="afterok$d" --kill-on-invalid-dep=yes); fi
      fi
      id=$(submit "${SB[@]}" "${opts[@]}" scripts/simulation/job.sh)
      ID_OF[$c]=$id
      if [ "$PRINT" != 1 ]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$BATCH" "$(date -u +%FT%TZ)" "$c" "$id" "$m" "$t" >> "$LEDGER"
        echo "  submitted $c as job $id$([ "$t" -gt 1 ] && echo " [array 1-$t]")"
      fi
    done
    [ "$PRINT" = 1 ] || echo "batch $BATCH: ${#ORDERED[@]} job(s) recorded in $LEDGER; outputs in outputs/simulation/$STAGE/"
    if [ "$THEN_GRADE" = 1 ]; then
      [ "$STAGE" = production ] || { echo "--then-grade applies to production only"; exit 2; }
      ids=$(printf '%s\n' "${ID_OF[@]}" | paste -sd: -)
      submit_grade "$ids"
    fi ;;

  grade)
    parse_opts "$@"
    [ "$PRINT" = 1 ] || need_setup
    if [ "$AFTER_PROD" = 1 ]; then
      L=outputs/simulation/logs/production/jobs.tsv
      [ -f "$L" ] || { echo "--after-production: no production batch recorded in $L"; exit 1; }
      last=$(awk -F'\t' 'NR > 1 {b = $1} END {print b}' "$L")
      AFTER=$(awk -F'\t' -v b="$last" 'NR > 1 && $1 == b && $4 ~ /^[0-9]+$/ {print $4}' "$L" | paste -sd: -)
      [ -n "$AFTER" ] || { echo "--after-production: production batch $last in $L has no job ids"; exit 1; }
      echo "the summaries wait for every job of production batch $last: $(echo "$AFTER" | tr ':' '\n' | grep -c .) job(s)"
    fi
    submit_grade "$AFTER" ;;

  local)
    S="${1:?local prep|production <cell> [<task>]}"; C="${2:?cell}"; K="${3:-}"
    need_setup; T=$(cells_of "$S"); install_kappa
    n=$(col "$T" "$C" 7 1); [ -n "$n" ] || { echo "unknown cell '$C' in $T"; exit 2; }
    if [ "$n" -gt 1 ] && [ -z "$K" ]; then echo "$C has $n array tasks: give a task number 1..$n"; exit 2; fi
    SIM_ROOT="$ROOT" SIM_STAGE="$S" SIM_CELL="$C" SIM_CELLS="$ROOT/$T" SIM_TASK="$K" bash scripts/simulation/job.sh ;;

  *) sed -n '2,16p' "$0"; exit 2 ;;
esac
