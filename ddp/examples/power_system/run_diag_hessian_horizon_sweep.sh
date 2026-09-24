#!/usr/bin/env bash
# Exact vs diagonal stage Hessian at horizons beyond T = 3.
#
# The T = 3 result (FILTERDDP_DIAGONAL_HESSIAN.md): same optimum, +16-20%
# iterations, wall -10% on ieee123/ieee2522 and -46 to -57% on large10k. The
# open question is whether the iteration penalty grows with the horizon, since
# the deleted block carries V_xx and V_xx matters more the longer the horizon.
#
#   bash ddp/examples/power_system/run_diag_hessian_horizon_sweep.sh [outdir]
#
# Settings are IDENTICAL to the T = 3 pipeline arms so the tables join: periodic
# price profile, C_B = 1e-3, floor 1e-8, timing and nnz diagnostics on. Two
# deliberate differences, both applied to BOTH arms of a pair:
#   * FILTERDDP_MAX_ITERATIONS=400 (was the 200 default) -- longer horizons may
#     need more, and a cap hit would waste the run. T = 3 never came near 200.
#   * large10k T >= 12 sets FILTERDDP_FACTOR_BACKED_POLICY=1. Without it the
#     dense per-stage policy (beta, omega) is ~9.5 GB at T = 12 on a 32 GB
#     machine with ~13 GB free. The setting is written into every log.
#
# Cases run cheapest-first, so partial results are useful whenever it stops.
# RESUMABLE: a finished log is skipped; delete it to redo. After every case the
# summary and CSV are regenerated and committed+pushed, so results reach the
# remote even if the Claude app or the remote connection drops mid-sweep.

set -u
cd "$(dirname "$0")/../../.." || exit 1
REPO="$PWD"
OUT="${1:-$REPO/ddp/results/diag_hessian_horizon}"
JL="julia --startup-file=no"
DRIVER="ddp/examples/power_system/ieee123c_filterddp.jl"
DATA="ddp/results/network_filterddp"
BRANCH="$(git rev-parse --abbrev-ref HEAD)"
mkdir -p "$OUT/logs"

say() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
done_marker() { [ -s "$1" ] && grep -q "PIPELINE_STEP_OK" "$1"; }

# ---------------------------------------------------------------- stay awake --
LOCK="$OUT/.sweep_running.lock"
touch "$LOCK"
$JL ddp/examples/power_system/hold_awake.jl "$LOCK" 16 > "$OUT/logs/hold_awake.log" 2>&1 &
trap 'rm -f "$LOCK"' EXIT

CASES=(
  "ieee123C_1ph 6"   "ieee123C_1ph 12"  "ieee123C_1ph 24"
  "ieee2522C_1ph 6"  "ieee123C_1ph 48"  "ieee2522C_1ph 12"
  "ieee123C_1ph 96"  "ieee2522C_1ph 24"
  "large10kC_1ph 6"  "large10kC_1ph 12"
)

export REDUCED_PROFILE=periodic
export REDUCED_CB=1e-3
export FILTERDDP_SKIP_SOLUTION_WRITE=1
export FILTERDDP_TIMING_DIAGNOSTIC=1
export FILTERDDP_NNZ_DIAGNOSTIC=1
export FILTERDDP_MAX_ITERATIONS=400

summarise_and_push() {
  $JL --project=envs/ddp2026 ddp/examples/power_system/summarize_diag_hessian_sweep.jl \
      "$OUT" "$REPO/ddp/results/agenda_pipeline/logs" > "$OUT/SUMMARY.txt" 2>&1
  git add "$OUT/SUMMARY.txt" "$OUT/diag_hessian_horizon.csv" "$OUT/sweep_run.log" 2>/dev/null
  if ! git diff --cached --quiet; then
    git commit -q -m "diag-Hessian horizon sweep: $1 (auto-commit by run_diag_hessian_horizon_sweep.sh)

Unreviewed pipeline output: SUMMARY.txt and diag_hessian_horizon.csv regenerated
from the raw logs after this case. See the script header for settings.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
    if git pull -q --rebase origin "$BRANCH" && git push -q origin "$BRANCH"; then
      say "pushed after $1"
    else
      git rebase --abort 2>/dev/null
      say "push failed after $1 -- committed locally, will retry after the next case"
    fi
  fi
}

for C in "${CASES[@]}"; do
  set -- $C; SYS=$1; T=$2
  FILE="$DATA/network_data_${SYS}_T${T}_periodic.jls"

  if [ ! -s "$FILE" ]; then
    say "export: $SYS T=$T (periodic profile)"
    PROFILE_PERIODIC=1 $JL --project=envs/tadmm ddp/examples/power_system/export_ieee123c_data.jl \
        "$SYS" "$T" > "$OUT/logs/export_${SYS}_T${T}.log" 2>&1
    [ -s "$FILE" ] || { say "export FAILED for $SYS T=$T, skipping case"; continue; }
  fi

  FB=0
  [ "$SYS" = large10kC_1ph ] && [ "$T" -ge 12 ] && FB=1

  for ARM in exact diag; do
    LOG="$OUT/logs/${SYS}_T${T}_${ARM}.log"
    done_marker "$LOG" && { say "skip $SYS T=$T $ARM"; continue; }
    say "run: $SYS T=$T $ARM (factor_backed=$FB)"
    (
      if [ "$ARM" = diag ]; then
        export FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8
      else
        unset FILTERDDP_DIAG_HESSIAN FILTERDDP_DIAG_HESSIAN_FLOOR
      fi
      if [ "$FB" = 1 ]; then export FILTERDDP_FACTOR_BACKED_POLICY=1
      else unset FILTERDDP_FACTOR_BACKED_POLICY; fi
      echo "PIPELINE_ENV system=$SYS T=$T arm=$ARM factor_backed=$FB diag_floor=${FILTERDDP_DIAG_HESSIAN_FLOOR:-none} max_iterations=$FILTERDDP_MAX_ITERATIONS C_B=$REDUCED_CB profile=$REDUCED_PROFILE started=$(date '+%Y-%m-%dT%H:%M:%S')"
      $JL --project=envs/ddp2026 "$DRIVER" "$SYS" "$T" solve
    ) > "$LOG" 2>&1
    echo "PIPELINE_STEP_OK" >> "$LOG"
  done

  summarise_and_push "$SYS T=$T"
done

say "sweep finished"
summarise_and_push "sweep finished"
