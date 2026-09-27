#!/usr/bin/env bash
# Capture one representative stage-1 KKT system per (system, Hessian arm) for the
# MA57 / ordering benchmark (ddp/notes/KKT_ORDERING_AND_MA57.md).
#
# Matched protocol: periodic profile, C_B = 1e-3, soft terminal SOC, T = 3.
# FILTERDDP_CAPTURE_KKT rewrites its file at every backward pass of stage
# FILTERDDP_CAPTURE_STAGE, so stopping at FILTERDDP_MAX_ITERATIONS=$ITER leaves the
# stage-1 K and its full (n_x+1)-column RHS from iteration $ITER -- mid-solve,
# with the barrier parameter and active set away from the cold start.
#
# Captures are large (~1.6 GB at large10k) and gitignored.
#
#   bash ddp/examples/power_system/capture_ordering_kkt.sh [iteration]

set -u
cd "$(dirname "$0")/../../.." || exit 1
ITER="${1:-20}"
OUT=ddp/results/kkt_ordering
mkdir -p "$OUT/captures" "$OUT/logs"

export REDUCED_PROFILE=periodic REDUCED_CB=1e-3 TERMINAL_SOC_SOFT=1
export FILTERDDP_SKIP_SOLUTION_WRITE=1 FILTERDDP_CAPTURE_STAGE=1
export FILTERDDP_MAX_ITERATIONS="$ITER"
export FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1

for SYS in ieee123C_1ph ieee2522C_1ph large10kC_1ph; do
  for ARM in diag exact; do
    CAP="$OUT/captures/kkt_${SYS}_T3_${ARM}_iter${ITER}_stage1.jls"
    LOG="$OUT/logs/capture_${SYS}_T3_${ARM}.log"
    if [ -s "$CAP" ] && grep -q "CAPTURE_DONE" "$LOG" 2>/dev/null; then
      echo "skip $SYS $ARM"; continue
    fi
    (
      if [ "$ARM" = diag ]; then
        export FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8
      else
        unset FILTERDDP_DIAG_HESSIAN FILTERDDP_DIAG_HESSIAN_FLOOR
      fi
      export FILTERDDP_CAPTURE_KKT="$CAP"
      echo "CAPTURE_ENV system=$SYS T=3 arm=$ARM iteration=$ITER stage=1 C_B=$REDUCED_CB profile=$REDUCED_PROFILE started=$(date '+%Y-%m-%dT%H:%M:%S')"
      julia --startup-file=no --project=envs/ddp2026 \
        ddp/examples/power_system/ieee123c_filterddp.jl "$SYS" 3 solve
      echo "CAPTURE_DONE exit=$? finished=$(date '+%Y-%m-%dT%H:%M:%S')"
    ) > "$LOG" 2>&1
    echo "$SYS $ARM: $(stat -c %s "$CAP" 2>/dev/null) bytes"
  done
done
