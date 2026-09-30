#!/usr/bin/env bash
# Stage-1 KKT systems for the entry-magnitude study (kkt_magnitude_analysis.jl),
# in the paper's current formulation: T = 6, periodic profile, per-system C_B,
# soft terminal SOC, diagonal Hessian with the exact assembly rewrites.
# FILTERDDP_CAPTURE_KKT rewrites its file at every backward pass of stage 1, so
# stopping at FILTERDDP_MAX_ITERATIONS = k keeps iteration k's system. Three
# points per system: early (5), mid-solve (40) and near the Table II
# near-optimality iteration (ieee123 67, med2522 68, large10k 100).
# Captures are gitignored (large10k ~1.6 GB each).
#
#   bash ddp/examples/power_system/capture_kkt_magnitude.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
CAPDIR=ddp/results/kkt_ordering/captures/magnitude
LOGDIR=ddp/results/kkt_magnitude/logs_capture
mkdir -p "$CAPDIR" "$LOGDIR"

export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
export FILTERDDP_SKIP_SOLUTION_WRITE=1 FILTERDDP_CAPTURE_STAGE=1
export FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
export FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_BLOCKED_SOLVE=16 JULIA_NUM_THREADS=8
unset FILTERDDP_FACTOR_BACKED_POLICY FILTERDDP_NEAR_OPT_REFERENCE

for cell in ieee123C_1ph:5 ieee123C_1ph:40 ieee123C_1ph:67 \
            ieee2522C_1ph:5 ieee2522C_1ph:40 ieee2522C_1ph:68 \
            large10kC_1ph:5 large10kC_1ph:40 large10kC_1ph:100; do
  S=${cell%%:*}; K=${cell##*:}
  CAP="$CAPDIR/kkt_${S}_T6_diag_cbsys_iter${K}_stage1.jls"
  LOG="$LOGDIR/capture_${S}_T6_iter${K}.log"
  [ -s "$CAP" ] && grep -q "solve complete" "$LOG" 2>/dev/null && { echo "skip $S iter $K"; continue; }
  FILTERDDP_CAPTURE_KKT="$CAP" FILTERDDP_MAX_ITERATIONS="$K" \
    julia --startup-file=no --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl "$S" 6 solve \
    > "$LOG" 2>&1
  echo "[$(date '+%H:%M:%S')] $S iter $K: $(grep -E 'solve complete' "$LOG" | tail -1) | $(ls -la "$CAP" 2>/dev/null | awk '{print $5}') bytes"
done
echo "MAGNITUDE_CAPTURE_DONE"
