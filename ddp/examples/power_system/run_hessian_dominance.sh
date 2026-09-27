#!/usr/bin/env bash
# How diagonal is the battery-power block of the exact stage Hessian, and how
# much of its diagonal is C_B? (hessian_dominance_analysis.jl.) One untimed
# EXACT-Hessian FilterDDP solve per (system, C_B) at T=6 in the paper's
# formulation, capturing every stage each STRIDE-th iteration and analysing the
# captures as they are written (then deleting them). C_B = 0, the per-system
# value ("system") and the earlier uniform 1e-3; large10k is sampled every 10th
# iteration and run at 0 and its per-system value only, for time.
#
#   bash ddp/examples/power_system/run_hessian_dominance.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
OUT=ddp/results/hessian_dominance
CAPROOT=ddp/results/kkt_ordering/captures/hessian_dominance
mkdir -p "$OUT" "$CAPROOT"
JL="julia --startup-file=no"
export REDUCED_PROFILE=periodic TERMINAL_SOC_SOFT=1
export FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
export FILTERDDP_MAX_ITERATIONS="${MAXIT:-150}" FILTERDDP_SKIP_SOLUTION_WRITE=1
unset FILTERDDP_DIAG_HESSIAN FILTERDDP_DIRECT_DIAG_HESSIAN FILTERDDP_FACTOR_BACKED_POLICY \
      FILTERDDP_BLOCKED_SOLVE FILTERDDP_NEAR_OPT_REFERENCE

for cell in ieee123C_1ph:0:1 ieee123C_1ph:1e-3:1 \
            ieee2522C_1ph:0:1 ieee2522C_1ph:system:1 ieee2522C_1ph:1e-3:1 \
            large10kC_1ph:0:10 large10kC_1ph:system:10; do
  IFS=: read -r S CB STRIDE <<< "$cell"; T=6
  TAG="${S}_T${T}_cb${CB}"
  CSV="$OUT/dominance_${TAG}.csv"
  [ -s "$CSV" ] && { echo "skip $TAG"; continue; }
  CAP="$CAPROOT/$TAG"; rm -rf "$CAP"; mkdir -p "$CAP"; DONE="$CAP/.done"
  KKT_EVOL_DELETE=1 OPENBLAS_NUM_THREADS=1 $JL --project=envs/ddp2026 \
    ddp/examples/power_system/hessian_dominance_analysis.jl "$CAP" \
    "ddp/results/network_filterddp/network_data_${S}_T${T}_periodic.jls" "$S" "$CB" "$CSV" "$DONE" \
    > "$OUT/analysis_${TAG}.log" 2>&1 &
  APID=$!
  echo "PIPELINE_ENV system=$S T=$T cb=$CB hessian=exact stride=$STRIDE started=$(date '+%Y-%m-%dT%H:%M:%S')" > "$OUT/fddp_${TAG}.log"
  REDUCED_CB="$CB" FILTERDDP_PERIODIC_CAPTURE_DIR="$CAP" FILTERDDP_PERIODIC_CAPTURE_STRIDE="$STRIDE" \
    $JL --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl "$S" "$T" solve \
    >> "$OUT/fddp_${TAG}.log" 2>&1
  touch "$DONE"; wait "$APID"
  echo "[$(date '+%H:%M:%S')] $TAG: $(grep -E 'solve complete' "$OUT/fddp_${TAG}.log" | tail -1) | $(tail -1 "$OUT/analysis_${TAG}.log")"
done
echo "DOMINANCE_DONE"
