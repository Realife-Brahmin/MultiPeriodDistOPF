#!/usr/bin/env bash
# KKT matrix evolution between FilterDDP iterations (kkt_evolution_analysis.jl).
# One untimed FilterDDP run per system at T=6, in the paper's formulation
# (periodic profile, per-system C_B, soft terminal SOC, diagonal Hessian and the
# exact assembly rewrites), capturing every stage's KKT system each
# STRIDE-th iteration through DDP4OPF's periodic capture hook. The factor-backed
# policy and the blocked solve are left off: they change how the KKT system is
# solved, not the matrix, and the capture hook needs the explicit policy.
# No near-optimality stop, so each run follows FilterDDP to its own tolerance
# (at most MAXIT iterations).
#
# The analysis follows each capture directory while FilterDDP writes it and
# deletes med2522/large10k captures once processed (a large10k stage snapshot
# is ~2.4 GB). ieee123's are kept. Solutions are not written, so no stored
# solution file is replaced.
#
#   bash ddp/examples/power_system/run_kkt_evolution.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
OUT=ddp/results/kkt_evolution
CAPROOT=ddp/results/kkt_ordering/captures/kkt_evolution
mkdir -p "$OUT" "$CAPROOT"
JL="julia --startup-file=no"
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
export FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
export FILTERDDP_MAX_ITERATIONS="${MAXIT:-200}" FILTERDDP_SKIP_SOLUTION_WRITE=1
unset FILTERDDP_FACTOR_BACKED_POLICY FILTERDDP_BLOCKED_SOLVE FILTERDDP_NEAR_OPT_REFERENCE

for cell in ieee123C_1ph:6:1:0 ieee2522C_1ph:6:1:1 large10kC_1ph:6:10:1; do
  IFS=: read -r S T STRIDE DEL <<< "$cell"
  CSV="$OUT/kkt_evolution_${S}_T${T}.csv"
  [ -s "$CSV" ] && { echo "skip $S T=$T"; continue; }
  CAP="$CAPROOT/${S}_T${T}"; rm -rf "$CAP"; mkdir -p "$CAP"; DONE="$CAP/.done"
  KKT_EVOL_DELETE=$DEL OPENBLAS_NUM_THREADS=1 $JL --project=envs/ddp2026 \
    ddp/examples/power_system/kkt_evolution_analysis.jl "$CAP" \
    "ddp/results/network_filterddp/network_data_${S}_T${T}_periodic.jls" "$S" "$CSV" "$DONE" \
    > "$OUT/analysis_${S}_T${T}.log" 2>&1 &
  APID=$!
  echo "PIPELINE_ENV system=$S T=$T stride=$STRIDE cb=system started=$(date '+%Y-%m-%dT%H:%M:%S')" > "$OUT/fddp_${S}_T${T}.log"
  FILTERDDP_PERIODIC_CAPTURE_DIR="$CAP" FILTERDDP_PERIODIC_CAPTURE_STRIDE="$STRIDE" \
    $JL --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl "$S" "$T" solve \
    >> "$OUT/fddp_${S}_T${T}.log" 2>&1
  touch "$DONE"; wait "$APID"
  echo "[$(date '+%H:%M:%S')] $S T=$T: $(grep -E 'solve complete' "$OUT/fddp_${S}_T${T}.log" | tail -1) | $(tail -1 "$OUT/analysis_${S}_T${T}.log")"
  [ -s "$CSV" ] || { echo "KKT_EVOLUTION_FAILED at $S"; exit 1; }
done
echo "KKT_EVOLUTION_DONE"
