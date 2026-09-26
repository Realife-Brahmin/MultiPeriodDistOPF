#!/usr/bin/env bash
# Full FilterDDP runs: FilterDDP's ldiv!(F, rhs) against the blocked
# multi-right-hand-side solve over the same UMFPACK factors
# (FILTERDDP_BLOCKED_SOLVE=<w>). Matched protocol, full per-iteration logging,
# baseline and blocked alternated, background load logged, solutions kept for
# comparison.
#
#   bash ddp/examples/power_system/run_blocked_solve_fullrun.sh <system> <T> <arm> <w> [repeats]
#
# Other FILTERDDP_* switches set by the caller (e.g. the exact assembly rewrites
# FILTERDDP_DIRECT_DIAG_HESSIAN, FILTERDDP_TRIPLET_SECOND_DERIVATIVES,
# FILTERDDP_CACHE_KKT_PATTERN) apply to both arms; give such runs a distinct
# RUN_TAG_SUFFIX so their logs do not collide with plain ones.

set -u
cd "$(dirname "$0")/../../.." || exit 1
SYS=$1; T=$2; ARM=$3; W=$4; REP="${5:-1}"
SUFFIX="${RUN_TAG_SUFFIX:-}"
VARIANTS="${VARIANTS:-baseline blocked_w$W}"      # VARIANTS=blocked_w16 runs the blocked arm only
REWRITES="direct_diag=${FILTERDDP_DIRECT_DIAG_HESSIAN:-0} triplet=${FILTERDDP_TRIPLET_SECOND_DERIVATIVES:-0} kkt_pattern_cache=${FILTERDDP_CACHE_KKT_PATTERN:-0} factor_backed=${FILTERDDP_FACTOR_BACKED_POLICY:-0}"
OUT=ddp/results/kkt_ordering/fullrun_blocked
SOLDIR=ddp/results/kkt_ordering/captures/solutions
mkdir -p "$OUT" "$SOLDIR"
REF=$(grep -oE "CENTRAL_IPOPT .*" "ddp/results/matched_ipopt_race/logs/ipopt_${SYS}_T${T}.log" | \
      grep -oE " objective=[-0-9.eE+]+" | cut -d= -f2)
[ -n "$REF" ] || { echo "no matched Ipopt objective for $SYS T=$T"; exit 1; }
case $SYS in ieee2522C_1ph) PRIMAL=1e-5;; large10kC_1ph) PRIMAL=1e-4;; *) PRIMAL=1e-6;; esac

export REDUCED_PROFILE=periodic REDUCED_CB=1e-3 TERMINAL_SOC_SOFT=1
export FILTERDDP_MAX_ITERATIONS=400
export FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
export FILTERDDP_NEAR_OPT_REFERENCE="$REF" FILTERDDP_NEAR_OPT_GAP=0.005 FILTERDDP_NEAR_OPT_PRIMAL="$PRIMAL"
if [ "$ARM" = diag ]; then export FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8; fi
# BLAS/OpenMP pinned to one thread by default. PIN_BLAS_THREADS=0 leaves them
# unset, which is how the matched race (Table II of the paper) ran; the thread
# count Julia actually uses is recorded either way.
if [ "${PIN_BLAS_THREADS:-1}" = 1 ]; then export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
else unset OMP_NUM_THREADS OPENBLAS_NUM_THREADS; fi
JL="julia --startup-file=no"
BLAS_THREADS=$($JL -e 'using LinearAlgebra; print(BLAS.get_num_threads())')
LOADJL=ddp/examples/power_system/sample_background_load.jl
SOL=ddp/results/network_filterddp/filterddp_solution_${SYS}_T${T}_periodic_CB1e-3.jls

for r in $(seq 1 "$REP"); do
  for VARIANT in $VARIANTS; do
    TAG="${SYS}_T${T}_${ARM}${SUFFIX}_${VARIANT}_r${r}"
    LOG="$OUT/fddp_${TAG}.log"
    [ -s "$LOG" ] && grep -q "solve complete" "$LOG" && { echo "skip $TAG"; continue; }
    QW=$($JL "$LOADJL" wait 10 1.5 2>/dev/null | grep QUIET_WAIT)
    STOP="$OUT/.sampling_$$"; touch "$STOP"
    $JL "$LOADJL" sample "$STOP" "$OUT/load_${TAG}.csv" 15 > /dev/null 2>&1 &
    SPID=$!
    (
      if [ "$VARIANT" = baseline ]; then unset FILTERDDP_BLOCKED_SOLVE
      else export FILTERDDP_BLOCKED_SOLVE="$W"; fi
      echo "$QW"
      echo "PIPELINE_ENV system=$SYS T=$T arm=$ARM variant=$VARIANT repeat=$r blocked_solve=${FILTERDDP_BLOCKED_SOLVE:-off} $REWRITES blas_threads=$BLAS_THREADS julia_threads=${JULIA_NUM_THREADS:-1} near_opt_reference=$REF near_opt_primal=$PRIMAL started=$(date '+%Y-%m-%dT%H:%M:%S')"
      $JL --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl "$SYS" "$T" solve
    ) > "$LOG" 2>&1
    rm -f "$STOP"; wait "$SPID" 2>/dev/null
    cp "$SOL" "$SOLDIR/sol_${TAG}.jls" 2>/dev/null
    $JL --project=envs/ddp2026 ddp/examples/power_system/extract_filterddp_feasibility_trace.jl \
        "$LOG" "$OUT/trace_${TAG}.csv" "$REF" >> "$LOG" 2>&1
    echo "[$(date '+%H:%M:%S')] $TAG: $(grep -E 'solve complete' "$LOG" | tail -1) | $(grep -oE 'FILTERDDP_NEAR_OPT iteration=[0-9]+ elapsed_s=[0-9.]+' "$LOG" | head -1) | $(grep -oE 'FilterDDP objective=[-0-9.eE+]+' "$LOG")"
  done
done
