#!/usr/bin/env bash
# Full FilterDDP runs: baseline lu(K) against UMFPACK's symmetric strategy
# (FILTERDDP_UMFPACK_STRATEGY=3), the only configuration with a clear isolated
# gain in kkt_ordering_benchmark.csv (ieee123 T=3, exact Hessian: 0.79x on
# factor + wide solve). Matched protocol, full per-iteration logging, runs
# alternated baseline/candidate so machine drift hits both arms.
#
#   bash ddp/examples/power_system/run_umfpack_strategy_fullrun.sh [system] [T] [repeats]

set -u
cd "$(dirname "$0")/../../.." || exit 1
SYS="${1:-ieee123C_1ph}"; T="${2:-3}"; REP="${3:-2}"
OUT=ddp/results/kkt_ordering/fullrun
mkdir -p "$OUT" ddp/results/kkt_ordering/captures/solutions
REF=$(grep -oE "CENTRAL_IPOPT .*" "ddp/results/matched_ipopt_race/logs/ipopt_${SYS}_T${T}.log" | \
      grep -oE " objective=[-0-9.eE+]+" | cut -d= -f2)
[ -n "$REF" ] || { echo "no matched Ipopt objective for $SYS T=$T"; exit 1; }
case $SYS in ieee2522C_1ph) PRIMAL=1e-5;; large10kC_1ph) PRIMAL=1e-4;; *) PRIMAL=1e-6;; esac

export REDUCED_PROFILE=periodic REDUCED_CB=1e-3 TERMINAL_SOC_SOFT=1
export FILTERDDP_MAX_ITERATIONS=400
export FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
export FILTERDDP_NEAR_OPT_REFERENCE="$REF" FILTERDDP_NEAR_OPT_GAP=0.005 FILTERDDP_NEAR_OPT_PRIMAL="$PRIMAL"
JL="julia --startup-file=no"
SOL=ddp/results/network_filterddp/filterddp_solution_${SYS}_T${T}_periodic_CB1e-3.jls

for ARM in exact diag; do
  for r in $(seq 1 "$REP"); do
    for VARIANT in baseline umf_sym; do
      TAG="${SYS}_T${T}_${ARM}_${VARIANT}_r${r}"
      LOG="$OUT/fddp_${TAG}.log"
      (
        if [ "$ARM" = diag ]; then export FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8
        else unset FILTERDDP_DIAG_HESSIAN FILTERDDP_DIAG_HESSIAN_FLOOR; fi
        if [ "$VARIANT" = umf_sym ]; then export FILTERDDP_UMFPACK_STRATEGY=3
        else unset FILTERDDP_UMFPACK_STRATEGY; fi
        echo "PIPELINE_ENV system=$SYS T=$T arm=$ARM variant=$VARIANT repeat=$r umfpack_strategy=${FILTERDDP_UMFPACK_STRATEGY:-default} near_opt_reference=$REF near_opt_primal=$PRIMAL started=$(date '+%Y-%m-%dT%H:%M:%S')"
        $JL --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl "$SYS" "$T" solve
      ) > "$LOG" 2>&1
      cp "$SOL" "ddp/results/kkt_ordering/captures/solutions/sol_${TAG}.jls" 2>/dev/null
      $JL --project=envs/ddp2026 ddp/examples/power_system/extract_filterddp_feasibility_trace.jl \
          "$LOG" "$OUT/trace_${TAG}.csv" "$REF" >> "$LOG" 2>&1
      echo "$TAG: $(grep -E 'solve complete' "$LOG" | tail -1)  $(grep -oE 'FilterDDP objective=[-0-9.eE+]+' "$LOG")  $(grep -E 'NEAR_OPT' "$LOG" | head -1 | cut -c1-120)"
    done
  done
done
