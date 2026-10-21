#!/usr/bin/env bash
# The battery value recursion by pairwise combination of the periods
# (pint_scan.jl), checked and timed inside an ordinary FilterDDP run with the
# exact battery-block Hessian. Current configuration otherwise (run_parsim.sh):
# screened model, power-flow start, exact step, tree solver, lean path, halved
# block. Untimed as a solve: no warm-up, and the scanned passes carry the
# diagnostic's own work.
#
#   bash spar_ddp/examples/power_system/run_pint_scan.sh <system> <T> [stride] [offset] [max_iterations] [tag]
#
# The scan runs on the backward passes k with k % stride == offset. The run
# stops at near-optimality or after max_iterations, whichever comes first.
# Logs go to spar_ddp/results/pint_scan/logs/.
set -u
cd "$(dirname "$0")/../../.." || exit 1
unset NO_COLOR FORCE_COLOR
SYS=$1; T=$2; STRIDE=${3:-4}; OFFSET=${4:-2}; MAXIT=${5:-400}; TAG=${6:-scan}
OUT=spar_ddp/results/pint_scan/logs; mkdir -p "$OUT"
H=spar_ddp/results/ipopt_hsl/logs
case $SYS in
  ieee123C_1ph)  PRIMAL=1e-6; REF=$H/ipopt_ma57_screen_lf_${SYS}_T$T.log; THREADS=1;;
  ieee2522C_1ph) PRIMAL=1e-5; REF=$H/ipopt_ma57_screen_lf_${SYS}_T$T.log; THREADS=8;;
  large10kC_1ph) PRIMAL=1e-5; REF=$H/ipopt_ma57_screen_lf_pv_${SYS}_T$T.log; THREADS=8;;
esac
OBJ=$(grep -oE "CENTRAL_IPOPT .*" "$REF" | grep -oE " objective=[-0-9.eE+]+" | cut -d= -f2)
[ -n "$OBJ" ] || { echo "no Ipopt reference in $REF"; exit 1; }
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
export FILTERDDP_MAX_ITERATIONS=$MAXIT FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
export FILTERDDP_NEAR_OPT_REFERENCE="$OBJ" FILTERDDP_NEAR_OPT_GAP=0.005 FILTERDDP_NEAR_OPT_PRIMAL="$PRIMAL"
unset FILTERDDP_DIAG_HESSIAN FILTERDDP_DIAG_HESSIAN_FLOOR FILTERDDP_WARMUP
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_TYPED_EQUATIONS=1 FILTERDDP_STRUCTURED_DYNAMICS=1
export FILTERDDP_TREE_KKT=1 FILTERDDP_BLOCKED_SOLVE=16 FILTERDDP_SKIP_SOLUTION_WRITE=1
export FILTERDDP_SCREEN=substation,vupper,ell,psubs FILTERDDP_LOADFLOW_START=1 FILTERDDP_AFFINE_LINESEARCH=1
export FILTERDDP_LEAN_VALUE=1 FILTERDDP_HALF_BLOCK=1 FILTERDDP_PARSIM=1
export FILTERDDP_PINT_SCAN=1 FILTERDDP_PINT_SCAN_STRIDE=$STRIDE FILTERDDP_PINT_SCAN_OFFSET=$OFFSET
export JULIA_NUM_THREADS=${JULIA_NUM_THREADS:-$THREADS}
LOG="$OUT/pint_${SYS}_T${T}_${TAG}.log"
julia --startup-file=no --project=envs/spar_ddp spar_ddp/examples/power_system/ieee123c_filterddp.jl "$SYS" "$T" solve > "$LOG" 2>&1
echo "[$(date '+%H:%M:%S')] $SYS T=$T $TAG: $(grep -E 'solve complete' "$LOG" | tail -1) | $(grep -oE 'FILTERDDP_NEAR_OPT iteration=[0-9]+' "$LOG" | head -1) | scans=$(grep -c '^PINT_SCAN ' "$LOG")"
grep -E "^PINT_SCAN|ERROR|rror:" "$LOG" | cut -c1-2000 | head -40
