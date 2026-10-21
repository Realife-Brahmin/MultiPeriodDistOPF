#!/usr/bin/env bash
# OpenDSS check of the dispatches behind the clean batch
# (run_clean_hp_batch.sh). The timed runs do not write their solutions, so each
# case is solved once more here, untimed, in the same configuration and to
# the same tolerance (1e-5), with the solution kept, and that dispatch is handed
# to OpenDSS (run_opendss_check.sh). Run it only when no timed run is going.
#
#   bash spar_ddp/examples/power_system/run_clean_hp_opendss.sh
#   CELLS="ieee123C_1ph:6" ARMS="exact" bash spar_ddp/examples/power_system/run_clean_hp_opendss.sh
#
# A near-optimality stop can be checked instead of the converged point by
# setting NEAR_OPT_GAP (e.g. 0.005): the run then stops at the first iterate
# within that gap of the matched Ipopt objective and under the system's primal
# threshold, and the tag carries the gap.
set -u
cd "$(dirname "$0")/../../.." || exit 1
unset NO_COLOR FORCE_COLOR
CELLS=${CELLS:-"ieee123C_1ph:6 ieee123C_1ph:24 ieee123C_1ph:96 ieee2522C_1ph:6 ieee2522C_1ph:24 ieee2522C_1ph:96 large10kC_1ph:6 large10kC_1ph:24 large10kC_1ph:48 large10kC_1ph:96"}
ARMS=${ARMS:-"diag exact"}
GAP=${NEAR_OPT_GAP:-}
SOLDIR=spar_ddp/results/kkt_ordering/captures/solutions; mkdir -p "$SOLDIR"
OUT=spar_ddp/results/clean_hp; mkdir -p "$OUT"
H=spar_ddp/results/ipopt_hsl/logs
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
export FILTERDDP_MAX_ITERATIONS=400 FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
export FILTERDDP_OPTIMALITY_TOLERANCE=1e-5
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_TYPED_EQUATIONS=1 FILTERDDP_STRUCTURED_DYNAMICS=1
export FILTERDDP_TREE_KKT=1 FILTERDDP_BLOCKED_SOLVE=16
export FILTERDDP_SCREEN=substation,vupper,ell,psubs FILTERDDP_LOADFLOW_START=1 FILTERDDP_AFFINE_LINESEARCH=1
export FILTERDDP_LEAN_VALUE=1 FILTERDDP_HALF_BLOCK=1
unset FILTERDDP_SKIP_SOLUTION_WRITE FILTERDDP_WARMUP FILTERDDP_PARSIM FILTERDDP_PINT_SCAN
for ARM in $ARMS; do
  if [ "$ARM" = diag ]; then export FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8; unset FILTERDDP_LEAN_EXACT
  else unset FILTERDDP_DIAG_HESSIAN FILTERDDP_DIAG_HESSIAN_FLOOR; export FILTERDDP_LEAN_EXACT=1; fi
  for cell in $CELLS; do
    S=${cell%%:*}; T=${cell##*:}
    case $S in
      ieee123C_1ph)  PRIMAL=1e-6; PV="";    REF=$H/ipopt_ma57_screen_lf_hp_${S}_T$T.log;    export JULIA_NUM_THREADS=1;;
      ieee2522C_1ph) PRIMAL=1e-5; PV="";    REF=$H/ipopt_ma57_screen_lf_hp_${S}_T$T.log;    export JULIA_NUM_THREADS=8;;
      large10kC_1ph) PRIMAL=1e-5; PV="_pv"; REF=$H/ipopt_ma57_screen_lf_pv_hp_${S}_T$T.log; export JULIA_NUM_THREADS=8;;
    esac
    unset FILTERDDP_NEAR_OPT_REFERENCE FILTERDDP_NEAR_OPT_GAP FILTERDDP_NEAR_OPT_PRIMAL
    STOP=tol1e5
    if [ -n "$GAP" ]; then
      OBJ=$(grep -oE "CENTRAL_IPOPT .*" "$REF" | grep -oE " objective=[-0-9.eE+]+" | cut -d= -f2)
      [ -n "$OBJ" ] || { echo "no Ipopt reference in $REF"; continue; }
      export FILTERDDP_NEAR_OPT_REFERENCE="$OBJ" FILTERDDP_NEAR_OPT_GAP="$GAP" FILTERDDP_NEAR_OPT_PRIMAL="$PRIMAL"
      STOP=gap$GAP
    fi
    RUNTAG="${S}_T${T}_${ARM}_hp_${STOP}${PV}_opendss"
    [ -s "spar_ddp/results/opendss_check/$RUNTAG.log" ] && grep -q OPENDSS_SUMMARY "spar_ddp/results/opendss_check/$RUNTAG.log" && { echo "skip $RUNTAG"; continue; }
    LOG="$OUT/fddp_${RUNTAG}.log"
    MARK="$OUT/.mark_$$"; touch "$MARK"
    julia --startup-file=no --project=envs/spar_ddp spar_ddp/examples/power_system/ieee123c_filterddp.jl "$S" "$T" solve > "$LOG" 2>&1
    SOL=spar_ddp/results/network_filterddp/filterddp_solution_${S}_T${T}_periodic_CBsystem.jls
    # the solution must be the one this run wrote, not an older file of the same name
    grep -q "solve complete" "$LOG" && [ "$SOL" -nt "$MARK" ] || { echo "solve failed or no new solution: $RUNTAG"; rm -f "$MARK"; continue; }
    rm -f "$MARK"
    cp "$SOL" "$SOLDIR/sol_${RUNTAG}.jls"
    echo "[$(date '+%H:%M:%S')] $RUNTAG: $(grep -E 'solve complete' "$LOG" | tail -1) | $(grep -oE 'BOUND_CHECK.*' "$LOG" | cut -c1-110)"
    bash spar_ddp/examples/power_system/run_opendss_check.sh "$S" "$T" "$SOLDIR/sol_${RUNTAG}.jls" "$RUNTAG" | grep -E "voltage:|OPENDSS_SUMMARY|FAILED" | cut -c1-330
  done
done
echo "CLEAN_HP_OPENDSS_DONE"
