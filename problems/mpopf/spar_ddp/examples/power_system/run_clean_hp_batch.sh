#!/usr/bin/env bash
# Clean batch on one footing (user's go of 2026-10-10): both solvers on the
# same screened model and power-flow start, under the Windows power plan that
# is active when it runs (recorded; High performance when this was written),
# with every iteration logged so that any stopping rule can be applied to
# both afterwards.
#
#   Ipopt    MA57, MA97, MUMPS; tol = 1e-6; IPOPT_TRACE=1 (clock, objective,
#            primal and dual infeasibility per iteration; ipopt_trace_posthoc.py)
#   FilterDDP  tolerance 1e-5, no near-optimality stop; diagonal Hessian, and
#            exact battery-block Hessian with the lean step B
#            (FILTERDDP_LEAN_EXACT=1); lean path, halved block, tree solver;
#            full warm-up; one-worker-per-period accounting (parsim_from_log.py),
#            both primal infeasibilities (near_opt_posthoc.py [--primal=original])
#
# One case at a time, each after a quiet wait, background load logged.
# Resumable: a case whose log is complete is skipped.
#
#   bash problems/mpopf/spar_ddp/examples/power_system/run_clean_hp_batch.sh [ipopt|ddp|all]
set -u
cd "$(dirname "$0")/../../../../.." || exit 1
unset NO_COLOR FORCE_COLOR
PART=${1:-all}
TAG=${CLEAN_TAG:-hp_tol1e5}            # CLEAN_TAG, CLEAN_ARMS, CLEAN_CELLS: a later batch in another configuration
OUT=problems/mpopf/spar_ddp/results/clean_hp; mkdir -p "$OUT"
JL="julia --startup-file=no"
LOADJL=problems/mpopf/spar_ddp/examples/power_system/sample_background_load.jl
SMALL="ieee123C_1ph:6 ieee123C_1ph:24 ieee123C_1ph:96 ieee2522C_1ph:6 ieee2522C_1ph:24 ieee2522C_1ph:96"
LARGE="large10kC_1ph:6 large10kC_1ph:24 large10kC_1ph:48 large10kC_1ph:96"
echo "CLEAN_HP_BATCH part=$PART start=$(date '+%Y-%m-%dT%H:%M:%S') plan=$(powercfg //getactivescheme | grep -oE '\(.*\)')" | tee -a "$OUT/batch.log"

if [ "$PART" = ipopt ] || [ "$PART" = all ]; then
  export IPOPT_SCREEN=substation,vupper,ell,psubs IPOPT_LOADFLOW_START=1 IPOPT_TRACE=1 IPOPT_MORE_OPTIONS="tol=1e-6"
  for LS in ma57 ma97 mumps; do
    CELLS="$SMALL" SOLVERS=$LS IPOPT_TAG_SUFFIX=_screen_lf_hp bash problems/mpopf/spar_ddp/examples/power_system/run_ipopt_hsl.sh
    # MUMPS at large10k T=96 has not been run before and may not fit in memory: last, below.
    [ $LS = mumps ] && C="${LARGE% *}" || C="$LARGE"
    CELLS="$C" SOLVERS=$LS IPOPT_TAG_SUFFIX=_screen_lf_pv_hp bash problems/mpopf/spar_ddp/examples/power_system/run_ipopt_hsl.sh
  done 2>&1 | tee -a "$OUT/batch.log"
  unset IPOPT_SCREEN IPOPT_LOADFLOW_START IPOPT_TRACE IPOPT_MORE_OPTIONS
fi

if [ "$PART" = ddp ] || [ "$PART" = all ]; then
  export FILTERDDP_LEAN_VALUE=1 FILTERDDP_HALF_BLOCK=1 STRICT=1 FILTERDDP_OPTIMALITY_TOLERANCE=1e-5
  unset FILTERDDP_PINT_SCAN
  for ARM in ${CLEAN_ARMS:-diag exact}; do
    case $ARM in exact) export FILTERDDP_LEAN_EXACT=1;; *) unset FILTERDDP_LEAN_EXACT;; esac
    for cell in ${CLEAN_CELLS:-$SMALL $LARGE}; do
      S=${cell%%:*}; T=${cell##*:}
      LOG=problems/mpopf/spar_ddp/results/parallel_in_time/logs/fddp_${S}_T${T}_${ARM}_${TAG}.log
      [ -s "$LOG" ] && grep -q "solve complete" "$LOG" && { echo "skip $S T=$T $ARM"; continue; }
      case $S in large10kC_1ph) REF=problems/mpopf/spar_ddp/results/ipopt_hsl/logs/ipopt_ma57_screen_lf_pv_hp_${S}_T$T.log;;
                 *)             REF=problems/mpopf/spar_ddp/results/ipopt_hsl/logs/ipopt_ma57_screen_lf_hp_${S}_T$T.log;; esac
      QW=$($JL "$LOADJL" wait 10 1.5 2>/dev/null | grep QUIET_WAIT)
      STOP="$OUT/.sampling_$$"; touch "$STOP"
      $JL "$LOADJL" sample "$STOP" "$OUT/load_fddp_${S}_T${T}_${ARM}_${TAG}.csv" 15 > /dev/null 2>&1 &
      SPID=$!
      IPOPT_REF_LOG="$REF" bash problems/mpopf/spar_ddp/examples/power_system/run_parsim.sh "$S" "$T" "$ARM" "$TAG" | head -1
      rm -f "$STOP"; wait "$SPID" 2>/dev/null
      echo "  $QW"
    done
  done 2>&1 | tee -a "$OUT/batch.log"
fi

if [ "$PART" = all ]; then
  export IPOPT_SCREEN=substation,vupper,ell,psubs IPOPT_LOADFLOW_START=1 IPOPT_TRACE=1 IPOPT_MORE_OPTIONS="tol=1e-6"
  CELLS="large10kC_1ph:96" SOLVERS=mumps IPOPT_TAG_SUFFIX=_screen_lf_pv_hp bash problems/mpopf/spar_ddp/examples/power_system/run_ipopt_hsl.sh 2>&1 | tee -a "$OUT/batch.log"
fi
echo "CLEAN_HP_BATCH_DONE part=$PART end=$(date '+%Y-%m-%dT%H:%M:%S')" | tee -a "$OUT/batch.log"
