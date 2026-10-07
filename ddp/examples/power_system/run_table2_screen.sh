#!/usr/bin/env bash
# The Table II cells with the three changes of ddp/notes/CONSTRAINT_SCREENING.md,
# for both solvers on the same model and from the same starting point:
#   - exact bound screening (substation voltage, upper voltage limits the
#     rules prove cannot bind, ell >= 0, P_Subs >= 0),
#   - the power-flow start (batteries idle),
#   - FilterDDP only: the exact-step line search for affine dynamics.
#
# FilterDDP: each system in its fastest diagonal-Hessian configuration
# (ieee123: structured dynamics + UMFPACK on one thread; med2522 and large10k:
# the tree solver on eight), once as before (tag _screen, compilation inside
# the timed solve) and once as the second solve in the same process
# (FILTERDDP_WARMUP=full, tag _screen_warm: compilation excluded).
# Ipopt: MUMPS, MA57, MA97 on 1 and 8 threads (run_ipopt_hsl.sh), logs tagged
# _screen_lf in ddp/results/ipopt_hsl/logs/.
#
# Cells whose log is complete are skipped. Starts only when no julia.exe is
# running; every case waits for a quiet machine and logs the background load.
#
#   bash ddp/examples/power_system/run_table2_screen.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
until [ "$(tasklist //FI "IMAGENAME eq julia.exe" 2>/dev/null | grep -c julia.exe)" = "0" ]; do sleep 30; done
echo "[$(date '+%H:%M:%S')] queue start"
SCREEN=substation,vupper,ell,psubs
CELLS="${CELLS:-ieee123C_1ph:6 ieee123C_1ph:24 ieee123C_1ph:96 ieee2522C_1ph:6 ieee2522C_1ph:24 ieee2522C_1ph:96 large10kC_1ph:6 large10kC_1ph:24 large10kC_1ph:48}"

if [ "${SKIP_IPOPT:-0}" != 1 ]; then
  IPOPT_SCREEN=$SCREEN IPOPT_LOADFLOW_START=1 IPOPT_TAG_SUFFIX=_screen_lf CELLS="$CELLS" \
    bash ddp/examples/power_system/run_ipopt_hsl.sh
fi

S=ddp/examples/power_system/run_blocked_solve_fullrun.sh
export PIN_BLAS_THREADS=0 VARIANTS=blocked_w16 CB=system
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_TYPED_EQUATIONS=1 FILTERDDP_STRUCTURED_DYNAMICS=1
export FILTERDDP_SCREEN=$SCREEN FILTERDDP_LOADFLOW_START=1 FILTERDDP_AFFINE_LINESEARCH=1
for warm in 0 full; do
  [ $warm = 0 ] && { unset FILTERDDP_WARMUP; W=""; } || { export FILTERDDP_WARMUP=$warm; W=_warm; }
  for cell in $CELLS; do
    sys=${cell%%:*}; T=${cell##*:}
    case $sys in
      ieee123C_1ph) unset FILTERDDP_TREE_KKT; export JULIA_NUM_THREADS=1 RUN_TAG_SUFFIX=_tableII_cbsys_typed_sd_screen$W;;
      *)            export FILTERDDP_TREE_KKT=1 JULIA_NUM_THREADS=8 RUN_TAG_SUFFIX=_tableII_jt8_cbsys_typed_tree3_screen$W;;
    esac
    bash "$S" "$sys" "$T" diag 16 1
  done
done
echo "[$(date '+%H:%M:%S')] TABLE2_SCREEN_DONE"
