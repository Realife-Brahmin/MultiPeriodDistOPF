#!/usr/bin/env bash
# large10k on the corrected instance. Until 2026-10-05 the parser read every
# large10k PV unit as 0 kW (see parse_opendss.jl), so every earlier large10k
# result, for either solver, is on a system without its 357.7 MW of PV. The
# instances were re-exported (old ones kept as *_periodic_nopv.jls); every log
# written here carries the tag `pv`, so nothing collides with or is skipped
# because of an earlier run.
#
#   now     both solvers in the current configuration (screened model,
#           power-flow start; FilterDDP also the exact step), T = 6, 24, 48:
#           Ipopt with MUMPS, MA57, MA97; FilterDDP to near-optimality and to
#           its own tolerance, compilation excluded;
#   before  the earlier configuration on the same corrected instance, for a
#           like-for-like before/after: Ipopt unscreened from the flat start,
#           FilterDDP with the tree solver, flat start and halving line search.
#
# One run at a time, quiet machine, background load logged, full per-iteration
# FilterDDP logs. Completed logs are skipped. PARTS="now before" selects.
#
#   bash ddp/examples/power_system/run_large10k_pv_rerun.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
until [ "$(tasklist //FI "IMAGENAME eq julia.exe" 2>/dev/null | grep -c julia.exe)" = "0" ]; do sleep 30; done
echo "[$(date '+%H:%M:%S')] large10k pv rerun start"
PARTS="${PARTS:-now before}"
S=ddp/examples/power_system/run_blocked_solve_fullrun.sh
H=ddp/results/ipopt_hsl/logs
SCREEN=substation,vupper,ell,psubs
CELLS="large10kC_1ph:6 large10kC_1ph:24 large10kC_1ph:48"
export PIN_BLAS_THREADS=0 VARIANTS=blocked_w16 CB=system JULIA_NUM_THREADS=8
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_TYPED_EQUATIONS=1 FILTERDDP_STRUCTURED_DYNAMICS=1
export FILTERDDP_TREE_KKT=1

for part in $PARTS; do case $part in
now)
  IPOPT_SCREEN=$SCREEN IPOPT_LOADFLOW_START=1 IPOPT_TAG_SUFFIX=_screen_lf_pv SOLVERS="ma57 ma97 mumps" CELLS="$CELLS" \
    bash ddp/examples/power_system/run_ipopt_hsl.sh
  export FILTERDDP_SCREEN=$SCREEN FILTERDDP_LOADFLOW_START=1 FILTERDDP_AFFINE_LINESEARCH=1
  for T in 6 24 48; do
    REF=$H/ipopt_ma57_screen_lf_pv_large10kC_1ph_T$T.log
    FILTERDDP_WARMUP=full RUN_TAG_SUFFIX=_tableII_jt8_cbsys_typed_tree3_screen_warm_pv \
      IPOPT_REF_LOG=$REF bash "$S" large10kC_1ph $T diag 16 1
  done
  for T in 6 24 48; do
    REF=$H/ipopt_ma57_screen_lf_pv_large10kC_1ph_T$T.log
    [ $T = 6 ] && W=full || W=8
    STRICT=1 FILTERDDP_WARMUP=$W RUN_TAG_SUFFIX=_tableII_jt8_cbsys_typed_tree3_screen_warm_strict_pv \
      IPOPT_REF_LOG=$REF bash "$S" large10kC_1ph $T diag 16 1
  done
  unset FILTERDDP_SCREEN FILTERDDP_LOADFLOW_START FILTERDDP_AFFINE_LINESEARCH
  ;;
before)
  IPOPT_TAG_SUFFIX=_pv SOLVERS="ma57 ma97 mumps" CELLS="$CELLS" bash ddp/examples/power_system/run_ipopt_hsl.sh
  for T in 6 24 48; do
    REF=$H/ipopt_ma57_screen_lf_pv_large10kC_1ph_T$T.log
    FILTERDDP_WARMUP=full RUN_TAG_SUFFIX=_tableII_jt8_cbsys_typed_tree3_warm_pv \
      IPOPT_REF_LOG=$REF bash "$S" large10kC_1ph $T diag 16 1
  done
  ;;
esac; done
echo "[$(date '+%H:%M:%S')] LARGE10K_PV_DONE"
