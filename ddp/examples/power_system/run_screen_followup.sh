#!/usr/bin/env bash
# Follow-up to run_table2_screen.sh: every earlier iteration-dependent result
# was measured with the line search throttled by the substation voltage bound
# (ddp/notes/CONSTRAINT_SCREENING.md), so the cases that fed a conclusion are
# rerun in the new configuration (screening + power-flow start + exact step),
# compilation excluded by a warm-up (a full first solve on the Table II
# cells, 8 iterations at long horizons). One run at a time, quiet machine, full
# per-iteration logs. Completed logs are skipped.
#
#   A. longer horizons, diagonal Hessian, tree solver: med2522 T=192, 384 and
#      large10k T=96, 192, with Ipopt (MA57, MA97) on the same model and start;
#   B. exact Hessian through the tree solver on the Table II cells;
#   C. C_B = 0 with the diagonal Hessian, where it failed before
#      (large10k T=6, med2522 T=96);
#   D. med2522 T=1152, and T=1536 (failed before; no Ipopt reference, strict).
#
# Waits for the file given as $1 to contain TABLE2_SCREEN_DONE (or for no
# julia.exe if no argument). PARTS="A B" selects parts.
#
#   bash ddp/examples/power_system/run_screen_followup.sh [wait_file]

set -u
cd "$(dirname "$0")/../../.." || exit 1
if [ -n "${1:-}" ]; then until grep -q TABLE2_SCREEN_DONE "$1" 2>/dev/null; do sleep 30; done; fi
until [ "$(tasklist //FI "IMAGENAME eq julia.exe" 2>/dev/null | grep -c julia.exe)" = "0" ]; do sleep 30; done
echo "[$(date '+%H:%M:%S')] follow-up start"
PARTS="${PARTS:-A B C D}"
S=ddp/examples/power_system/run_blocked_solve_fullrun.sh
H=ddp/results/ipopt_hsl/logs
SCREEN=substation,vupper,ell,psubs
export PIN_BLAS_THREADS=0 VARIANTS=blocked_w16 CB=system
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_TYPED_EQUATIONS=1 FILTERDDP_STRUCTURED_DYNAMICS=1
export FILTERDDP_SCREEN=$SCREEN FILTERDDP_LOADFLOW_START=1 FILTERDDP_AFFINE_LINESEARCH=1
export FILTERDDP_TREE_KKT=1 JULIA_NUM_THREADS=8

for part in $PARTS; do case $part in
A)
  IPOPT_SCREEN=$SCREEN IPOPT_LOADFLOW_START=1 IPOPT_TAG_SUFFIX=_screen_lf SOLVERS="ma57 ma97" \
    CELLS="ieee2522C_1ph:192 large10kC_1ph:96 ieee2522C_1ph:384 large10kC_1ph:192" \
    bash ddp/examples/power_system/run_ipopt_hsl.sh
  export FILTERDDP_WARMUP=8 RUN_TAG_SUFFIX=_jt8_cbsys_typed_tree3_screen_warm
  IPOPT_REF_LOG=$H/ipopt_ma57_ieee2522C_1ph_T192.log      bash "$S" ieee2522C_1ph 192 diag 16 1
  IPOPT_REF_LOG=$H/ipopt_ma97_large10kC_1ph_T96.log       bash "$S" large10kC_1ph 96 diag 16 1
  IPOPT_REF_LOG=$H/ipopt_oom_ma57_ieee2522C_1ph_T384.log  bash "$S" ieee2522C_1ph 384 diag 16 1
  IPOPT_REF_LOG=$H/ipopt_ma97_large10kC_1ph_T192.log      bash "$S" large10kC_1ph 192 diag 16 1
  ;;
B)
  export FILTERDDP_WARMUP=full
  for cell in ieee123C_1ph:6 ieee123C_1ph:24 ieee123C_1ph:96 ieee2522C_1ph:6 ieee2522C_1ph:24 \
              ieee2522C_1ph:96 large10kC_1ph:6 large10kC_1ph:24; do
    sys=${cell%%:*}; T=${cell##*:}
    case $sys in ieee123C_1ph) export JULIA_NUM_THREADS=1 RUN_TAG_SUFFIX=_tableII_cbsys_typed_tree4_screen_warm;;
                 *)            export JULIA_NUM_THREADS=8 RUN_TAG_SUFFIX=_tableII_jt8_cbsys_typed_tree4_screen_warm;; esac
    bash "$S" "$sys" "$T" exact 16 1
  done
  export JULIA_NUM_THREADS=8
  ;;
C)
  export FILTERDDP_WARMUP=full RUN_TAG_SUFFIX=_tableII_jt8_cb0_typed_tree3_screen_warm
  CB=0 bash "$S" large10kC_1ph 6 diag 16 1
  CB=0 bash "$S" ieee2522C_1ph 96 diag 16 1
  ;;
D)
  export FILTERDDP_WARMUP=8 RUN_TAG_SUFFIX=_jt8_cbsys_typed_tree3_screen_warm
  IPOPT_REF_LOG=$H/ipopt_oom_ma57_ieee2522C_1ph_T1152.log bash "$S" ieee2522C_1ph 1152 diag 16 1
  LOG=ddp/results/kkt_ordering/fullrun_blocked/fddp_ieee2522C_1ph_T1536_diag_jt8_cbsys_typed_tree3_screen_warm_strict_r1.log
  if ! grep -q "solve complete" "$LOG" 2>/dev/null; then
    REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1 FILTERDDP_MAX_ITERATIONS=400 \
    FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1 FILTERDDP_DIAG_HESSIAN=1 \
    FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8 FILTERDDP_BLOCKED_SOLVE=16 FILTERDDP_SKIP_SOLUTION_WRITE=1 \
      julia --startup-file=no --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl \
      ieee2522C_1ph 1536 solve > "$LOG" 2>&1
    echo "[$(date '+%H:%M:%S')] ieee2522C_1ph_T1536 strict: $(grep -E 'solve complete' "$LOG") | $(grep -oE 'FilterDDP objective=\S+' "$LOG") | $(grep -oE 'maxrss_mib=\S+' "$LOG")"
  fi
  ;;
esac; done
echo "[$(date '+%H:%M:%S')] SCREEN_FOLLOWUP_DONE"
