#!/usr/bin/env bash
# What the screening rules other than the substation one are worth
# (ddp/notes/CONSTRAINT_SCREENING.md, Section 10). Untimed: iterations to
# near-optimality only, no warm-up. Four runs per system: the substation rule
# alone or all four rules, each with and without the power-flow start and the
# exact step. Otherwise the configuration of run_parsim.sh.
#
#   bash ddp/examples/power_system/run_screen_rule_ablation.sh
set -u
cd "$(dirname "$0")/../../.." || exit 1
unset NO_COLOR FORCE_COLOR
OUT=ddp/results/screen_ablation; mkdir -p "$OUT"
H=ddp/results/ipopt_hsl/logs
run() {  # system T variant screen loadflow affine
  local SYS=$1 T=$2 V=$3 SCR=$4 LF=$5 AFF=$6 PRIMAL REF THREADS
  case $SYS in
    ieee123C_1ph)  PRIMAL=1e-6; REF=$H/ipopt_ma57_screen_lf_${SYS}_T$T.log; THREADS=1;;
    ieee2522C_1ph) PRIMAL=1e-5; REF=$H/ipopt_ma57_screen_lf_${SYS}_T$T.log; THREADS=8;;
    large10kC_1ph) PRIMAL=1e-4; REF=$H/ipopt_ma57_screen_lf_pv_${SYS}_T$T.log; THREADS=8;;
  esac
  local OBJ; OBJ=$(grep -oE "CENTRAL_IPOPT .*" "$REF" | grep -oE " objective=[-0-9.eE+]+" | cut -d= -f2)
  [ -n "$OBJ" ] || { echo "no Ipopt reference in $REF"; return; }
  local LOG="$OUT/rules_${SYS}_T${T}_${V}.log"
  (
    export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
    export FILTERDDP_MAX_ITERATIONS=400 FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
    export FILTERDDP_NEAR_OPT_REFERENCE="$OBJ" FILTERDDP_NEAR_OPT_GAP=0.005 FILTERDDP_NEAR_OPT_PRIMAL="$PRIMAL"
    export FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8
    export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
    export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_TYPED_EQUATIONS=1 FILTERDDP_STRUCTURED_DYNAMICS=1
    export FILTERDDP_TREE_KKT=1 FILTERDDP_BLOCKED_SOLVE=16 FILTERDDP_SKIP_SOLUTION_WRITE=1
    export FILTERDDP_SCREEN="$SCR"
    unset FILTERDDP_LOADFLOW_START FILTERDDP_AFFINE_LINESEARCH FILTERDDP_WARMUP FILTERDDP_PARSIM
    [ "$LF" = 1 ] && export FILTERDDP_LOADFLOW_START=1
    [ "$AFF" = 1 ] && export FILTERDDP_AFFINE_LINESEARCH=1
    export JULIA_NUM_THREADS=$THREADS
    julia --startup-file=no --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl "$SYS" "$T" solve > "$LOG" 2>&1
  )
  echo "[$(date '+%H:%M:%S')] $SYS T=$T $V: $(grep -oE 'FILTERDDP_NEAR_OPT iteration=[0-9]+' "$LOG" | head -1) | $(grep -oE 'BOUND_CHECK.*' "$LOG" | cut -c1-110)"
}
ALL=substation,vupper,ell,psubs
for cell in ieee123C_1ph:24 ieee2522C_1ph:24 large10kC_1ph:6; do
  S=${cell%%:*}; T=${cell##*:}
  run "$S" "$T" full_PE   "$ALL"     1 1     # current configuration (control)
  run "$S" "$T" sub_PE    substation 1 1     # substation rule only, with start and step
  run "$S" "$T" sub_only  substation 0 0     # substation rule only
  run "$S" "$T" full_only "$ALL"     0 0     # all four rules, nothing else
done
echo SCREEN_RULE_ABLATION_DONE
