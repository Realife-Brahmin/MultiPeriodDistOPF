#!/usr/bin/env bash
# tp_test.sh <system> <T> <diag|exact> <threads> [tag]   (untimed correctness run)
set -u
R="/c/Users/Aryan Ritwajeet Jha/Documents/documents_general/MultiPeriodDistOPF/.claude/worktrees/new-session-a7b7a2"
cd "$R" || exit 1
unset NO_COLOR FORCE_COLOR
S=$1; T=$2; ARM=$3; TH=$4; TAG=${5:-tp}
OUT="$(dirname "$0")/logs"; mkdir -p "$OUT"
export TP_REPO="$R"
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
export FILTERDDP_MAX_ITERATIONS=400 FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
export FILTERDDP_OPTIMALITY_TOLERANCE=1e-5
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_TYPED_EQUATIONS=1 FILTERDDP_STRUCTURED_DYNAMICS=1
export FILTERDDP_TREE_KKT=1 FILTERDDP_BLOCKED_SOLVE=16 FILTERDDP_SKIP_SOLUTION_WRITE=1
export FILTERDDP_SCREEN=substation,vupper,ell,psubs FILTERDDP_LOADFLOW_START=1 FILTERDDP_AFFINE_LINESEARCH=1
export FILTERDDP_LEAN_VALUE=1 FILTERDDP_HALF_BLOCK=1
unset FILTERDDP_WARMUP FILTERDDP_PARSIM FILTERDDP_PINT_SCAN FILTERDDP_NEAR_OPT_REFERENCE FILTERDDP_DIAG_HESSIAN FILTERDDP_DIAG_HESSIAN_FLOOR FILTERDDP_LEAN_EXACT
case $ARM in
  diag)  export FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8;;
  exact) export FILTERDDP_LEAN_EXACT=1;;
esac
export JULIA_NUM_THREADS=$TH
LOG="$OUT/tp_${S}_T${T}_${ARM}_${TAG}.log"
julia --startup-file=no --project=envs/ddp2026 "$(dirname "$0")/threaded_test.jl" "$S" "$T" solve > "$LOG" 2>&1
echo "$S T=$T $ARM threads=$TH $TAG: $(grep -E 'solve complete' "$LOG" | tail -1) | $(grep -oE 'FilterDDP objective=[-0-9.eE+]+' "$LOG") | threaded passes $(grep -c FILTERDDP_THREADED_PASS "$LOG"), reuse passes $(grep -c FILTERDDP_REUSE_PASS "$LOG")"
grep -E "ERROR|Error|error:" "$LOG" | head -5
python ddp/examples/power_system/compare_runs.py ddp/results/parallel_in_time/logs/fddp_${S}_T${T}_${ARM}_hp_tol1e5.log "$LOG"
