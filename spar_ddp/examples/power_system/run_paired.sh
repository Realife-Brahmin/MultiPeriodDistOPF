#!/usr/bin/env bash
# One FilterDDP solve with the pairwise sweep inside it (paired_pass.jl), exact
# Hessian, tolerance 1e-5, factors reused after a barrier update, otherwise the
# configuration of run_parsim.sh. The dispatch goes to OpenDSS afterwards.
#
#   bash spar_ddp/examples/power_system/run_paired.sh <system> <T> <count|1|procs> [tag]
#
# count: sequential, every piece timed, for the one-worker count (clean_hp_table.py
#        --arms=exact,paired). Warm-up of PAIRED_WARMUP iterations (default 3).
# 1:     on Julia threads (JULIA_NUM_THREADS, default 8); correctness, not timing.
# procs: every period in a process of its own (procs_pass.jl, procs_main.jl), on
#        this machine: FILTERDDP_PROCS_WORKERS workers (default one per period)
#        and a coordinator. A full pass after a barrier update, not the
#        factor-reusing one. Warm-up of PAIRED_WARMUP iterations (default full).
# PAIRED_CHECK=1 also runs the ordinary sweep in every pass and prints the difference.
set -u
cd "$(dirname "$0")/../../.." || exit 1
unset NO_COLOR FORCE_COLOR
SYS=$1; T=$2; MODE=$3; TAG=${4:-${MODE}_hp_tol1e5}
[ "$MODE" = 1 ] && TAG=${4:-threads_hp_tol1e5}
[ "$MODE" = procs ] && TAG=${4:-procs_hp_tol1e5}
OUT=spar_ddp/results/parallel_in_time/logs
case $SYS in
  large10kC_1ph) REF=spar_ddp/results/ipopt_hsl/logs/ipopt_ma57_screen_lf_pv_hp_${SYS}_T$T.log;;
  *)             REF=spar_ddp/results/ipopt_hsl/logs/ipopt_ma57_screen_lf_hp_${SYS}_T$T.log;;
esac
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
export FILTERDDP_MAX_ITERATIONS=400 FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
# PAIRED_TOL: another tolerance (ieee123 T=24 ends above its 1e-6 infeasibility threshold at 1e-5).
export FILTERDDP_OPTIMALITY_TOLERANCE=${PAIRED_TOL:-1e-5}
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_TYPED_EQUATIONS=1 FILTERDDP_STRUCTURED_DYNAMICS=1
export FILTERDDP_TREE_KKT=1 FILTERDDP_BLOCKED_SOLVE=16
export FILTERDDP_SCREEN=substation,vupper,ell,psubs FILTERDDP_LOADFLOW_START=1 FILTERDDP_AFFINE_LINESEARCH=1
export FILTERDDP_LEAN_VALUE=1 FILTERDDP_HALF_BLOCK=1 FILTERDDP_LEAN_EXACT=1 FILTERDDP_REUSE_BARRIER=1
unset FILTERDDP_DIAG_HESSIAN FILTERDDP_DIAG_HESSIAN_FLOOR FILTERDDP_SKIP_SOLUTION_WRITE FILTERDDP_NEAR_OPT_REFERENCE FILTERDDP_THREADED FILTERDDP_PINT_SCAN
export FILTERDDP_PAIRED_SWEEP=$MODE FILTERDDP_PAIRED_CHECK=${PAIRED_CHECK:-0}
# PAIRED_INPLACE=0: the stages' functions and the combinations allocated afresh every
# pass, as before 2026-10-11 (same solve; checked on ieee123 T=6, 24 and med2522 T=6).
# PAIRED_MAIN: another script in place of the wrapper, e.g. alloc_profile.jl.
export FILTERDDP_PAIRED_INPLACE=${PAIRED_INPLACE:-1}
[ -n "${PAIRED_MAX_ITERATIONS:-}" ] && export FILTERDDP_MAX_ITERATIONS=$PAIRED_MAX_ITERATIONS
# PAIRED_ROUNDS=0: the ordinary sweep in the same pass (its count on the same accounting).
export FILTERDDP_PAIRED_ROUNDS=${PAIRED_ROUNDS:-1}
if [ "$MODE" = count ]; then
  export FILTERDDP_PARSIM=1 FILTERDDP_PARSIM_GCNET=1 FILTERDDP_WARMUP=${PAIRED_WARMUP:-3}
  case $SYS in ieee123C_1ph) export JULIA_NUM_THREADS=1;; *) export JULIA_NUM_THREADS=${JULIA_NUM_THREADS:-8};; esac
elif [ "$MODE" = procs ]; then
  unset FILTERDDP_PARSIM FILTERDDP_PARSIM_GCNET FILTERDDP_REUSE_BARRIER
  export FILTERDDP_WARMUP=${PAIRED_WARMUP:-full} JULIA_NUM_THREADS=1
  PAIRED_MAIN=${PAIRED_MAIN:-spar_ddp/examples/power_system/procs_main.jl}
else
  unset FILTERDDP_PARSIM; export FILTERDDP_WARMUP=0 JULIA_NUM_THREADS=${JULIA_NUM_THREADS:-8}
fi
LOG="$OUT/fddp_${SYS}_T${T}_paired_${TAG}.log"
MARK="$OUT/.mark_$$"; touch "$MARK"
julia --startup-file=no --project=envs/spar_ddp "${PAIRED_MAIN:-spar_ddp/examples/power_system/paired_main.jl}" "$SYS" "$T" solve > "$LOG" 2>&1
SOL=spar_ddp/results/network_filterspar_ddp/filterddp_solution_${SYS}_T${T}_periodic_CBsystem.jls
SOLDIR=spar_ddp/results/kkt_ordering/captures/solutions; mkdir -p "$SOLDIR"
RUNTAG="parsim_${SYS}_T${T}_paired_${TAG}"
if grep -q "solve complete" "$LOG" && [ "$SOL" -nt "$MARK" ]; then
  cp "$SOL" "$SOLDIR/sol_${RUNTAG}.jls"
  bash spar_ddp/examples/power_system/run_opendss_check.sh "$SYS" "$T" "$SOLDIR/sol_${RUNTAG}.jls" "$RUNTAG" 2>&1 | grep -E "OPENDSS_SUMMARY.*voltage|FAILED" | cut -c1-230
else
  echo "OPENDSS_CHECK: no solution written by this run ($RUNTAG)"
fi
rm -f "$MARK"
echo "[$(date '+%H:%M:%S')] $SYS T=$T paired $TAG: $(grep -E 'solve complete' "$LOG" | tail -1) | $(grep -oE 'FilterDDP objective=[-0-9.eE+]+' "$LOG") | paired passes $(grep -c FILTERDDP_PAIRED_PASS "$LOG"), procs passes $(grep -c FILTERDDP_PROCS_PASS "$LOG"), reuse passes $(grep -c FILTERDDP_REUSE_PASS "$LOG"), worst scan error $(grep -oE 'scan_err_Vxx=[0-9.e+-]+' "$LOG" | cut -d= -f2 | sort -g | tail -1)"
grep -E "ERROR|WARMUP_STREAM" "$LOG" | head -3
