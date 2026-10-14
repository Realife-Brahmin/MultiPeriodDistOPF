#!/usr/bin/env bash
# One FilterDDP run in the current configuration (screened model, power-flow
# start, exact step, tree solver) with the one-worker-per-period accounting
# (FILTERDDP_PARSIM=1, see backward_pass.jl and parsim_from_log.py). The run is
# the ordinary sequential one; only its timers differ.
#
#   bash ddp/examples/power_system/run_parsim.sh <system> <T> <diag|exact> [tag]
#
# PARSIM=0 runs the same configuration without the accounting (regression).
# STRICT=1 runs to the solver's own tolerance instead of near-optimality.
# Logs go to ddp/results/parallel_in_time/logs/.
set -u
cd "$(dirname "$0")/../../.." || exit 1
# NO_COLOR (or FORCE_COLOR) makes Julia wrap stdout in an IOContext. The warm-up
# solve prints to a plain IOStream, so every print statement would then be
# compiled again in the first pass of the timed solve: about 0.7 s, and 1.5-2.3 s
# with the FILTERDDP_PARSIM lines. Found 2026-10-06 in runs launched from a
# PowerShell-started queue, whose environment carries NO_COLOR.
unset NO_COLOR FORCE_COLOR
SYS=$1; T=$2; ARM=$3; TAG=${4:-parsim}
OUT=ddp/results/parallel_in_time/logs; mkdir -p "$OUT"
H=ddp/results/ipopt_hsl/logs
case $SYS in
  ieee123C_1ph)  PRIMAL=1e-6; REF=$H/ipopt_ma57_screen_lf_${SYS}_T$T.log; THREADS=1;;
  ieee2522C_1ph) PRIMAL=1e-5; REF=$H/ipopt_ma57_screen_lf_${SYS}_T$T.log; THREADS=8;;
  large10kC_1ph) PRIMAL=1e-5; REF=$H/ipopt_ma57_screen_lf_pv_${SYS}_T$T.log; THREADS=8;;
esac
REF=${IPOPT_REF_LOG:-$REF}
OBJ=$(grep -oE "CENTRAL_IPOPT .*" "$REF" | grep -oE " objective=[-0-9.eE+]+" | cut -d= -f2)
[ -n "$OBJ" ] || { echo "no Ipopt reference in $REF"; exit 1; }
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
export FILTERDDP_MAX_ITERATIONS=400 FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
[ "${STRICT:-0}" = 1 ] || export FILTERDDP_NEAR_OPT_REFERENCE="$OBJ" FILTERDDP_NEAR_OPT_GAP=0.005 FILTERDDP_NEAR_OPT_PRIMAL="$PRIMAL"
# Arms: diag (diagonal Hessian); exact (the caller sets FILTERDDP_LEAN_EXACT=1
# for the lean step); hybrid (diagonal until the barrier parameter has fallen
# to HYBRID_MU, exact with the lean step from then on).
case $ARM in
  diag)   export FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8;;
  hybrid) export FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8 FILTERDDP_LEAN_EXACT=1
          export FILTERDDP_HYBRID_MU=${HYBRID_MU:?set HYBRID_MU for the hybrid arm};;
esac
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_TYPED_EQUATIONS=1 FILTERDDP_STRUCTURED_DYNAMICS=1
export FILTERDDP_TREE_KKT=1 FILTERDDP_BLOCKED_SOLVE=16
# OPENDSS_CHECK=1 (default since 2026-10-11): the run keeps its dispatch (written
# after the solve's timer has stopped) and, once Julia has exited, that dispatch
# is handed to OpenDSS before anything else runs: the verification belongs to
# the run, not to a later batch. OPENDSS_CHECK=0 skips both, as before.
if [ "${OPENDSS_CHECK:-1}" = 1 ]; then unset FILTERDDP_SKIP_SOLUTION_WRITE; else export FILTERDDP_SKIP_SOLUTION_WRITE=1; fi
export FILTERDDP_SCREEN=substation,vupper,ell,psubs FILTERDDP_LOADFLOW_START=1 FILTERDDP_AFFINE_LINESEARCH=1
export FILTERDDP_WARMUP=${FILTERDDP_WARMUP:-full} FILTERDDP_PARSIM=${PARSIM:-1}
export JULIA_NUM_THREADS=${JULIA_NUM_THREADS:-$THREADS}
LOG="$OUT/fddp_${SYS}_T${T}_${ARM}_${TAG}.log"
MARK="$OUT/.mark_$$"; touch "$MARK"
julia --startup-file=no --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl "$SYS" "$T" solve > "$LOG" 2>&1
if [ "${OPENDSS_CHECK:-1}" = 1 ]; then
  # Summary to stderr: callers that keep only the first line of stdout still see it.
  SOL=ddp/results/network_filterddp/filterddp_solution_${SYS}_T${T}_periodic_CBsystem.jls
  SOLDIR=ddp/results/kkt_ordering/captures/solutions; mkdir -p "$SOLDIR"
  RUNTAG="parsim_${SYS}_T${T}_${ARM}_${TAG}"
  if grep -q "solve complete" "$LOG" && [ "$SOL" -nt "$MARK" ]; then
    cp "$SOL" "$SOLDIR/sol_${RUNTAG}.jls"
    bash ddp/examples/power_system/run_opendss_check.sh "$SYS" "$T" "$SOLDIR/sol_${RUNTAG}.jls" "$RUNTAG" 2>&1 \
      | grep -E "OPENDSS_SUMMARY|FAILED" | cut -c1-330 >&2
  else
    echo "OPENDSS_CHECK: no solution written by this run ($RUNTAG)" >&2
  fi
fi
rm -f "$MARK"
echo "[$(date '+%H:%M:%S')] $SYS T=$T $ARM $TAG: $(grep -E 'solve complete' "$LOG" | tail -1) | $(grep -oE 'FILTERDDP_NEAR_OPT iteration=[0-9]+ elapsed_s=[0-9.]+' "$LOG" | head -1) | $(grep -oE 'FilterDDP objective=[-0-9.eE+]+' "$LOG")"
[ "${PARSIM:-1}" = 1 ] && python ddp/examples/power_system/parsim_from_log.py "$LOG"
