#!/usr/bin/env bash
# Barrier schedule of FilterDDP, and the hybrid Hessian: iteration counts only.
#
# The solver lowers the barrier parameter by
#     mu <- max(tol / 10, min(kappa_mu * mu, mu ^ theta_mu))
# once the barrier problem is solved to kappa_eps * mu, starting from mu_init.
# Its defaults are mu_init = 1, kappa_mu = 0.2, theta_mu = 1.2, kappa_eps = 10
# (9 levels down to 1e-6); Ipopt's are 0.1, 0.2, 1.5, 10 (6 levels).
#
# A variant is name:mu_init:kappa_mu:theta_mu:kappa_eps[:hybrid_mu]. With the
# sixth field the run starts with the diagonal Hessian and switches to the
# exact one when mu first falls to hybrid_mu or below (FILTERDDP_HYBRID_MU;
# such a variant runs once, not per arm).
#
# Untimed: no warm-up, no accounting, several of these may run at once.
# Tolerance 1e-5, current configuration otherwise (run_parsim.sh), factors
# reused after a barrier update (FILTERDDP_REUSE_BARRIER=1, same iterates).
# Logs: spar_ddp/results/barrier_sweep/logs/bar_<system>_T<T>_<arm>_<name>.log;
# completed logs are skipped. Summary: barrier_sweep_summary.py.
#
#   bash spar_ddp/examples/power_system/run_barrier_sweep.sh
#   CELLS="ieee123C_1ph:24" ARMS=diag VARIANTS="t15:1:0.2:1.5:10" bash .../run_barrier_sweep.sh
set -u
cd "$(dirname "$0")/../../.." || exit 1
unset NO_COLOR FORCE_COLOR
OUT=spar_ddp/results/barrier_sweep/logs; mkdir -p "$OUT"
CELLS=${CELLS:-"ieee123C_1ph:24 ieee123C_1ph:96 ieee2522C_1ph:6 ieee2522C_1ph:24 large10kC_1ph:6 large10kC_1ph:24"}
ARMS=${ARMS:-"diag exact"}
VARIANTS=${VARIANTS:-"base:1:0.2:1.2:10 t15:1:0.2:1.5:10 m01:0.1:0.2:1.2:10 m01t15:0.1:0.2:1.5:10 k01t15:1:0.1:1.5:10 m01k01t20:0.1:0.1:2.0:10 m001t15:0.01:0.2:1.5:10 e30:1:0.2:1.2:30"}
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
export FILTERDDP_MAX_ITERATIONS=400 FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
export FILTERDDP_OPTIMALITY_TOLERANCE=1e-5
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_TYPED_EQUATIONS=1 FILTERDDP_STRUCTURED_DYNAMICS=1
export FILTERDDP_TREE_KKT=1 FILTERDDP_BLOCKED_SOLVE=16 FILTERDDP_SKIP_SOLUTION_WRITE=1
export FILTERDDP_SCREEN=substation,vupper,ell,psubs FILTERDDP_LOADFLOW_START=1 FILTERDDP_AFFINE_LINESEARCH=1
export FILTERDDP_LEAN_VALUE=1 FILTERDDP_HALF_BLOCK=1 FILTERDDP_REUSE_BARRIER=1
unset FILTERDDP_WARMUP FILTERDDP_PARSIM FILTERDDP_PINT_SCAN FILTERDDP_NEAR_OPT_REFERENCE
for V in $VARIANTS; do
  IFS=: read -r NAME MU0 KMU TMU KEPS HYB <<< "$V"
  export FILTERDDP_MU_INIT=$MU0 FILTERDDP_KAPPA_MU=$KMU FILTERDDP_THETA_MU=$TMU FILTERDDP_KAPPA_EPS=$KEPS
  [ -n "${HYB:-}" ] && RUNARMS=hybrid || RUNARMS=$ARMS
  for ARM in $RUNARMS; do
    unset FILTERDDP_DIAG_HESSIAN FILTERDDP_DIAG_HESSIAN_FLOOR FILTERDDP_LEAN_EXACT FILTERDDP_HYBRID_MU
    case $ARM in
      diag)   export FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8;;
      exact)  export FILTERDDP_LEAN_EXACT=1;;
      hybrid) export FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8 FILTERDDP_LEAN_EXACT=1 FILTERDDP_HYBRID_MU=$HYB;;
    esac
    for cell in $CELLS; do
      S=${cell%%:*}; T=${cell##*:}
      case $S in ieee123C_1ph) export JULIA_NUM_THREADS=1;; *) export JULIA_NUM_THREADS=${SWEEP_THREADS:-8};; esac
      LOG="$OUT/bar_${S}_T${T}_${ARM}_${NAME}.log"
      [ -s "$LOG" ] && grep -q "solve complete" "$LOG" && continue
      julia --startup-file=no --project=envs/spar_ddp spar_ddp/examples/power_system/ieee123c_filterddp.jl "$S" "$T" solve > "$LOG" 2>&1
      echo "[$(date '+%H:%M:%S')] $S T=$T $ARM $NAME: $(grep -E 'solve complete' "$LOG" | tail -1) | $(grep -oE 'FilterDDP objective=[-0-9.eE+]+' "$LOG") | $(grep -oE 'BARRIER_OPTIONS.*' "$LOG" | head -1) $(grep -oE 'FILTERDDP_HYBRID_SWITCH.*' "$LOG" | head -1)"
    done
  done
done
echo "BARRIER_SWEEP_DONE"
