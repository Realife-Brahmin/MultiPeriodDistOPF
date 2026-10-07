#!/usr/bin/env bash
# Tree solver (FILTERDDP_TREE_KKT=1 + FILTERDDP_STRUCTURED_DYNAMICS=1, typed
# residuals) beyond the Table II set, one run at a time on a quiet machine:
#   1. one Julia thread (like-for-like with single-threaded Ipopt);
#   2. repeats of the med2522 Table II cells;
#   3. longer horizons against the HSL references of ddp/results/ipopt_hsl
#      (same instances, per-system C_B), up to Ipopt's last horizon on this
#      machine (med2522 T=1152, large10k T=192);
#   4. med2522 T=1536, where MA57 and MA97 both ran out of memory: there is no
#      Ipopt objective, so FilterDDP runs to its own strict tolerance.
# Each run prints FILTERDDP_MEMORY_PEAK. Waits for the file given as $1 to
# contain VERIFY_DONE (or for no julia.exe if no argument).
#
#   bash ddp/examples/power_system/run_tree_long_horizon.sh [wait_file]

set -u
cd "$(dirname "$0")/../../.." || exit 1
if [ -n "${1:-}" ]; then until grep -q VERIFY_DONE "$1" 2>/dev/null; do sleep 30; done; fi
until [ "$(tasklist //FI "IMAGENAME eq julia.exe" 2>/dev/null | grep -c julia.exe)" = "0" ]; do sleep 30; done
echo "[$(date '+%H:%M:%S')] queue start"
S=ddp/examples/power_system/run_blocked_solve_fullrun.sh
H=ddp/results/ipopt_hsl/logs
export PIN_BLAS_THREADS=0 VARIANTS=blocked_w16 CB=system
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_TYPED_EQUATIONS=1
export FILTERDDP_STRUCTURED_DYNAMICS=1 FILTERDDP_TREE_KKT=1

# 1. one thread
export JULIA_NUM_THREADS=1 RUN_TAG_SUFFIX=_tableII_jt1_cbsys_typed_tree3
bash "$S" large10kC_1ph 6 diag 16 1
bash "$S" ieee2522C_1ph 24 diag 16 1
# 2. repeats of the med2522 cells
export JULIA_NUM_THREADS=8 RUN_TAG_SUFFIX=_tableII_jt8_cbsys_typed_tree3
for T in 6 24 96; do bash "$S" ieee2522C_1ph $T diag 16 2; done
# 3. longer horizons
export RUN_TAG_SUFFIX=_jt8_cbsys_typed_tree3
IPOPT_REF_LOG=$H/ipopt_ma57_ieee2522C_1ph_T192.log      bash "$S" ieee2522C_1ph 192 diag 16 1
IPOPT_REF_LOG=$H/ipopt_oom_ma57_ieee2522C_1ph_T384.log  bash "$S" ieee2522C_1ph 384 diag 16 1
IPOPT_REF_LOG=$H/ipopt_ma97_large10kC_1ph_T96.log       bash "$S" large10kC_1ph 96 diag 16 1
IPOPT_REF_LOG=$H/ipopt_oom_ma57_ieee2522C_1ph_T1152.log bash "$S" ieee2522C_1ph 1152 diag 16 1
IPOPT_REF_LOG=$H/ipopt_ma97_large10kC_1ph_T192.log      bash "$S" large10kC_1ph 192 diag 16 1
# 4. beyond Ipopt's memory limit: no reference, strict tolerance
LOG=ddp/results/kkt_ordering/fullrun_blocked/fddp_ieee2522C_1ph_T1536_diag_jt8_cbsys_typed_tree3_strict_r1.log
if ! grep -q "solve complete" "$LOG" 2>/dev/null; then
  REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1 FILTERDDP_MAX_ITERATIONS=400 \
  FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1 FILTERDDP_DIAG_HESSIAN=1 \
  FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8 FILTERDDP_BLOCKED_SOLVE=16 FILTERDDP_SKIP_SOLUTION_WRITE=1 \
    julia --startup-file=no --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl \
    ieee2522C_1ph 1536 solve > "$LOG" 2>&1
  echo "[$(date '+%H:%M:%S')] ieee2522C_1ph_T1536 strict: $(grep -E 'solve complete' "$LOG") | $(grep -oE 'FilterDDP objective=\S+' "$LOG") | $(grep -oE 'maxrss_mib=\S+' "$LOG")"
fi
echo "[$(date '+%H:%M:%S')] TREE_LONG_DONE"
