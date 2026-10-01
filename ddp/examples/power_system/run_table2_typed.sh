#!/usr/bin/env bash
# Re-run the FilterDDP cells of the paper's Tables II/III with the typed
# residual callback (FILTERDDP_TYPED_EQUATIONS=1), otherwise in exactly the
# configuration of the per-system-C_B runs: diagonal Hessian, blocked solve,
# factor-backed policy, exact assembly rewrites, ieee123 on one thread and
# med2522/large10k on eight (run_table2_blocked.sh). Cells whose log is
# complete are skipped; ieee123 T=6 and med2522 T=24 get a second repeat
# because their first ran while the machine was in use, and large10k T=24
# already has a clean repeat (r2).
#
# Finally, one battery-block run (FILTERDDP_BATTERY_SCHUR=1, MUMPS analysis
# reused) at large10k T=24, for a clean comparison against the typed run.
# Starts only when no julia.exe is running.
#
#   bash ddp/examples/power_system/run_table2_typed.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
until [ "$(tasklist //FI "IMAGENAME eq julia.exe" 2>/dev/null | grep -c julia.exe)" = "0" ]; do sleep 30; done
echo "[$(date '+%H:%M:%S')] queue start"
S=ddp/examples/power_system/run_blocked_solve_fullrun.sh
export PIN_BLAS_THREADS=0 VARIANTS=blocked_w16 CB=system
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_TYPED_EQUATIONS=1

# system T repeats  (ieee123 on one thread, as in Table II)
export JULIA_NUM_THREADS=1 RUN_TAG_SUFFIX=_tableII_cbsys_typed
for cell in "ieee123C_1ph 6 2" "ieee123C_1ph 24 1" "ieee123C_1ph 96 1"; do
  set -- $cell; bash "$S" "$1" "$2" diag 16 "$3"
done
export JULIA_NUM_THREADS=8 RUN_TAG_SUFFIX=_tableII_jt8_cbsys_typed
for cell in "ieee2522C_1ph 6 1" "ieee2522C_1ph 24 2" "ieee2522C_1ph 96 1" \
            "large10kC_1ph 6 1" "large10kC_1ph 48 1"; do
  set -- $cell; bash "$S" "$1" "$2" diag 16 "$3"
done

export FILTERDDP_BATTERY_SCHUR=1 RUN_TAG_SUFFIX=_tableII_jt8_cbsys_typed_schurreuse
bash "$S" large10kC_1ph 24 diag 16 1
echo "[$(date '+%H:%M:%S')] TABLE2_TYPED_DONE"
