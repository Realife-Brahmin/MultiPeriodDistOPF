#!/usr/bin/env bash
# Re-run the FilterDDP cells of the paper's Table II (tab:near_opt) with the
# blocked multi-column solve, in Table II's configuration: diagonal Hessian,
# exact assembly rewrites, factor-backed policy (as the matched race used),
# Julia's default BLAS threads, per-system near-optimality thresholds. Only the
# blocked arm runs; the Ipopt column is unchanged. large10k T=6 was already run
# in this configuration (fullrun_blocked/*_rewrites_blasdefault_*).
#
# Ordered so the cheap cells finish first and large10k T=48 (~3.7 h) runs last.
#
#   bash ddp/examples/power_system/run_table2_blocked.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
export PIN_BLAS_THREADS=0 VARIANTS=blocked_w16 RUN_TAG_SUFFIX=_tableII
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1
S=ddp/examples/power_system/run_blocked_solve_fullrun.sh

for cell in ieee123C_1ph:6 ieee123C_1ph:24 ieee123C_1ph:96 \
            ieee2522C_1ph:6 ieee2522C_1ph:24 large10kC_1ph:24 \
            ieee2522C_1ph:96 large10kC_1ph:48; do
  bash "$S" "${cell%%:*}" "${cell##*:}" diag 16 1
done
echo "TABLE2_DONE"
