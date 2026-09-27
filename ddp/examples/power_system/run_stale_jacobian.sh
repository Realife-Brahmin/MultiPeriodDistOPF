#!/usr/bin/env bash
# Can the Jacobian refresh be skipped? FilterDDP with the KKT matrix's
# constraint Jacobian refreshed only every p-th iteration
# (FILTERDDP_STALE_JACOBIAN_PERIOD, backward_pass.jl) and everything else
# current, in Table II's configuration at the per-system C_B. p = 1 is the
# Table II run itself (*_tableII_*cbsys_*). Iterations and time to
# near-optimality against those.
#
#   bash ddp/examples/power_system/run_stale_jacobian.sh [cells] [periods]
#   e.g. bash ... "ieee123C_1ph:6 ieee2522C_1ph:6" "2 5 10"

set -u
cd "$(dirname "$0")/../../.." || exit 1
CELLS="${1:-ieee123C_1ph:6 ieee2522C_1ph:6}"
PERIODS="${2:-2 5 10}"
export PIN_BLAS_THREADS=0 VARIANTS=blocked_w16 CB=system
export FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1
export FILTERDDP_FACTOR_BACKED_POLICY=1
S=ddp/examples/power_system/run_blocked_solve_fullrun.sh
for cell in $CELLS; do
  SYS=${cell%%:*}; T=${cell##*:}
  if [ "$SYS" = ieee123C_1ph ]; then TH=1; else TH=8; fi
  for p in $PERIODS; do
    JULIA_NUM_THREADS=$TH FILTERDDP_STALE_JACOBIAN_PERIOD=$p RUN_TAG_SUFFIX=_stalejac${p}_cbsys \
      bash "$S" "$SYS" "$T" diag 16 1 | grep -E "solve complete|skip|no matched"
  done
done
echo "STALE_JACOBIAN_DONE"
