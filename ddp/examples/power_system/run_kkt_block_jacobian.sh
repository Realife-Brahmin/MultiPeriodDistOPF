#!/usr/bin/env bash
# Block-diagonal Jacobian test (kkt_block_jacobian.jl) on the stage-1 captures
# of capture_kkt_magnitude.sh: mid-solve (iteration 40) and near-optimality.
#
#   bash ddp/examples/power_system/run_kkt_block_jacobian.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
CAPDIR=ddp/results/kkt_ordering/captures/magnitude
OUT=ddp/results/kkt_block_jacobian
mkdir -p "$OUT"
for cell in ieee123C_1ph:40 ieee123C_1ph:67 ieee2522C_1ph:40 ieee2522C_1ph:68 large10kC_1ph:40 large10kC_1ph:100; do
  S=${cell%%:*}; IT=${cell##*:}
  case $S in large10kC_1ph) KS=2,4,8,16,32,64,128;; *) KS=2,4,8,16,32,64;; esac
  OPENBLAS_NUM_THREADS=1 julia --startup-file=no --project=envs/ddp2026 ddp/examples/power_system/kkt_block_jacobian.jl \
    "$CAPDIR/kkt_${S}_T6_diag_cbsys_iter${IT}_stage1.jls" "ddp/results/network_filterddp/network_data_${S}_T6_periodic.jls" \
    "$S" "$IT" "$OUT/block_${S}_iter${IT}.csv" "$KS" 2>&1 | tail -1
done
first=1
for f in "$OUT"/block_*.csv; do
  if [ $first = 1 ]; then cat "$f" > "$OUT/block_jacobian.csv"; first=0; else tail -n +2 "$f" >> "$OUT/block_jacobian.csv"; fi
done
echo "BLOCK_JACOBIAN_DONE"
