#!/usr/bin/env bash
# KKT entry-magnitude study: analyse every capture of capture_kkt_magnitude.sh
# (kkt_magnitude_analysis.jl), plot full vs thresholded sparsity
# (plot_kkt_magnitude.py) and merge the per-capture CSVs.
#
#   bash ddp/examples/power_system/run_kkt_magnitude.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
CAPDIR=ddp/results/kkt_ordering/captures/magnitude
OUT=ddp/results/kkt_magnitude
PAT=$CAPDIR/patterns
mkdir -p "$OUT/csv" "$OUT/figures" "$PAT"
for f in "$CAPDIR"/kkt_*_T6_diag_cbsys_iter*_stage1.jls; do
  b=$(basename "$f"); S=$(echo "$b" | sed -E 's/kkt_(.*)_T6_diag.*/\1/'); IT=$(echo "$b" | sed -E 's/.*_iter([0-9]+)_stage1.jls/\1/')
  OPENBLAS_NUM_THREADS=1 julia --startup-file=no --project=envs/ddp2026 ddp/examples/power_system/kkt_magnitude_analysis.jl \
    "$f" "ddp/results/network_filterddp/network_data_${S}_T6_periodic.jls" "$S" "$IT" "$OUT/csv" "$PAT" 2>&1 | tail -1
done
python ddp/examples/power_system/plot_kkt_magnitude.py "$PAT" "$OUT/figures" \
  ieee123C_1ph:5,ieee123C_1ph:40,ieee123C_1ph:67,ieee2522C_1ph:5,ieee2522C_1ph:40,ieee2522C_1ph:68,large10kC_1ph:5,large10kC_1ph:40,large10kC_1ph:100
for kind in blocks threshold; do
  first=1
  for f in "$OUT"/csv/magnitude_${kind}_*.csv; do
    if [ $first = 1 ]; then cat "$f" > "$OUT/magnitude_${kind}.csv"; first=0; else tail -n +2 "$f" >> "$OUT/magnitude_${kind}.csv"; fi
  done
done
echo "KKT_MAGNITUDE_DONE"
