#!/usr/bin/env bash
# Run 4 more T=24 BF reps (runs 3, 4, 5, 6). Already have run1=287.8s, run2=272.3s.
# Backs up each result after the run so the next iteration doesn't clobber.
set -euo pipefail

cd /c/repos_addendum/MultiPeriodDistOPF

DIR=envs/tadmm/processedData/ieee2552C_1ph_T24

for REP in 3 4 5 6; do
  echo "=================================================================="
  echo "[$(date +%H:%M:%S)] Starting T=24 BF rep $REP"
  echo "=================================================================="
  SYSTEM_OVERRIDE=ieee2552C_1ph \
  T_OVERRIDE=24 \
    julia --project=envs/tadmm --threads=16 run_bf.jl
  cp "$DIR/results_socp_bf.txt" "$DIR/results_socp_bf_run${REP}.txt"
  cp "$DIR/ipopt_bf.log"        "$DIR/ipopt_bf_run${REP}.log"
  WT=$(grep 'Wall-clock' "$DIR/results_socp_bf.txt" | awk '{print $3}')
  echo "[$(date +%H:%M:%S)] Done rep $REP: WT=${WT}s"
done

echo
echo "=== Summary across all 6 runs ==="
for R in 1 2 3 4 5 6; do
  f="$DIR/results_socp_bf_run${R}.txt"
  [ -f "$f" ] && printf "run%d: %s s\n" "$R" "$(grep 'Wall-clock' "$f" | awk '{print $3}')"
done
