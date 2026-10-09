#!/usr/bin/env bash
# OpenDSS feasibility check (opendss_check.jl) of every saved FilterDDP solution
# whose run tag matches a pattern: the routine verification step for a batch of
# results. Solutions are the per-run copies kept by run_blocked_solve_fullrun.sh
# in ddp/results/kkt_ordering/captures/solutions/. Checks already done are
# skipped. Not a timed step; run it when no timed run is going.
#
#   bash ddp/examples/power_system/run_opendss_check_all.sh              # current configuration
#   bash ddp/examples/power_system/run_opendss_check_all.sh '_screen_warm_pv_'
set -u
cd "$(dirname "$0")/../../.." || exit 1
PATTERN="${1:-_screen_warm}"
SOL=ddp/results/kkt_ordering/captures/solutions
OUT=ddp/results/opendss_check
for f in "$SOL"/sol_*"$PATTERN"*.jls; do
  [ -e "$f" ] || continue
  tag=$(basename "$f" .jls); tag=${tag#sol_}
  sys=$(echo "$tag" | grep -oE '^[A-Za-z0-9]+_1ph')
  T=$(echo "$tag" | grep -oE '_T[0-9]+_' | head -1 | tr -d '_T')
  # large10k solutions from before the PV fix belong to an instance that no longer exists
  case "$tag" in large10kC_1ph*) case "$tag" in *_pv_*) ;; *) continue;; esac;; esac
  [ -s "$OUT/$tag.log" ] && grep -q OPENDSS_SUMMARY "$OUT/$tag.log" && continue
  bash ddp/examples/power_system/run_opendss_check.sh "$sys" "$T" "$f" "$tag" | grep -E "voltage:|FAILED" | cut -c1-260
done
echo "OPENDSS_CHECK_ALL_DONE"
