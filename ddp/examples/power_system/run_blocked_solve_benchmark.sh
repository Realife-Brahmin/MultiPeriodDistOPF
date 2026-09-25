#!/usr/bin/env bash
# Run blocked_multirhs_solve_benchmark.jl on every stage-1 capture: sequential
# (one Julia thread) for all six, then a separately labelled 8-thread run for
# med2522 and large10k. Background load is logged during every case.
#
#   bash ddp/examples/power_system/run_blocked_solve_benchmark.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
OUT=ddp/results/kkt_ordering
JL="julia --startup-file=no"
LOADJL=ddp/examples/power_system/sample_background_load.jl
mkdir -p "$OUT/blocked" "$OUT/logs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

run_case() {  # $1 system, $2 arm, $3 julia threads, $4 repeats
  local SYS=$1 ARM=$2 TH=$3 R=$4
  local TAG="${SYS}_${ARM}_t${TH}"
  local CSV="$OUT/blocked/blocked_${TAG}.csv" LOG="$OUT/logs/blocked_${TAG}.log"
  [ -s "$CSV" ] && { echo "skip $TAG"; return; }
  $JL "$LOADJL" wait 10 1.5 2>/dev/null | grep QUIET_WAIT > "$LOG"
  local STOP="$OUT/logs/.sampling_blk_$$"; touch "$STOP"
  $JL "$LOADJL" sample "$STOP" "$OUT/logs/load_blocked_${TAG}.csv" 15 > /dev/null 2>&1 &
  local SPID=$!
  $JL -t "$TH" --project=envs/ddp2026 ddp/examples/power_system/blocked_multirhs_solve_benchmark.jl \
      "$OUT/captures/kkt_${SYS}_T3_${ARM}_iter20_stage1.jls" "$SYS" "$ARM" "$CSV" "$R" >> "$LOG" 2>&1
  rm -f "$STOP"; wait "$SPID" 2>/dev/null
  echo "[$(date '+%H:%M:%S')] $TAG: $(grep -c blocked "$LOG") blocked rows"
}

for SYS in ieee123C_1ph ieee2522C_1ph large10kC_1ph; do
  case $SYS in ieee123C_1ph) R=21;; ieee2522C_1ph) R=9;; large10kC_1ph) R=3;; esac
  for ARM in diag exact; do run_case "$SYS" "$ARM" 1 "$R"; done
done
for SYS in ieee2522C_1ph large10kC_1ph; do
  case $SYS in ieee2522C_1ph) R=9;; large10kC_1ph) R=3;; esac
  for ARM in diag exact; do run_case "$SYS" "$ARM" 8 "$R"; done
done

MERGED=$OUT/blocked_multirhs_solve.csv
first=1
for f in "$OUT"/blocked/blocked_*.csv; do
  if [ $first = 1 ]; then cat "$f" > "$MERGED"; first=0; else tail -n +2 "$f" >> "$MERGED"; fi
done
echo "merged $(($(wc -l < "$MERGED") - 1)) rows into $MERGED"
