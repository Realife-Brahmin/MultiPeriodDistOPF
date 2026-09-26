#!/usr/bin/env bash
# HSL MA57 / HSL_MA97 against UMFPACK and the blocked solve, on the diagonal-
# Hessian stage-1 captures of all three systems (hsl_kkt_benchmark.jl).
#   t1: everything on one thread (OpenMP, OpenBLAS and Julia threads = 1)
#   t8: the parallel solvers only -- MA97 with 8 OpenMP threads and the blocked
#       solve with 8 Julia threads; OpenBLAS stays at 1 so the parallelism is
#       the solvers' own.
# Needs the locally built HSL libraries in HSL_LIB_DIR (not in this repo).
#
#   bash ddp/examples/power_system/run_hsl_benchmark.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
OUT=ddp/results/kkt_ordering
CAP=$OUT/captures
JL="julia --startup-file=no"
LOADJL=ddp/examples/power_system/sample_background_load.jl
mkdir -p "$OUT/hsl" "$OUT/logs"

run() {  # $1 system $2 threads $3 repeats $4 families
  local SYS=$1 TH=$2 R=$3 FAM=$4
  local CSV="$OUT/hsl/hsl_${SYS}_diag_t${TH}.csv" LOG="$OUT/logs/hsl_${SYS}_diag_t${TH}.log"
  [ -s "$CSV" ] && { echo "skip $SYS t$TH"; return; }
  $JL "$LOADJL" wait 10 1.5 2>/dev/null | grep QUIET_WAIT > "$LOG"
  local STOP="$OUT/logs/.sampling_hsl_$$"; touch "$STOP"
  $JL "$LOADJL" sample "$STOP" "$OUT/logs/load_hsl_${SYS}_diag_t${TH}.csv" 15 > /dev/null 2>&1 &
  local SPID=$!
  OMP_NUM_THREADS=$TH OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=$TH \
    $JL --project=envs/ddp2026 ddp/examples/power_system/hsl_kkt_benchmark.jl \
    "$CAP/kkt_${SYS}_T3_diag_iter20_stage1.jls" "$SYS" diag "$CSV" "$R" "$FAM" >> "$LOG" 2>&1
  rm -f "$STOP"; wait "$SPID" 2>/dev/null
  echo "[$(date '+%H:%M:%S')] $SYS t$TH: $(tail -1 "$LOG")"
}

for SYS in ieee123C_1ph ieee2522C_1ph large10kC_1ph; do
  case $SYS in ieee123C_1ph) R=21;; ieee2522C_1ph) R=9;; large10kC_1ph) R=3;; esac
  run "$SYS" 1 "$R" umfpack,blocked,ma57,ma97
  run "$SYS" 8 "$R" blocked,ma97
done

MERGED=$OUT/hsl_kkt_benchmark.csv
first=1
for f in "$OUT"/hsl/hsl_*.csv; do
  if [ $first = 1 ]; then cat "$f" > "$MERGED"; first=0; else tail -n +2 "$f" >> "$MERGED"; fi
done
echo "HSL_DONE merged $(($(wc -l < "$MERGED") - 1)) rows"
