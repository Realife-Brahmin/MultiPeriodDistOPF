#!/usr/bin/env bash
# Run kkt_ordering_benchmark.jl on every captured stage-1 KKT system and merge
# the results into ddp/results/kkt_ordering/kkt_ordering_benchmark.csv.
#
# Sequential throughout (BLAS, OpenMP, OpenBLAS pinned to one thread). The machine
# is shared: before each case it waits up to 10 min for other processes to drop
# below 1.5 busy cores, and during each case it logs the cores used by everything
# else (load_bench_<system>_<arm>.csv), as the matched race does.
#
#   bash ddp/examples/power_system/run_kkt_ordering_benchmark.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
OUT=ddp/results/kkt_ordering
CAP=$OUT/captures
JL="julia --startup-file=no"
LOADJL=ddp/examples/power_system/sample_background_load.jl
mkdir -p "$OUT/logs" "$OUT/bench"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

for SYS in ieee123C_1ph ieee2522C_1ph large10kC_1ph; do
  case $SYS in ieee123C_1ph) R=21;; ieee2522C_1ph) R=9;; large10kC_1ph) R=3;; esac
  for ARM in diag exact; do
    IN="$CAP/kkt_${SYS}_T3_${ARM}_iter20_stage1.jls"
    CSV="$OUT/bench/bench_${SYS}_${ARM}.csv"
    LOG="$OUT/logs/bench_${SYS}_${ARM}.log"
    [ -s "$CSV" ] && { echo "skip $SYS $ARM"; continue; }
    QW=$($JL "$LOADJL" wait 10 1.5 2>/dev/null | grep QUIET_WAIT)
    echo "[$(date '+%H:%M:%S')] $SYS $ARM repeats=$R $QW"
    echo "$QW" > "$LOG"
    STOP="$OUT/logs/.sampling_$$"; touch "$STOP"
    $JL "$LOADJL" sample "$STOP" "$OUT/logs/load_bench_${SYS}_${ARM}.csv" 15 > /dev/null 2>&1 &
    SPID=$!
    $JL --project=envs/ddp2026 ddp/examples/power_system/kkt_ordering_benchmark.jl \
        "$IN" "$SYS" "$ARM" "$CSV" "$R" >> "$LOG" 2>&1
    rm -f "$STOP"; wait "$SPID" 2>/dev/null
    tail -1 "$LOG"
  done
done

MERGED=$OUT/kkt_ordering_benchmark.csv
first=1
for f in "$OUT"/bench/bench_*.csv; do
  if [ $first = 1 ]; then cat "$f" > "$MERGED"; first=0; else tail -n +2 "$f" >> "$MERGED"; fi
done
echo "merged $(($(wc -l < "$MERGED") - 1)) rows into $MERGED"
