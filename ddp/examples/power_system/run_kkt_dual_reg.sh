#!/usr/bin/env bash
# Does an Ipopt-style constraint-block regularization, K = [H cu'; cu -δ_c I],
# rescue MA57/MA97 on FilterDDP's stage KKTs? Same captures and timing as
# run_hsl_benchmark.sh (diagonal Hessian, one thread), with KKT_DUAL_REG set;
# δ_c = 0 is the existing hsl_kkt_benchmark.csv. Each row records the delayed
# pivots, factor size and times, the residual against the regularized K, and
# dev_wide = ||X - X_exact|| / ||X_exact||, how far the regularized solution
# moves from the exact one.
#
#   bash ddp/examples/power_system/run_kkt_dual_reg.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
OUT=ddp/results/kkt_ordering
CAP=$OUT/captures
JL="julia --startup-file=no"
LOADJL=ddp/examples/power_system/sample_background_load.jl
mkdir -p "$OUT/dual_reg" "$OUT/logs"

for SYS in ieee123C_1ph ieee2522C_1ph large10kC_1ph; do
  case $SYS in ieee123C_1ph) R=21;; ieee2522C_1ph) R=9;; large10kC_1ph) R=3;; esac
  for D in 1e-12 1e-10 1e-8 1e-6; do
    CSV="$OUT/dual_reg/dualreg_${SYS}_diag_d${D}.csv"; LOG="$OUT/logs/dualreg_${SYS}_diag_d${D}.log"
    [ -s "$CSV" ] && { echo "skip $SYS d$D"; continue; }
    $JL "$LOADJL" wait 10 1.5 2>/dev/null | grep QUIET_WAIT > "$LOG"
    STOP="$OUT/logs/.sampling_dr_$$"; touch "$STOP"
    $JL "$LOADJL" sample "$STOP" "$OUT/logs/load_dualreg_${SYS}_diag_d${D}.csv" 15 > /dev/null 2>&1 &
    SPID=$!
    KKT_DUAL_REG=$D OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 \
      $JL --project=envs/ddp2026 ddp/examples/power_system/hsl_kkt_benchmark.jl \
      "$CAP/kkt_${SYS}_T3_diag_iter20_stage1.jls" "$SYS" diag "$CSV" "$R" umfpack,ma57,ma97 >> "$LOG" 2>&1
    rm -f "$STOP"; wait "$SPID" 2>/dev/null
    echo "[$(date '+%H:%M:%S')] $SYS d$D: $(tail -1 "$LOG")"
  done
done

MERGED=$OUT/dual_reg_benchmark.csv
first=1
for f in "$OUT"/dual_reg/dualreg_*.csv; do
  if [ $first = 1 ]; then cat "$f" > "$MERGED"; first=0; else tail -n +2 "$f" >> "$MERGED"; fi
done
echo "DUALREG_DONE merged $(($(wc -l < "$MERGED") - 1)) rows"
