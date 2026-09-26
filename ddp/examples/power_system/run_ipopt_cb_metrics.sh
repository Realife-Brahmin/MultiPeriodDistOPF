#!/usr/bin/env bash
# Centralized Ipopt on every Table II instance WITH the battery term
# (C_B = 1e-3), exactly as the matched race ran it. Rerun to record battery
# utilization (CENTRAL_IPOPT_BATTERY), which the race logs predate, as the
# counterpart of run_ipopt_no_cb.sh (C_B = 0). Writes only to
# ddp/results/ipopt_cb_1e-3/; the race logs are not touched.
#
# Background load is logged per case; one case at a time.
#
#   bash ddp/examples/power_system/run_ipopt_cb_metrics.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
OUT=ddp/results/ipopt_cb_1e-3
mkdir -p "$OUT/logs"
export REDUCED_PROFILE=periodic REDUCED_CB=1e-3 TERMINAL_SOC_SOFT=1
JL="julia --startup-file=no"
LOADJL=ddp/examples/power_system/sample_background_load.jl

for cell in ieee123C_1ph:6 ieee123C_1ph:24 ieee123C_1ph:96 \
            ieee2522C_1ph:6 ieee2522C_1ph:24 ieee2522C_1ph:96 \
            large10kC_1ph:6 large10kC_1ph:24 large10kC_1ph:48; do
  S=${cell%%:*}; T=${cell##*:}
  LOG="$OUT/logs/ipopt_cb1e-3_${S}_T${T}.log"
  [ -s "$LOG" ] && grep -q "CENTRAL_IPOPT " "$LOG" && { echo "skip $S T=$T"; continue; }
  $JL "$LOADJL" wait 10 1.5 2>/dev/null | grep QUIET_WAIT > "$LOG"
  STOP="$OUT/logs/.sampling_$$"; touch "$STOP"
  $JL "$LOADJL" sample "$STOP" "$OUT/logs/load_cb1e-3_${S}_T${T}.csv" 15 > /dev/null 2>&1 &
  SPID=$!
  $JL --project=envs/ddp2026 ddp/examples/power_system/centralized_ipopt_matched.jl \
      "$S" "$T" "$OUT/logs/ipopt_cb1e-3_${S}_T${T}_ipoptlog.txt" >> "$LOG" 2>&1
  rm -f "$STOP"; wait "$SPID" 2>/dev/null
  echo "[$(date '+%H:%M:%S')] $(grep -oE 'CENTRAL_IPOPT system=\S+ T=[0-9]+ .*status=\S+ iterations=[0-9]+ objective=\S+ solve_time_s=\S+' "$LOG")"
done
echo "CB_DONE"
