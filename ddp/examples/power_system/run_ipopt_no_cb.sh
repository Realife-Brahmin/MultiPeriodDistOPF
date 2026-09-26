#!/usr/bin/env bash
# Centralized Ipopt on every Table II instance WITHOUT the battery term
# (C_B = 0), everything else as in the matched race: periodic profile, soft
# terminal SOC with the per-system gamma, same exported instances, same Ipopt
# options. Question (meeting of 2026-09-25): does the C_B * P_B^2 term make the
# problem easier for the solvers? The C_B = 1e-3 counterparts are the race logs
# in ddp/results/matched_ipopt_race/logs/ipopt_<system>_T<T>.log.
#
# Background load is logged per case; one case at a time.
#
#   bash ddp/examples/power_system/run_ipopt_no_cb.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
OUT=ddp/results/ipopt_no_cb
mkdir -p "$OUT/logs"
export REDUCED_PROFILE=periodic REDUCED_CB=0 TERMINAL_SOC_SOFT=1
JL="julia --startup-file=no"
LOADJL=ddp/examples/power_system/sample_background_load.jl

for cell in ieee123C_1ph:6 ieee123C_1ph:24 ieee123C_1ph:96 \
            ieee2522C_1ph:6 ieee2522C_1ph:24 ieee2522C_1ph:96 \
            large10kC_1ph:6 large10kC_1ph:24 large10kC_1ph:48; do
  S=${cell%%:*}; T=${cell##*:}
  LOG="$OUT/logs/ipopt_nocb_${S}_T${T}.log"
  [ -s "$LOG" ] && grep -q "CENTRAL_IPOPT " "$LOG" && { echo "skip $S T=$T"; continue; }
  $JL "$LOADJL" wait 10 1.5 2>/dev/null | grep QUIET_WAIT > "$LOG"
  STOP="$OUT/logs/.sampling_$$"; touch "$STOP"
  $JL "$LOADJL" sample "$STOP" "$OUT/logs/load_nocb_${S}_T${T}.csv" 15 > /dev/null 2>&1 &
  SPID=$!
  $JL --project=envs/ddp2026 ddp/examples/power_system/centralized_ipopt_matched.jl \
      "$S" "$T" "$OUT/logs/ipopt_nocb_${S}_T${T}_ipoptlog.txt" >> "$LOG" 2>&1
  rm -f "$STOP"; wait "$SPID" 2>/dev/null
  echo "[$(date '+%H:%M:%S')] $(grep -oE 'CENTRAL_IPOPT system=\S+ T=[0-9]+ .*status=\S+ iterations=[0-9]+ objective=\S+ solve_time_s=\S+' "$LOG")"
done
echo "NOCB_DONE"
