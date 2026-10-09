#!/usr/bin/env bash
# Gurobi on the matched instances (periodic profile, per-system C_B, soft
# terminal SOC): the same JuMP model as centralized_ipopt_matched.jl, solved
# with CENTRAL_SOLVER=gurobi from envs/tadmm (Gurobi.jl with its bundled
# Gurobi 13). For reference only; the paper's comparison is against Ipopt.
# Each cell runs with Gurobi's default thread count and with one thread
# (Ipopt's MUMPS, MA57 and MA97 figures are single-threaded). The model is
# left unscreened: Gurobi needs the lower bounds on v and ell to recognise the
# cone, and it applies its own presolve (statistics are in each log).
# One case at a time, background load logged, completed logs skipped.
#
#   bash ddp/examples/power_system/run_gurobi_matched.sh
#   CELLS="large10kC_1ph:96" THREADS="1" bash ddp/examples/power_system/run_gurobi_matched.sh

set -u
cd "$(dirname "$0")/../../.." || exit 1
OUT=ddp/results/gurobi_matched
mkdir -p "$OUT/logs"
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1 CENTRAL_SOLVER=gurobi
JL="julia --startup-file=no"
LOADJL=ddp/examples/power_system/sample_background_load.jl
CELLS=${CELLS:-"ieee123C_1ph:6 ieee123C_1ph:24 ieee123C_1ph:96 ieee2522C_1ph:6 ieee2522C_1ph:24 ieee2522C_1ph:96 large10kC_1ph:6 large10kC_1ph:24 large10kC_1ph:48"}
THREADS=${THREADS:-"0 1"}                      # 0 = Gurobi's default (all cores)
for cell in $CELLS; do
  S=${cell%%:*}; T=${cell##*:}
  for TH in $THREADS; do
    [ "$TH" = 0 ] && TAG=tdefault || TAG=t$TH
    LOG="$OUT/logs/gurobi_${TAG}_${S}_T${T}.log"
    [ -s "$LOG" ] && grep -q "CENTRAL_GUROBI " "$LOG" && { echo "skip $TAG $S T=$T"; continue; }
    $JL "$LOADJL" wait 10 1.5 2>/dev/null | grep QUIET_WAIT > "$LOG"
    STOP="$OUT/logs/.sampling_$$"; touch "$STOP"
    $JL "$LOADJL" sample "$STOP" "$OUT/logs/load_${TAG}_${S}_T${T}.csv" 15 > /dev/null 2>&1 &
    SPID=$!
    GUROBI_EXTRA_OPTIONS="Threads=$TH" \
      $JL --project=envs/tadmm ddp/examples/power_system/centralized_ipopt_matched.jl \
      "$S" "$T" "$OUT/logs/gurobi_${TAG}_${S}_T${T}_gurobilog.txt" >> "$LOG" 2>&1
    rm -f "$STOP"; wait "$SPID" 2>/dev/null
    echo "[$(date '+%H:%M:%S')] gurobi_$TAG $(grep -oE 'CENTRAL_GUROBI system=\S+ T=[0-9]+ .*status=\S+ iterations=[0-9]+ objective=\S+ solve_time_s=\S+' "$LOG") | $(grep -oE 'Presolve removed [0-9]+ rows and [0-9]+ columns' "$LOG" | head -1)"
  done
done
echo "GUROBI_MATCHED_DONE"
