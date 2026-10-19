#!/usr/bin/env bash
# Gurobi on the matched instances, on the footing of the clean batch of
# 2026-10-11 (run_clean_hp_batch.sh): the same JuMP model as the Ipopt runs
# (centralized_ipopt_matched.jl, CENTRAL_SOLVER=gurobi, envs/tadmm), with
#   - the reduced problem in the form Gurobi can read: it recognises the cone
#     P^2 + Q^2 <= v ell only while v and ell are declared nonnegative, so those
#     two stay and the other screened limits go (the driver does this when the
#     solver is Gurobi; BOUND_SCREEN in the log says what was dropped);
#   - the power-flow start values passed (IPOPT_LOADFLOW_START=1). Gurobi's
#     barrier method takes no starting point for a conic problem, so this is
#     not expected to change anything; START=0 runs without, to show it;
#   - the dispatch kept and handed to OpenDSS after every solve.
# Default threads (all cores) and one thread. One case at a time, after a quiet
# wait, background load logged, completed logs skipped. For reference: the
# paper's comparison is against Ipopt.
#
#   bash spar_ddp/examples/power_system/run_gurobi_hp.sh
#   CELLS="ieee123C_1ph:24" THREADS="1" START=0 bash spar_ddp/examples/power_system/run_gurobi_hp.sh
set -u
cd "$(dirname "$0")/../../.." || exit 1
unset NO_COLOR FORCE_COLOR
OUT=spar_ddp/results/gurobi_matched
mkdir -p "$OUT/logs"
SOLDIR=spar_ddp/results/kkt_ordering/captures/solutions; mkdir -p "$SOLDIR"
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1 CENTRAL_SOLVER=gurobi
export IPOPT_SCREEN=substation,vupper,ell,psubs
START=${START:-1}
[ "$START" = 1 ] && export IPOPT_LOADFLOW_START=1 || unset IPOPT_LOADFLOW_START
JL="julia --startup-file=no"
LOADJL=spar_ddp/examples/power_system/sample_background_load.jl
CELLS=${CELLS:-"ieee123C_1ph:6 ieee123C_1ph:24 ieee123C_1ph:96 ieee2522C_1ph:6 ieee2522C_1ph:24 ieee2522C_1ph:96 large10kC_1ph:6 large10kC_1ph:24 large10kC_1ph:48 large10kC_1ph:96"}
THREADS=${THREADS:-"0 1"}                      # 0 = Gurobi's default (all cores)
echo "GUROBI_HP start=$(date '+%Y-%m-%dT%H:%M:%S') plan=$(powercfg //getactivescheme | grep -oE '\(.*\)') start_values=$START"
for cell in $CELLS; do
  S=${cell%%:*}; T=${cell##*:}
  case $S in large10kC_1ph) PV=_pv;; *) PV="";; esac
  for TH in $THREADS; do
    [ "$TH" = 0 ] && TT=tdefault || TT=t$TH
    [ "$START" = 1 ] && ST=screen_lf || ST=screen
    TAG="gurobi_${TT}_${ST}${PV}_hp_${S}_T${T}"
    LOG="$OUT/logs/$TAG.log"
    [ -s "$LOG" ] && grep -q "CENTRAL_GUROBI " "$LOG" && { echo "skip $TAG"; continue; }
    $JL "$LOADJL" wait 10 1.5 2>/dev/null | grep QUIET_WAIT > "$LOG"
    STOP="$OUT/logs/.sampling_$$"; touch "$STOP"
    $JL "$LOADJL" sample "$STOP" "$OUT/logs/load_$TAG.csv" 15 > /dev/null 2>&1 &
    SPID=$!
    GUROBI_EXTRA_OPTIONS="Threads=$TH${GUROBI_MORE_OPTIONS:+;$GUROBI_MORE_OPTIONS}" IPOPT_SAVE_DISPATCH="$SOLDIR/sol_$TAG.jls" \
      $JL --project=envs/tadmm spar_ddp/examples/power_system/centralized_ipopt_matched.jl \
      "$S" "$T" "$OUT/logs/${TAG}_gurobilog.txt" >> "$LOG" 2>&1
    rm -f "$STOP"; wait "$SPID" 2>/dev/null
    echo "[$(date '+%H:%M:%S')] $TAG $(grep -oE 'CENTRAL_GUROBI .*status=\S+ iterations=[0-9]+ objective=\S+ solve_time_s=\S+' "$LOG" | grep -oE 'status=.*') | $(grep -oE 'Presolve removed [0-9]+ rows and [0-9]+ columns' "$LOG" | head -1)"
    if [ -s "$SOLDIR/sol_$TAG.jls" ]; then
      bash spar_ddp/examples/power_system/run_opendss_check.sh "$S" "$T" "$SOLDIR/sol_$TAG.jls" "$TAG" 2>&1 | grep -E "OPENDSS_SUMMARY|FAILED" | cut -c1-330
    else
      echo "OPENDSS_CHECK: no dispatch saved for $TAG"
    fi
  done
done
echo "GUROBI_HP_DONE"
