#!/usr/bin/env bash
# Push centralized Ipopt on med2522 (per-system C_B, periodic profile) to
# longer horizons until each linear solver runs out of memory: by default HSL
# MA57 and HSL MA97 on one thread (SOLVERS="mumps ma57 ma97" adds MUMPS,
# Ipopt's default). A solver is dropped from the
# sweep after its first failure. A watchdog kills any julia.exe whose private
# memory passes LIMIT_GB (default 27 of the PC's 32 GB) so that an allocation
# spiral ends as a recorded out-of-memory instead of paging the machine to a
# halt; such a run is logged as OOM_KILLED with the size it reached.
#
#   bash ddp/examples/power_system/run_ipopt_oom_sweep.sh [horizons]
#   e.g. bash ... "384 576 768 1152 1536"

set -u
cd "$(dirname "$0")/../../.." || exit 1
S=ieee2522C_1ph
HORIZONS="${1:-384 576 768 1152 1536}"
LIMIT_GB="${LIMIT_GB:-27}"
HSL="${HSL_LIB_DIR:-C:/Users/Aryan Ritwajeet Jha/Documents/hsl/build}"
OUT=ddp/results/ipopt_hsl
mkdir -p "$OUT/logs"
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export PATH="$HSL:$PATH"
JL="julia --startup-file=no"
ALIVE="${SOLVERS:-ma57 ma97}"

watchdog() {  # kill any julia.exe above LIMIT_GB private memory; note it in $1
  while [ -f "$1.running" ]; do
    powershell -NoProfile -Command "Get-Process julia -ErrorAction SilentlyContinue | Where-Object { \$_.PrivateMemorySize64 -gt ${LIMIT_GB}GB } | ForEach-Object { 'OOM_KILLED pid=' + \$_.Id + ' private_gb=' + [math]::Round(\$_.PrivateMemorySize64/1GB,1); Stop-Process -Id \$_.Id -Force }" >> "$1" 2>/dev/null
    sleep 10
  done
}

for T in $HORIZONS; do
  [ -z "$ALIVE" ] && break
  DATA="ddp/results/network_filterddp/network_data_${S}_T${T}_periodic.jls"
  if [ ! -s "$DATA" ]; then
    PROFILE_PERIODIC=1 $JL --project=envs/tadmm ddp/examples/power_system/export_ieee123c_data.jl "$S" "$T" \
      > "$OUT/logs/export_${S}_T${T}.log" 2>&1
    [ -s "$DATA" ] || { echo "export failed for T=$T"; break; }
  fi
  NEXT=""
  for LS in $ALIVE; do
    case $LS in
      mumps) OPTS="linear_system_scaling=none";;
      ma57)  OPTS="linear_solver=ma57;hsllib=$HSL/libma57.dll;linear_system_scaling=none";;
      ma97)  OPTS="linear_solver=ma97;hsllib=$HSL/libhsl_ma97.dll;linear_system_scaling=none";;
    esac
    LOG="$OUT/logs/ipopt_oom_${LS}_${S}_T${T}.log"
    touch "$LOG.running"; watchdog "$LOG" & WPID=$!
    IPOPT_EXTRA_OPTIONS="$OPTS" $JL --project=envs/ddp2026 ddp/examples/power_system/centralized_ipopt_matched.jl \
      "$S" "$T" "$OUT/logs/ipopt_oom_${LS}_${S}_T${T}_ipoptlog.txt" >> "$LOG" 2>&1
    rm -f "$LOG.running"; wait "$WPID" 2>/dev/null
    LINE=$(grep -oE 'status=\S+ iterations=[0-9]+ objective=\S+ solve_time_s=\S+' "$LOG")
    MEM=$(grep -oE 'maxrss_mib=\S+' "$LOG")
    if echo "$LINE" | grep -q "LOCALLY_SOLVED"; then
      echo "[$(date '+%H:%M:%S')] T=$T $LS: $LINE $MEM"; NEXT="$NEXT $LS"
    else
      echo "[$(date '+%H:%M:%S')] T=$T $LS: FAILED ($(grep -oE 'OOM_KILLED.*|status=\S+' "$LOG" | head -1); last: $(grep -iE 'error|memory|INFO' "$LOG" | tail -1 | cut -c1-120))"
    fi
  done
  ALIVE=$NEXT
done
echo "OOM_SWEEP_DONE"
