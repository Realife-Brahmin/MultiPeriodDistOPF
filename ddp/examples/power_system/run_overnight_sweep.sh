#!/usr/bin/env bash
# Overnight reduced-vs-full-space sweep over horizon and battery cost.
#
# Every cell runs BOTH formulations at the SAME C_B, because changing C_B
# changes the problem: a reduced solve compared against a reference at a
# different C_B is meaningless. Full-space runs first so its solution is on disk
# (tagged with C_B) and the reduced run reports the objective gap directly.
#
# Phases are ordered cheapest-and-most-informative first, so an interrupted run
# still leaves a usable table. Each cell has its own timeout and a failure is
# logged and stepped over rather than stopping the sweep.
#
# Results append to ddp/results/reduced_space/overnight/sweep.csv after every
# run, so nothing is lost if the machine or the shell dies.

set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO" || exit 1
OUT="$REPO/ddp/results/reduced_space/overnight"
LOGS="$OUT/logs"
CSV="$OUT/sweep.csv"
mkdir -p "$LOGS"

DEADLINE=$(( $(date +%s) + 11*3600 ))   # stop launching new work after 11 h

if [ ! -f "$CSV" ]; then
  echo "system,horizon,C_B,formulation,status,iterations,wall_s,objective_usd,obj_gap_rel,inner_solves,hessian_solves,inner_time_s,inner_pct,dual_inf,note" > "$CSV"
fi

log()  { echo "[$(date +%H:%M:%S)] $*"; }
left() { echo $(( DEADLINE - $(date +%s) )); }

# ---------------------------------------------------------------- exports ---
ensure_data() {
  local sys=$1 T=$2
  local f="$REPO/ddp/results/network_filterddp/network_data_${sys}_T${T}.jls"
  [ -f "$f" ] && return 0
  log "EXPORT $sys T=$T"
  timeout 3600 julia --startup-file=no --project=envs/tadmm \
    ddp/examples/power_system/export_ieee123c_data.jl "$sys" "$T" \
    > "$LOGS/export_${sys}_T${T}.log" 2>&1
  [ -f "$f" ] || { log "EXPORT FAILED $sys T=$T"; return 1; }
  return 0
}

# ------------------------------------------------------------- one cell ----
# run_cell <system> <T> <C_B> <timeout_full> <timeout_reduced>
run_cell() {
  local sys=$1 T=$2 cb=$3 tf=$4 tr=$5
  local tag="${sys}_T${T}_CB${cb}"

  ensure_data "$sys" "$T" || {
    echo "$sys,$T,$cb,export,FAILED,,,,,,,,,,no network data" >> "$CSV"; return 1; }

  # ---- full space ----
  if [ "$(left)" -gt 0 ]; then
    log "FULL  $tag  (timeout ${tf}s, $(( $(left)/60 )) min budget left)"
    local fl="$LOGS/full_${tag}.log"
    REDUCED_CB="$cb" FILTERDDP_MAX_ITERATIONS=200 \
      timeout "$tf" julia --startup-file=no --project=envs/ddp2026 \
      ddp/examples/power_system/ieee123c_filterddp.jl "$sys" "$T" solve quiet \
      > "$fl" 2>&1
    local rc=$?
    local line; line=$(grep -m1 "solve complete" "$fl" 2>/dev/null)
    if [ -n "$line" ]; then
      local w it st dual
      w=$(sed -E 's/.*solve complete: ([0-9.]+) s.*/\1/'   <<< "$line")
      it=$(sed -E 's/.*iterations=([0-9]+).*/\1/'          <<< "$line")
      st=$(sed -E 's/.*status=([^ ,]+).*/\1/'              <<< "$line")
      dual=$(grep -m1 "final residuals" "$fl" | sed -E 's/.*dual=([0-9.e+-]+).*/\1/')
      echo "$sys,$T,$cb,full,$st,$it,$w,,,,,,,$dual," >> "$CSV"
      log "  full: ${w}s it=$it status=$st"
    else
      echo "$sys,$T,$cb,full,TIMEOUT_OR_CRASH,,,,,,,,,,rc=$rc" >> "$CSV"
      log "  full: FAILED rc=$rc"
    fi
  fi

  # ---- reduced space ----
  if [ "$(left)" -gt 0 ]; then
    log "RED   $tag  (timeout ${tr}s)"
    local rl="$LOGS/reduced_${tag}.log"
    REDUCED_CB="$cb" REDUCED_HESS_RANK=25 REDUCED_HESS_OVERSAMPLE=10 \
      FILTERDDP_QUIET=1 FILTERDDP_MAX_ITERATIONS=200 \
      timeout "$tr" julia --startup-file=no --project=envs/ddp2026 \
      ddp/examples/power_system/reduced_space_filterddp.jl "$sys" "$T" lowrank \
      > "$rl" 2>&1
    local rc=$?
    local line; line=$(grep -m1 "^solve:" "$rl" 2>/dev/null)
    if [ -n "$line" ]; then
      local w it st obj gap is hs itime ipct dual
      w=$(sed -E 's/^solve: ([0-9.]+) s.*/\1/'      <<< "$line")
      it=$(sed -E 's/.*iterations=([0-9]+).*/\1/'   <<< "$line")
      st=$(sed -E 's/.*status=([^ ,]+).*/\1/'       <<< "$line")
      obj=$(grep -m1 "reduced-space objective" "$rl" | sed -E 's/.*: ([0-9.]+) USD.*/\1/')
      gap=$(grep -m1 "objective gap" "$rl" | sed -E 's/.*rel ([0-9.e+-]+).*/\1/')
      is=$(grep -m1 "inner solves=" "$rl" | sed -E 's/.*inner solves=([0-9]+).*/\1/')
      hs=$(grep -m1 "hessian solves=" "$rl" | sed -E 's/.*hessian solves=([0-9]+).*/\1/')
      itime=$(grep -m1 "inner time=" "$rl" | sed -E 's/.*inner time=([0-9.]+) s.*/\1/')
      ipct=$(grep -m1 "inner time=" "$rl" | sed -E 's/.*\(([0-9]+)% of wall\).*/\1/')
      dual=$(grep -m1 "^residuals:" "$rl" | sed -E 's/.*dual=([0-9.e+-]+).*/\1/')
      echo "$sys,$T,$cb,reduced,$st,$it,$w,$obj,$gap,$is,$hs,$itime,$ipct,$dual," >> "$CSV"
      log "  reduced: ${w}s it=$it status=$st gap=$gap inner=${ipct}%"
    else
      echo "$sys,$T,$cb,reduced,TIMEOUT_OR_CRASH,,,,,,,,,,rc=$rc" >> "$CSV"
      log "  reduced: FAILED rc=$rc"
    fi
  fi
}

# Smoke mode: one cheap ieee123 cell to prove the parsing before a long night
# is committed to it. A CSV row of empty fields means the regexes drifted.
if [ "${SWEEP_SMOKE:-0}" = "1" ]; then
  log "=== SMOKE TEST: ieee123 T=3 ==="
  run_cell ieee123C_1ph 3 1e-4 900 900
  column -s, -t "$CSV" 2>/dev/null || cat "$CSV"
  exit 0
fi

log "=== overnight sweep start; deadline in 11 h ==="

# Phase A -- ieee2522, C_B sweep at the two cheap horizons. This is where the
# "is C_B T-invariant?" question gets answered, so it runs first.
for T in 12 24; do
  for cb in 1e-5 1e-4 1e-3; do
    [ "$(left)" -gt 0 ] || break 2
    run_cell ieee2522C_1ph "$T" "$cb" 3600 5400
  done
done

# Phase B -- ieee2522, long horizons at the geometric-mean C_B only.
for T in 36 48; do
  [ "$(left)" -gt 0 ] || break
  run_cell ieee2522C_1ph "$T" 1e-4 5400 9000
done

# Phase C -- large10k. Most expensive, so last; T=12 first so at least one
# large-system matched pair lands.
for T in 12 24; do
  [ "$(left)" -gt 0 ] || break
  run_cell large10kC_1ph "$T" 1e-4 21600 14400
done

log "=== sweep finished or deadline reached ==="
column -s, -t "$CSV" 2>/dev/null || cat "$CSV"
