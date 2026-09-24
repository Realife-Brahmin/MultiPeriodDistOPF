#!/usr/bin/env bash
# Goal-post race: centralized Ipopt vs FilterDDP (diagonal Hessian), both on the
# IDENTICAL problem, at horizons long enough for Ipopt to hit its knee.
#
# Why (2026-09-18): the paper's centralized sweep shows large10k Ipopt jumping
# from 573 s at T = 24 to 8262 s at T = 48 (iterations 112 -> 698), and large10k
# T = 48 is the goal post FilterDDP should beat. That sweep ran on a different
# instance family (C_B ~ 8.8e-8, old price sampling, free terminal SOC), so it
# cannot be raced against. Here every pair is one problem:
#   - same exported instance: periodic profile (REDUCED_PROFILE=periodic)
#   - same C_B = 1e-3 (REDUCED_CB), fixed for every system and horizon
#   - same SOFT terminal SOC, gamma * sum (B^T - B^0)^2 with the per-system,
#     horizon-independent gamma in terminal_soc_penalty.jl (TERMINAL_SOC_SOFT=1)
# and it is PROVEN one problem: phase 1 stops everything unless the two solvers'
# objectives agree. (Checked by hand 2026-09-18: 4.1e-8 at T = 3, 4.2e-8 at T = 6.)
#
# NOT comparable with ddp/results/diag_hessian_horizon/, which predates the soft
# terminal SOC. The race therefore runs its own FilterDDP arms, logged here.
#
#   bash ddp/examples/power_system/run_matched_ipopt_race.sh [outdir]
#   bash ddp/examples/power_system/run_matched_ipopt_race.sh OUTDIR SYSTEM T
# The three-argument form runs/resumes exactly one matched case and stops.
#
# Phases, resumable (finished logs skipped), summary pushed after every case:
#   0  WAIT for run_diag_hessian_horizon_sweep.sh to finish (no timing overlap).
#   1  GATE  ieee123 T = 3, 6: Ipopt and FilterDDP-diag objectives within 1e-5.
#   2  pairs ieee123 T = 12 24 48 96
#   3  pairs ieee2522 T = 3 6 12 24 48 96
#   4  Ipopt large10k T = 3 6 12 24 48 (re-measures the goal post on this problem)
#   5  Ipopt ieee2522 T = 144 192 288 (look for ieee2522's knee)
#   6  FilterDDP large10k T = 48 (the goal post), then 24, 12, 6, 3
#   7  FilterDDP ieee2522 T = 144 192 288
# FilterDDP policy storage: FILTERDDP_FACTOR_BACKED_POLICY=1 for large10k T >= 12
# and ieee2522 T >= 144, where the dense per-stage policy would not fit in RAM.
#
# SHARED MACHINE (2026-09-19): the user runs their own analyses on this machine.
# Before each timed run the pipeline waits up to 10 min for other processes to
# drop below 1.5 busy cores, then runs regardless. During every run it records
# the cores used by processes that are NOT this pipeline (load_<run>.csv, via
# sample_background_load.jl), and the summary flags any run whose median
# background exceeded 1.5 cores. Idle floor on this machine is ~0.4 cores.

set -u
cd "$(dirname "$0")/../../.." || exit 1
REPO="$PWD"
OUT="${1:-$REPO/ddp/results/matched_ipopt_race}"
SWEEP="$REPO/ddp/results/diag_hessian_horizon"
JL="julia --startup-file=no"
DATA="ddp/results/network_filterddp"
BRANCH="$(git rev-parse --abbrev-ref HEAD)"
mkdir -p "$OUT/logs"

export REDUCED_PROFILE=periodic
export REDUCED_CB=1e-3
export TERMINAL_SOC_SOFT=1
export FILTERDDP_SKIP_SOLUTION_WRITE=1
export FILTERDDP_TIMING_DIAGNOSTIC=1
export FILTERDDP_MAX_ITERATIONS=400

say() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
LOADJL="ddp/examples/power_system/sample_background_load.jl"
quiet_wait() { $JL "$LOADJL" wait 10 1.5 2>/dev/null | grep QUIET_WAIT; }
# $1 = csv path; runs the rest of the arguments while sampling background load
sampled() {
  local csv=$1; shift
  local stop="$OUT/logs/.sampling_$$"
  touch "$stop"
  $JL "$LOADJL" sample "$stop" "$csv" 15 > /dev/null 2>&1 &
  local spid=$!
  "$@"
  rm -f "$stop"; wait "$spid" 2>/dev/null
}
done_marker() { [ -s "$1" ] && grep -q "PIPELINE_STEP_OK" "$1"; }

# ------------------------------------------------------------------- 0 WAIT --
waited=0
while [ -e "$SWEEP/.sweep_running.lock" ]; do
  [ "$waited" -eq 0 ] && say "waiting for the horizon sweep to finish"
  sleep 300; waited=$((waited + 300))
  if [ "$waited" -ge $((16 * 3600)) ]; then
    say "ABORT: horizon sweep still holds its lock after 16 h; not starting (would contaminate timings)"
    exit 1
  fi
done
say "horizon sweep finished; starting"

LOCK="$OUT/.race_running.lock"
touch "$LOCK"
$JL ddp/examples/power_system/hold_awake.jl "$LOCK" 40 > "$OUT/logs/hold_awake.log" 2>&1 &
trap 'rm -f "$LOCK"' EXIT

summarise_and_push() {
  $JL --project=envs/ddp2026 ddp/examples/power_system/near_opt_from_logs.jl \
      race "$OUT" > "$OUT/NEAR_OPT_SUMMARY.txt" 2>&1
  $JL --project=envs/ddp2026 ddp/examples/power_system/summarize_matched_race.jl \
      "$OUT" > "$OUT/RACE_SUMMARY.txt" 2>&1
  git add "$OUT/RACE_SUMMARY.txt" "$OUT/NEAR_OPT_SUMMARY.txt" "$OUT/near_opt.csv" \
      "$OUT/matched_race.csv" "$OUT/race_run.log" 2>/dev/null
  # Also stage every TRACKED file these pipelines touch (e.g. the horizon sweep's
  # log gets one last line after its final commit). An unstaged tracked change
  # makes `git pull --rebase` refuse, and every push after it would fail.
  git add -u -- "$OUT" "$SWEEP" 2>/dev/null
  if ! git diff --cached --quiet; then
    git commit -q -m "matched Ipopt race: $1 (auto-commit by run_matched_ipopt_race.sh)

Unreviewed pipeline output, regenerated from raw logs after this case.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
    if git pull -q --rebase origin "$BRANCH" && git push -q origin "$BRANCH"; then
      say "pushed after $1"
    else
      git rebase --abort 2>/dev/null
      say "push failed after $1 -- committed locally, retry next case"
    fi
  fi
}

ensure_export() {  # $1 system, $2 T
  local f="$DATA/network_data_$1_T$2_periodic.jls"
  [ -s "$f" ] && return 0
  say "export: $1 T=$2 (periodic profile)"
  PROFILE_PERIODIC=1 $JL --project=envs/tadmm ddp/examples/power_system/export_ieee123c_data.jl \
      "$1" "$2" > "$OUT/logs/export_$1_T$2.log" 2>&1
  [ -s "$f" ]
}

run_ipopt() {  # $1 system, $2 T
  local LOG="$OUT/logs/ipopt_$1_T$2.log"
  done_marker "$LOG" && { say "skip ipopt $1 T=$2"; return 0; }
  ensure_export "$1" "$2" || { say "export FAILED $1 T=$2"; return 1; }
  local QW; QW=$(quiet_wait)
  say "ipopt: $1 T=$2 ($QW)"
  echo "$QW" > "$LOG"
  sampled "$OUT/logs/load_ipopt_$1_T$2.csv" \
    $JL --project=envs/ddp2026 ddp/examples/power_system/centralized_ipopt_matched.jl \
        "$1" "$2" "$OUT/logs/ipopt_$1_T$2_ipoptlog.txt" >> "$LOG" 2>&1
  echo "PIPELINE_STEP_OK" >> "$LOG"
  summarise_and_push "ipopt $1 T=$2"
}

run_fddp() {  # $1 system, $2 T  -- diagonal-Hessian arm
  local SYS=$1 T=$2 LOG="$OUT/logs/fddp_diag_$1_T$2.log" FB=0
  done_marker "$LOG" && { say "skip filterddp $SYS T=$T"; return 0; }
  ensure_export "$SYS" "$T" || { say "export FAILED $SYS T=$T"; return 1; }
  [ "$SYS" = large10kC_1ph ] && [ "$T" -ge 12 ] && FB=1
  [ "$SYS" = ieee2522C_1ph ] && [ "$T" -ge 144 ] && FB=1
  [ "${RACE_FB_ALL:-0}" = 1 ] && FB=1
  local QW; QW=$(quiet_wait)
  say "filterddp diag: $SYS T=$T (factor_backed=$FB; $QW)"
  echo "$QW" > "$LOG"
  sampled "$OUT/logs/load_fddp_${SYS}_T${T}.csv" run_fddp_job "$SYS" "$T" "$FB" >> "$LOG" 2>&1
  local REF
  REF=$(grep -oE "CENTRAL_IPOPT .*" "$OUT/logs/ipopt_${SYS}_T${T}.log" | \
        grep -oE " objective=[-0-9.eE+]+" | cut -d= -f2)
  $JL --project=envs/ddp2026 ddp/examples/power_system/extract_filterddp_feasibility_trace.jl \
      "$LOG" "$OUT/fddp_diag_${SYS}_T${T}_trace.csv" "$REF" >> "$LOG" 2>&1
  echo "PIPELINE_STEP_OK" >> "$LOG"
  summarise_and_push "filterddp diag $SYS T=$T"
}

run_fddp_job() {  # $1 system, $2 T, $3 factor_backed
  local SYS=$1 T=$2 FB=$3
  local REF
  REF=$(grep -oE "CENTRAL_IPOPT .*" "$OUT/logs/ipopt_${SYS}_T${T}.log" | \
        grep -oE " objective=[-0-9.eE+]+" | cut -d= -f2)
  [ -n "$REF" ] || { echo "missing centralized objective for $SYS T=$T"; return 1; }
  (
    # Per-system primal infeasibility threshold for near-optimality.
    # These are the paper's reporting thresholds (Table V); do not change
    # without updating lean_results.tex and CLAUDE.md.
    local PRIMAL=1e-6
    [ "$SYS" = ieee2522C_1ph ] && PRIMAL=1e-5
    [ "$SYS" = large10kC_1ph ] && PRIMAL=1e-4
    export FILTERDDP_DIAG_HESSIAN=1 FILTERDDP_DIAG_HESSIAN_FLOOR=1e-8
    export FILTERDDP_NEAR_OPT_REFERENCE="$REF"
    export FILTERDDP_NEAR_OPT_GAP=0.005
    export FILTERDDP_NEAR_OPT_PRIMAL="$PRIMAL"
    export FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
    if [ "$FB" = 1 ]; then export FILTERDDP_FACTOR_BACKED_POLICY=1
    else unset FILTERDDP_FACTOR_BACKED_POLICY; fi
    echo "PIPELINE_ENV system=$SYS T=$T arm=diag factor_backed=$FB diag_floor=1e-8 terminal_soc_soft=$TERMINAL_SOC_SOFT max_iterations=$FILTERDDP_MAX_ITERATIONS C_B=$REDUCED_CB profile=$REDUCED_PROFILE near_opt_reference=$REF near_opt_gap=$FILTERDDP_NEAR_OPT_GAP near_opt_primal=$FILTERDDP_NEAR_OPT_PRIMAL feasibility_diagnostic=1 started=$(date '+%Y-%m-%dT%H:%M:%S')"
    $JL --project=envs/ddp2026 ddp/examples/power_system/ieee123c_filterddp.jl "$SYS" "$T" solve
  )
}

rel_obj_diff() {  # $1 system, $2 T -> relative objective difference, or nan
  local a b
  a=$(grep -oE "CENTRAL_IPOPT .*" "$OUT/logs/ipopt_$1_T$2.log" | grep -oE " objective=[-0-9.eE+]+" | cut -d= -f2)
  b=$(grep -oE "FilterDDP objective=[-0-9.eE+]+" "$OUT/logs/fddp_diag_$1_T$2.log" | cut -d= -f2)
  [ -n "$a" ] && [ -n "$b" ] || { echo nan; return; }
  python -c "a,b=float('$a'),float('$b'); print('%.3e' % (abs(a-b)/abs(a)))"
}

# One-case mode is used for supervised expensive runs. The centralized arm is
# skipped when its validated PIPELINE_STEP_OK log already exists.
if [ "$#" -ge 3 ]; then
  CASE_SYS=$2; CASE_T=$3
  run_ipopt "$CASE_SYS" "$CASE_T" || exit 1
  run_fddp "$CASE_SYS" "$CASE_T" || exit 1
  say "single case finished: $CASE_SYS T=$CASE_T"
  summarise_and_push "single case finished $CASE_SYS T=$CASE_T"
  exit 0
fi

# ------------------------------------------------------------------- 1 GATE --
for T in 3 6; do run_ipopt ieee123C_1ph $T; run_fddp ieee123C_1ph $T; done
g3=$(rel_obj_diff ieee123C_1ph 3); g6=$(rel_obj_diff ieee123C_1ph 6)
say "GATE: Ipopt vs FilterDDP objective, relative diff: T=3 $g3, T=6 $g6"
ok=$(python -c "import math
v=[float(x) for x in ('$g3','$g6')]
print('yes' if all(not math.isnan(x) and x < 1e-5 for x in v) else 'no')")
if [ "$ok" != yes ]; then
  say "ABORT: GATE FAILED -- Ipopt and FilterDDP do not solve the same problem; nothing below would be comparable"
  summarise_and_push "GATE FAILED"
  exit 1
fi

# ---------------------------------------------------------------- 2-3 PAIRS --
for T in 12 24 48 96;        do run_ipopt ieee123C_1ph  $T; run_fddp ieee123C_1ph  $T; done
for T in 3 6 12 24 48 96;    do run_ipopt ieee2522C_1ph $T; run_fddp ieee2522C_1ph $T; done
# ------------------------------------------------------------ 4-5 IPOPT ONLY --
for T in 3 6 12 24 48;       do run_ipopt large10kC_1ph $T; done
for T in 144 192 288;        do run_ipopt ieee2522C_1ph $T; done
# ---------------------------------------------------------- 6-7 FILTERDDP ONLY --
for T in 48 24 12 6 3;       do run_fddp large10kC_1ph $T; done
for T in 144 192 288;        do run_fddp ieee2522C_1ph $T; done

say "race finished"
summarise_and_push "race finished"
