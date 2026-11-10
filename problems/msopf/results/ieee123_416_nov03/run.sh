#!/usr/bin/env bash
# ieee123_s14 at 4.16 kV (user, 2026-11-03), on the 309 PC. Usage: run.sh <label> <ratings 0|1> <script:phase:dir> ...
#   unrated A/AD give the currents for the blanket rating (1.2 x the largest); rated runs redo A/AD and certify (C, CD).
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="C:/Users/Aryan Ritwajeet Jha/Documents/documents_general/MultiPeriodDistOPF/.claude/worktrees/opendss-voltage-bases-check-162c57"
JL="C:/Users/Aryan Ritwajeet Jha/.julia/juliaup/julia-1.12.7+0.x64.w64.mingw32/bin/julia.exe"
ENVP="C:/Users/Aryan Ritwajeet Jha/msopf_envs/sdp_hsl"
HSL="C:/Users/Aryan Ritwajeet Jha/Documents/hsl/build"
OUT="$HERE/$1/${SYS:-ieee123_s14}"; mkdir -p "$OUT/logs"; RAT=$2; shift 2
export MULTIPERIOD_OUT_DIR="$OUT" MULTIPERIOD_SOLVE_LOG=0 FB_SYSTEM=${SYS:-ieee123_s14} FB_RATINGS=$RAT FB_IDLE="$OUT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 MKL_NUM_THREADS=1
export IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=$HSL/libma57.dll;ma57_pivot_order=2;ma57_pre_alloc=3.0;linear_system_scaling=none"
unset NO_COLOR FORCE_COLOR
for job in "$@"; do
    IFS=: read -r sc ph d <<< "$job"; log="$OUT/logs/${ph}_$d"
    { echo "system=${SYS:-ieee123_s14} ratings=$RAT phase=$ph dir=$d host=$(hostname) start=$(date -Is)"; echo "commit=$(git -C "$REPO" rev-parse HEAD) dirty=$(git -C "$REPO" status --porcelain test_networks/ieee123_5poi_1ph | wc -l)"; } > "$log.meta"
    (cd "$REPO" && FB_PHASE=$ph FB_DIR=$d "$JL" --startup-file=no --project="$ENVP" problems/msopf/root_level/$sc > "$log.log" 2> "$log.err" < /dev/null)
    rc=$?; echo "exit=$rc end=$(date -Is)" >> "$log.meta"; echo "$ph $d exit $rc: $(tail -1 "$OUT/fb_${ph}_$d.csv" 2>/dev/null)"
done
