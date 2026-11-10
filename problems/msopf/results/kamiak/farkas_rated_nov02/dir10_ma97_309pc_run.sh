#!/usr/bin/env bash
# med3980_s123 rated, direction 10, on the 309 PC (2026-11-03): phase A with Ipopt on HSL MA97 (on MA57 its hour-3 solves
# near 1 deg returned LOCALLY_SOLVED but failed the balance check, so the day's reach came out 1.0 deg against a 5.57 deg
# hour boundary), then phase C as on Kamiak (default settings). Same settings as ~/msopf/jobs/farkas_rated_nov02/farkas_job.sh;
# ratings of ieee9500xf unchanged since e1ca54ec (0 switches raised by 9e310532/335a30e3).
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="C:/Users/Aryan Ritwajeet Jha/Documents/documents_general/MultiPeriodDistOPF/.claude/worktrees/opendss-voltage-bases-check-162c57"
JL="C:/Users/Aryan Ritwajeet Jha/.julia/juliaup/julia-1.12.7+0.x64.w64.mingw32/bin/julia.exe"
ENVP="C:/Users/Aryan Ritwajeet Jha/msopf_envs/sdp_hsl"
HSL="C:/Users/Aryan Ritwajeet Jha/Documents/hsl/build"
OUT="$HERE/med3980_s123"; mkdir -p "$OUT/logs"
export MULTIPERIOD_OUT_DIR="$OUT" MULTIPERIOD_SOLVE_LOG=0 FB_SYSTEM=med3980_s123 FB_DIR=10 FB_RATINGS=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 MKL_NUM_THREADS=1
unset NO_COLOR FORCE_COLOR
run() {
    local p=$1 log="$OUT/logs/$1_10"
    { echo "system=med3980_s123 phase=$p dir=10 host=$(hostname) start=$(date -Is)"; echo "ipopt_extra=$IPOPT_EXTRA_OPTIONS"
      echo "commit=$(git -C "$REPO" rev-parse HEAD)"; } > "$log.meta"
    (cd "$REPO" && FB_PHASE=$p "$JL" --startup-file=no --project="$ENVP" problems/msopf/root_level/farkas_bracketing.jl > "$log.log" 2> "$log.err" < /dev/null)
    local rc=$?; echo "exit=$rc end=$(date -Is)" >> "$log.meta"; return $rc
}
export IPOPT_EXTRA_OPTIONS="linear_solver=ma97;hsllib=$HSL/libhsl_ma97.dll;ma97_order=amd;linear_system_scaling=none"
run A && [ "$(wc -l < "$OUT/fb_A_10.csv")" -ge 2 ] && echo "A: $(tail -1 "$OUT/fb_A_10.csv")" || { echo "A failed"; exit 1; }
export IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=$HSL/libma57.dll;ma57_pivot_order=2;ma57_pre_alloc=3.0;linear_system_scaling=none"
run C && [ "$(wc -l < "$OUT/fb_C_10.csv")" -ge 2 ] && echo "C: $(tail -1 "$OUT/fb_C_10.csv")" || { echo "C failed"; exit 1; }
