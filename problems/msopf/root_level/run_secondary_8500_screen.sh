#!/usr/bin/env bash
set -euo pipefail
unset NO_COLOR FORCE_COLOR
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MULTIPERIOD_LEAN=1
export MULTI_VMARGIN=20000 MULTIPERIOD_MAX_WALL=180 MULTIPERIOD_SOLVE_LOG=1
export SECONDARY_LIMIT="${SECONDARY_LIMIT:-40}"
export SECONDARY_T="${SECONDARY_T:-24}"
hsl="/c/Users/Aryan Ritwajeet Jha/Documents/hsl/build"
if [[ -f "$hsl/libma57.dll" ]]; then
    export PATH="$hsl:$PATH"
    export IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=C:/Users/Aryan Ritwajeet Jha/Documents/hsl/build/libma57.dll;linear_system_scaling=none"
fi
root="${SECONDARY_ROOT:-problems/msopf/results/secondary_8500_oct25}"
if [[ $# == 0 ]]; then set -- secondary2 primary2 secondary3 primary3 secondary4 primary4; fi
for label in "$@"; do
    out="$root/$label"
    [[ ! -e "$out/metadata.txt" ]] || { echo "Refusing to overwrite $out" >&2; exit 2; }
    mkdir -p "$out"
    export SECONDARY_OUT="$out" MULTIPERIOD_OUT_DIR="$out"
    {
        date -u '+start_utc=%FT%TZ'
        git rev-parse HEAD
        git status --short
        hostname
        julia --version
        env | sort | grep -E '^(SECONDARY_|JULIA_|OPENBLAS|OMP_|MULTI|IPOPT_)'
        sha256sum problems/msopf/{Project.toml,Manifest.toml,full_angle_pf.jl,multiperiod_full_angle.jl} problems/msopf/root_level/{multiperiod_study.jl,secondary_8500_study.jl,run_secondary_8500_screen.sh}
        sha256sum test_networks/ieee8500{,xf}_1ph/*.dss
    } > "$out/metadata.txt"
    git diff -- problems/msopf > "$out/source.patch"
    t=$SECONDS
    status=0
    julia --startup-file=no --project=problems/msopf problems/msopf/root_level/secondary_8500_study.jl "$label" > "$out/run.log" 2> "$out/run.err" || status=$?
    echo "exit_code=$status" >> "$out/metadata.txt"
    echo "wall_seconds=$((SECONDS-t))" >> "$out/metadata.txt"
    echo "$label: exit=$status, wall=$((SECONDS-t)) s"
    [[ $status == 0 ]] || exit "$status"
done
