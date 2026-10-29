#!/usr/bin/env bash
set -euo pipefail
unset NO_COLOR FORCE_COLOR
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MULTIPERIOD_LEAN=1
export MULTI_VMARGIN=20000 MULTIPERIOD_MAX_WALL="${MATCHED_WALL:-360}" MULTIPERIOD_SOLVE_LOG=1 MULTIPERIOD_FAILED_SNAPSHOT="${MATCHED_FAILED_SNAPSHOT:-0}" SECONDARY_T=24 SECONDARY_LIMIT=60
hsl='/c/Users/Aryan Ritwajeet Jha/Documents/hsl/build'
export PATH="$hsl:$PATH"
export IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=C:/Users/Aryan Ritwajeet Jha/Documents/hsl/build/libma57.dll;linear_system_scaling=none;${MATCHED_IPOPT_OPTIONS:-}"
export MATCHED_ROOT="${MATCHED_ROOT:-problems/msopf/results/matched_secondary_oct25}"
label="$1"; phase="$2"
folder="${phase//,/x}"
export SECONDARY_OUT="$MATCHED_ROOT/$label/$folder" MULTIPERIOD_OUT_DIR="$MATCHED_ROOT/$label/$folder"
[[ ! -e "$SECONDARY_OUT/metadata.txt" ]] || { echo "Refusing overwrite $SECONDARY_OUT";exit 2; }
mkdir -p "$SECONDARY_OUT"
{
 date -u '+start_utc=%FT%TZ'; git rev-parse HEAD; hostname; julia --version
 env | sort | grep -E '^(MATCHED_|SECONDARY_|JULIA_|OPENBLAS|OMP_|MULTI|IPOPT_)'
 sha256sum problems/msopf/{Project.toml,Manifest.toml,full_angle_pf.jl,multiperiod_full_angle.jl} problems/msopf/root_level/{matched_secondary_study.jl,secondary_8500_study.jl,multiperiod_study.jl,run_matched_secondary.sh}
 sha256sum test_networks/ieee9500{,xf}_1ph/*.dss test_networks/smartds_p1u_xf_1ph/*.dss
} > "$SECONDARY_OUT/metadata.txt"
git diff -- problems/msopf/root_level > "$SECONDARY_OUT/source.patch"
t=$SECONDS; code=0
julia --startup-file=no --project=problems/msopf problems/msopf/root_level/matched_secondary_study.jl "$label" "$phase" > "$SECONDARY_OUT/run.log" 2> "$SECONDARY_OUT/run.err" || code=$?
echo "exit_code=$code" >> "$SECONDARY_OUT/metadata.txt"
echo "wall_seconds=$((SECONDS-t))" >> "$SECONDARY_OUT/metadata.txt"
echo "$label $phase exit=$code wall=$((SECONDS-t)) s"
exit "$code"
