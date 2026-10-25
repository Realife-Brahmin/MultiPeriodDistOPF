#!/usr/bin/env bash
# Pinned, isolated matched-secondary campaign. No DDP checkout is used.
set -euo pipefail
ROOT="$HOME/msopf/jobs/matched_secondary_oct25"
REPO="$ROOT/reference_repo"
JULIA="$HOME/opt/julia-1.12.7/bin/julia"
HSL="$HOME/opt/hsl-private/install/ma57-3.11.3/lib"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1 MULTIPERIOD_LEAN=1
export LD_PRELOAD="$HSL/libfakemetis.so"
export IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=$HSL/libma57.so;linear_system_scaling=none;ma57_pivot_order=2"
export MULTI_VMARGIN=20000 MULTIPERIOD_MAX_WALL=360 MULTIPERIOD_SOLVE_LOG=1 MULTIPERIOD_FAILED_SNAPSHOT=0
export SECONDARY_T=24 SECONDARY_LIMIT=60 MATCHED_ROOT="$ROOT/results"
label="$1";phase="reference"
if [[ "$phase" == array ]]; then
 phase=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$ROOT/tasks_${label}.txt")
 [[ -n "$phase" ]] || { echo 'Missing task mapping' >&2;exit 2; }
fi
folder="${phase//,/x}"
export SECONDARY_OUT="$MATCHED_ROOT/$label/$folder" MULTIPERIOD_OUT_DIR="$MATCHED_ROOT/$label/$folder"
mkdir -p "$SECONDARY_OUT"
[[ ! -e "$SECONDARY_OUT/metadata.txt" ]] || { echo 'Refusing to overwrite evidence';exit 2; }
cd "$REPO"
{
 date -u '+start_utc=%FT%TZ';echo "label=$label phase=$phase job=${SLURM_JOB_ID:-none}";hostname
 git rev-parse HEAD;git status --short
 env | sort | grep -E '^(MATCHED_|SECONDARY_|JULIA_|OPENBLAS|OMP_|MULTI|IPOPT_|SLURM_CPUS|SLURM_MEM|SLURM_JOB_ID)'
 sha256sum problems/msopf/{Project.toml,Manifest.toml,full_angle_pf.jl,multiperiod_full_angle.jl} problems/msopf/root_level/{matched_secondary_study.jl,secondary_8500_study.jl,multiperiod_study.jl}
 sha256sum test_networks/ieee9500{,xf}_1ph/*.dss test_networks/smartds_p1u_xf_1ph/*.dss
} > "$SECONDARY_OUT/metadata.txt"
t=$SECONDS;code=0
"$JULIA" --startup-file=no --project=problems/msopf problems/msopf/root_level/matched_banded_reference.jl "$label" > "$SECONDARY_OUT/run.log" 2> "$SECONDARY_OUT/run.err" || code=$?

echo "exit_code=$code" >> "$SECONDARY_OUT/metadata.txt"
echo "wall_seconds=$((SECONDS-t))" >> "$SECONDARY_OUT/metadata.txt"
echo "$label $phase exit=$code"
exit "$code"
