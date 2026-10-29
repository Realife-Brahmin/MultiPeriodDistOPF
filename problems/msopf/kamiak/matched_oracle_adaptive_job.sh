#!/usr/bin/env bash
# Immutable copy per submitted campaign; reuses but never modifies original caches.
set -euo pipefail
ROOT="$HOME/msopf/jobs/matched_secondary_oct25"
REPO="$ROOT/oracle_adaptive_repo"
JULIA="$HOME/opt/julia-1.12.7/bin/julia"
HSL="$HOME/opt/hsl-private/install/ma57-3.11.3/lib"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MULTIPERIOD_LEAN=1
export LD_PRELOAD="$HSL/libfakemetis.so"
export IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=$HSL/libma57.so;linear_system_scaling=none;ma57_pivot_order=2;ma57_pivtol=1e-4;mu_strategy=adaptive;acceptable_tol=1e-4;acceptable_constr_viol_tol=1e-7;acceptable_iter=1"
export MULTIPERIOD_MAX_WALL=180 MULTIPERIOD_SOLVE_LOG=1 MULTIPERIOD_FAILED_SNAPSHOT=1
export SECONDARY_T=24 SECONDARY_LIMIT=60 MATCHED_ROOT="$ROOT/results"
if [[ "$1" == array ]]; then
 read -r label case direction < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$ROOT/oracle_adaptive_tasks.txt" | tr -d '\r')
else
 label="$1";case="$2";direction="$3"
fi
export SECONDARY_OUT="$ROOT/oracle_adaptive_results/$label/${case}_${direction}" MULTIPERIOD_OUT_DIR="$ROOT/oracle_adaptive_results/$label/${case}_${direction}"
mkdir -p "$SECONDARY_OUT"
[[ ! -e "$SECONDARY_OUT/metadata.txt" ]] || { echo 'Refusing to overwrite evidence';exit 2; }
cd "$REPO"
{
 date -u '+start_utc=%FT%TZ';hostname;git rev-parse HEAD;git status --short
 env | sort | grep -E '^(MATCHED_|SECONDARY_|MULTI|IPOPT_|SLURM_JOB_ID|JULIA_|OPENBLAS|OMP_)'
 sha256sum problems/msopf/{Project.toml,Manifest.toml,multiperiod_full_angle.jl} problems/msopf/root_level/{matched_feasibility_search.jl,matched_secondary_study.jl}
} > "$SECONDARY_OUT/metadata.txt"
t=$SECONDS;code=0
"$JULIA" --startup-file=no --project=problems/msopf problems/msopf/root_level/matched_feasibility_search.jl "$label" "$case" "$direction" > "$SECONDARY_OUT/run.log" 2> "$SECONDARY_OUT/run.err" || code=$?
echo "exit_code=$code" >> "$SECONDARY_OUT/metadata.txt"
echo "wall_seconds=$((SECONDS-t))" >> "$SECONDARY_OUT/metadata.txt"
exit "$code"
