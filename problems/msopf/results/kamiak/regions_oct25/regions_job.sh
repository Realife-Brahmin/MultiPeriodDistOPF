#!/usr/bin/env bash
# One task of a region campaign on Kamiak (array job). Usage (from the job folder):
#   sbatch --array=1-N [slurm options] regions_job.sh tasks.txt
# Each line of tasks.txt is "<stage> <case> <directions>", e.g. "four_smartds_p1u_typeiii cyclic 1:2" or
# "three_ieee9500xf nobatt 0:1"; directions are colon-separated. Four-source stages use the 78 Fibonacci
# directions, three-source stages the 16 of stage "three". Settings as in the paper's campaigns: T = 24,
# MA57, import margin with voltage room (MULTI_VMARGIN=20000), the reference solve uncapped and later solves
# capped at 8x its time (retry factor MULTIPERIOD_RETRY, default 4). Writes out_<stage>_<case>_d<dirs>_<jobid>/.
set -uo pipefail
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$1" | tr -d '\r')
read -r stage case dirs <<< "$line"
JOBDIR=$(pwd)
REPO=$HOME/msopf/repo
JULIA=$HOME/opt/julia-1.12.7/bin/julia
HSL=$HOME/opt/hsl-private/install/ma57-3.11.3/lib
OUT=$JOBDIR/out_${stage}_${case}_d${dirs//:/_}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$OUT"
export MULTIPERIOD_OUT_DIR=$OUT LARGE_T=24 MULTIPERIOD_LEAN=1 MULTI_VMARGIN=20000
export MULTI_DIRSET=fib78 MULTI_DIRS=$dirs MULTI_CASES=$case MULTI_SAMPLES=0 THREE_DIRS=$dirs THREE_CASES=$case
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
export LD_PRELOAD=$HSL/libfakemetis.so
export IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=$HSL/libma57.so;linear_system_scaling=none;ma57_pivot_order=2"
unset NO_COLOR FORCE_COLOR
{
  echo "stage=$stage case=$case dirs=$dirs job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} host=$(hostname) start=$(date -Is)"
  env | grep -E '^(LARGE_|MULTIPERIOD_|MULTI_|THREE_|IPOPT_)' | sort
  (cd "$REPO" && echo "commit=$(git rev-parse HEAD)" && git status --porcelain | head)
  lscpu | grep -E 'Model name'
} > "$OUT/metadata.txt" 2>&1
cd "$REPO"
t0=$(date +%s)
"$JULIA" --startup-file=no --project=problems/msopf problems/msopf/root_level/multiperiod_study.jl "$stage" > "$OUT/run.log" 2> "$OUT/run.err"
code=$?
echo "exit=$code end=$(date -Is) wall_s=$(( $(date +%s) - t0 ))" >> "$OUT/metadata.txt"
exit $code
