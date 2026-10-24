#!/usr/bin/env bash
# MSOPF job on Kamiak. Usage (from the job folder):
#   sbatch [slurm options] kamiak_msopf_job.sh <label> <julia script and args ...>
# Environment it reads: SOLVER=mumps|ma57 (default mumps), and anything the driver reads
# (LARGE_AT, LARGE_T, ...). Writes its own output folder out_<label>_<jobid>/ and logs/.
set -uo pipefail
label=$1; shift
JOBDIR=$(pwd)
REPO=$HOME/msopf/repo
JULIA=$HOME/opt/julia-1.12.7/bin/julia
HSL=$HOME/opt/hsl-private/install/ma57-3.11.3/lib
OUT=$JOBDIR/out_${label}_${SLURM_JOB_ID:-local}
mkdir -p "$OUT" "$JOBDIR/logs"
export MULTIPERIOD_OUT_DIR=$OUT
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
unset NO_COLOR FORCE_COLOR
if [ "${SOLVER:-mumps}" = ma57 ]; then
  export LD_PRELOAD=$HSL/libfakemetis.so
  export IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=$HSL/libma57.so;linear_system_scaling=none;ma57_pivot_order=2"
fi
{
  echo "label=$label job=${SLURM_JOB_ID:-local} host=$(hostname) start=$(date -Is)"
  echo "solver=${SOLVER:-mumps} IPOPT_EXTRA_OPTIONS=${IPOPT_EXTRA_OPTIONS:-}"
  echo "args=$*"
  env | grep -E '^(LARGE_|MULTIPERIOD_)' | sort
  (cd "$REPO" && echo "commit=$(git rev-parse HEAD)" && git status --porcelain | head)
  lscpu | grep -E 'Model name|^CPU\(s\)'
  [ -n "${SLURM_JOB_ID:-}" ] && scontrol show job "$SLURM_JOB_ID" | grep -E 'NodeList|TRES|TimeLimit|Partition'
} > "$OUT/metadata.txt" 2>&1
cd "$REPO"
t0=$(date +%s)
"$JULIA" --startup-file=no --project=problems/msopf "$@" > "$OUT/run.log" 2> "$OUT/run.err"
code=$?
echo "exit=$code end=$(date -Is)" >> "$OUT/metadata.txt"
echo "wall_s=$(( $(date +%s) - t0 ))" >> "$OUT/metadata.txt"
[ -n "${SLURM_JOB_ID:-}" ] && sstat -n -j "$SLURM_JOB_ID.batch" --format=MaxRSS,MaxVMSize >> "$OUT/metadata.txt" 2>&1
exit $code
