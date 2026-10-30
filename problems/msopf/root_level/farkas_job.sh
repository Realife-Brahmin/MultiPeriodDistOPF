#!/usr/bin/env bash
# One task of Farkas bracketing on Kamiak (root_level/farkas_bracketing.jl).
# Usage (from ~/msopf/jobs/<study>): sbatch [--array=1-N] [slurm options] farkas_job.sh <system> <A|B|C>
# All phases of a system share the folder <system>/ (phase B reads phase A's files, phase C reads B's envelope);
# chain them with --dependency. The array index is the direction (FB_DIR). Logs: <system>/logs/<phase>_<dir>.log.
set -uo pipefail
system=$1; phase=$2; dir=${SLURM_ARRAY_TASK_ID:-0}
REPO=$HOME/msopf/repo
OUT=$(pwd)/$system
mkdir -p "$OUT/logs"
export MULTIPERIOD_OUT_DIR=$OUT MULTIPERIOD_SOLVE_LOG=${MULTIPERIOD_SOLVE_LOG:-0} FB_SYSTEM=$system FB_PHASE=$phase FB_DIR=$dir
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
unset NO_COLOR FORCE_COLOR
# Ipopt's linear solver: HSL MA57 by default (FB_LINSOLVE=ma57), as in the DDP runs on Kamiak, from the private
# build in ~/opt/hsl-private (AMD ordering; its dummy METIS is preloaded so the library loads); FB_LINSOLVE=mumps
# keeps Ipopt's own MUMPS, which the first jobs of 2026-10-29 used. ma57_pre_alloc=3 (default 1.05) allocates MA57's
# workspace with room to spare: several sf6937 tasks of 2026-10-29 died in MA57's factorization (signal 11) after 1-2 h.
HSL=$HOME/opt/hsl-private/install/ma57-3.11.3/lib
if [ "${FB_LINSOLVE:-ma57}" = ma57 ]; then
  export LD_PRELOAD=$HSL/libfakemetis.so
  export IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=$HSL/libma57.so;ma57_pivot_order=2;ma57_pre_alloc=3;linear_system_scaling=none"
fi
log=$OUT/logs/${phase}_${dir}
{ echo "system=$system phase=$phase dir=$dir job=${SLURM_JOB_ID:-local} host=$(hostname) start=$(date -Is)"
  echo "ipopt_extra=${IPOPT_EXTRA_OPTIONS:-}"
  (cd "$REPO" && echo "commit=$(git rev-parse HEAD)"); } > "$log.meta"
cd "$REPO"
$HOME/opt/julia-1.12.7/bin/julia --startup-file=no --project=$HOME/msopf/envs/sdp_hsl problems/msopf/root_level/farkas_bracketing.jl > "$log.log" 2> "$log.err"
code=$?
echo "exit=$code end=$(date -Is)" >> "$log.meta"
exit $code
