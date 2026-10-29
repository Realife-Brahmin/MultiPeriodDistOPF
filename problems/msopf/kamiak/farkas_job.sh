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
export MULTIPERIOD_OUT_DIR=$OUT MULTIPERIOD_SOLVE_LOG=0 FB_SYSTEM=$system FB_PHASE=$phase FB_DIR=$dir
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
unset NO_COLOR FORCE_COLOR
log=$OUT/logs/${phase}_${dir}
{ echo "system=$system phase=$phase dir=$dir job=${SLURM_JOB_ID:-local} host=$(hostname) start=$(date -Is)"
  (cd "$REPO" && echo "commit=$(git rev-parse HEAD)"); } > "$log.meta"
cd "$REPO"
$HOME/opt/julia-1.12.7/bin/julia --startup-file=no --project=$HOME/msopf/envs/sdp_hsl problems/msopf/root_level/farkas_bracketing.jl > "$log.log" 2> "$log.err"
code=$?
echo "exit=$code end=$(date -Is)" >> "$log.meta"
exit $code
