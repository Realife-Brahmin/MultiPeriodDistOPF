#!/usr/bin/env bash
# Certificate pilot on Kamiak (certificate_secondary_pilot.jl), reusing the committed hour boundaries and envelope.
# Usage (from a job folder): sbatch [slurm options] certificate_job.sh <qdldl|ma57> [floor] [cuts]
# CERT_WINDOW=1 (exported) runs the certified-window search instead, one direction per array task:
# CERT_THETA = 22.5 x (SLURM_ARRAY_TASK_ID - 1).
# Environment ~/msopf/envs/sdp_hsl: problems/msopf/environments/sdp plus HSL.jl; MA57 comes from
# ~/opt/hsl-private/julia_libhsl (a libhsl wrapper around the private MA57/MA97 builds with METIS).
set -uo pipefail
ls=$1; floor=${2:-1e-4}; cuts=${3:-kcl0,jbound}
REPO=$HOME/msopf/repo
tag=${SLURM_ARRAY_TASK_ID:+_d$SLURM_ARRAY_TASK_ID}
OUT=$(pwd)/out_${ls}${tag}_${SLURM_JOB_ID:-local}
[ -n "${SLURM_ARRAY_TASK_ID:-}" ] && export CERT_THETA=$(awk -v i=$SLURM_ARRAY_TASK_ID 'BEGIN {printf "%.1f", 22.5 * (i - 1)}')
mkdir -p "$OUT"
cp "$REPO"/problems/msopf/results/multiperiod/certificate_secondary_{hours,rstar,envelope}.csv "$OUT"/
export MULTIPERIOD_OUT_DIR=$OUT MULTIPERIOD_SOLVE_LOG=0 CERT_REUSE=1 CERT_FLOOR=$floor CERT_CUTS=$cuts CERT_LINSOLVE=$ls CERT_TAG=_$ls
export JULIA_HSL_LIBRARY_PATH=$HOME/opt/hsl-private/julia_libhsl OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
unset NO_COLOR FORCE_COLOR
{ echo "linsolve=$ls floor=$floor cuts=$cuts job=${SLURM_JOB_ID:-local} host=$(hostname) start=$(date -Is)"
  (cd "$REPO" && echo "commit=$(git rev-parse HEAD)"); } > "$OUT/metadata.txt"
cd "$REPO"
$HOME/opt/julia-1.12.7/bin/julia --startup-file=no --project=$HOME/msopf/envs/sdp_hsl problems/msopf/root_level/certificate_secondary_pilot.jl > "$OUT/run.log" 2> "$OUT/run.err"
code=$?
echo "exit=$code end=$(date -Is)" >> "$OUT/metadata.txt"
exit $code
