#!/usr/bin/env bash
# Exploratory (not a campaign row): large10k, 63 worker processes, full warm-up,
# collector once a pass, OpenDSS check after every solve. Arms per T:
#   rounds   dense value functions, pairwise rounds over all periods (current code)
#   chunked  dense, FILTERDDP_PROCS_CHUNKED=1 (chunk reduce, rounds over 63 summaries, close)
#   blocks   FILTERDDP_VALUE_BLOCKS=1 (feeder blocks; exact on large10k only: 102 feeders at the slack bus)
# Submitted from ~/kamiak_jobs/block_value_nov17 with:
#   sbatch --partition=kamiak --account=dubey --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=300G \
#     --time=05:00:00 --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' \
#     --job-name=block_value --output=slurm-%j.out job.sh "<T:arm list>"
set -uo pipefail
C=$SLURM_SUBMIT_DIR; cd "$C/repo" || exit 2
CASES=${1:-"144:chunked 144:rounds 288:chunked 288:rounds"}
export PATH=$HOME/opt/julia-1.12.7/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:/home/aryan.jha/kamiak_jobs/dense_merge_core128_nov05/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE
unset NO_COLOR FORCE_COLOR OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
M="$C/jobs/$SLURM_JOB_ID"; mkdir -p "$M"
{ date -Is; hostname; julia --version; lscpu; scontrol show job "$SLURM_JOB_ID"; } > "$M/hardware.txt"
grep -q '6438Y+' "$M/hardware.txt" || { echo "WRONG_NODE"; exit 2; }
squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%i %u %C' > "$M/node_jobs.txt"
P=problems/mpopf/spar_ddp
mkdir -p $P/results/network_filterddp $P/results/parallel_in_time/logs $P/results/opendss_check $P/results/kkt_ordering/captures/solutions
if [ ! -f "$C/env.ok" ]; then
  julia --startup-file=no --project=$P -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' > "$M/env.log" 2>&1 || { echo ENV_FAILED; exit 3; }
  touch "$C/env.ok"
fi
export FILTERDDP_PROCS_WORKERS=63 FILTERDDP_PROCS_GC=pass PAIRED_WARMUP=full PAIRED_TOL=1e-5
for c in $CASES; do
  T=${c%%:*}; ARM=${c##*:}
  case $ARM in rounds) BV=0 CH=0;; chunked) BV=0 CH=1;; blocks) BV=1 CH=0;; *) echo "bad arm $ARM"; continue;; esac
  D=$P/results/network_filterddp/network_data_large10kC_1ph_T${T}_periodic.jls
  if [ ! -s "$D" ]; then
    PROFILE_PERIODIC=1 julia --startup-file=no --project=$P $P/examples/power_system/export_ieee123c_data.jl large10kC_1ph $T > "$M/export_T$T.log" 2>&1 || { echo "EXPORT_FAILED T=$T" >> "$M/summary.txt"; continue; }
    sha256sum "$D" >> "$M/inputs.sha256"
  fi
  TAG=procs_${ARM}_w63_j$SLURM_JOB_ID
  echo "BEGIN T=$T arm=$ARM $(date -Is)" >> "$M/summary.txt"
  t0=$(date +%s)
  FILTERDDP_VALUE_BLOCKS=$BV FILTERDDP_PROCS_CHUNKED=$CH bash $P/examples/power_system/run_paired.sh large10kC_1ph $T procs $TAG > "$M/run_T${T}_$ARM.out" 2>&1
  L=$P/results/parallel_in_time/logs/fddp_large10kC_1ph_T${T}_paired_${TAG}.log
  cp "$L" "$M/"; cp $P/results/opendss_check/parsim_large10kC_1ph_T${T}_paired_${TAG}.log "$M/" 2>/dev/null
  echo "END T=$T arm=$ARM wall=$(( $(date +%s) - t0 ))s | $(grep -E 'solve complete' "$L" | tail -1) | $(grep -oE 'FilterDDP objective=[-0-9.eE+]+ max_equality_residual=[0-9.eE+-]+' "$L") | $(grep -oE 'below limit: [0-9]+ of [0-9]+|above limit: [0-9]+' "$M/parsim_large10kC_1ph_T${T}_paired_${TAG}.log" 2>/dev/null | tr '\n' ' ')" >> "$M/summary.txt"
done
date -Is > "$M/finished"
