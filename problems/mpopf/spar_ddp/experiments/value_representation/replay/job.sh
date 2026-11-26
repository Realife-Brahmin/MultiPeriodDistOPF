#!/usr/bin/env bash
set -euo pipefail
S=$(realpath "$1"); cd "$S/source"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
J=/home/aryan.jha/opt/julia-1.12.7/bin/julia
test -x "$J"
{ hostname; "$J" --version; lscpu; taskset -pc $$; scontrol show job "$SLURM_JOB_ID"; } > "$S/hardware-$SLURM_JOB_ID.txt"
grep -q 'Gold 6438Y+' "$S/hardware-$SLURM_JOB_ID.txt"
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OCCUPANCY_GUARD other_jobs=$OTHER"; exit 2; }
sha256sum -c "$S/source.sha256"
timeout --signal=TERM --kill-after=30s 10m "$J" --startup-file=no spar_ddp/experiments/value_representation/replay/replay.jl /home/aryan.jha/kamiak_jobs/stage_b_replay_oct15/capture /home/aryan.jha/kamiak_jobs/stage_b_replay_oct15/repo/ddp/results/network_filterddp/network_data_large10kC_1ph_T96_periodic.jls "$S/results" 8
date -u +%FT%TZ > "$S/COMPLETE"
