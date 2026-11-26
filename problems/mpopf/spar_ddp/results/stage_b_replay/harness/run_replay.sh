#!/usr/bin/env bash
set -euo pipefail
export REPLAY_ROOT=/home/aryan.jha/kamiak_jobs/stage_b_replay_oct15
cd "$REPLAY_ROOT/repo"
export PATH="$HOME/bin:$PATH"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE
export FILTERDDP_TREE_KKT=1
export REPLAY_MULTINODE=$1
date -Is
lscpu
scontrol show job "$SLURM_JOB_ID"
for node in $(scontrol show hostnames "$SLURM_JOB_NODELIST"); do
    others=$(squeue -h -w "$node" -t RUNNING -o '%A' | sort -u | grep -v "^${SLURM_JOB_ID}$" || true)
    test -z "$others" || { echo "Other active jobs on $node: $others"; exit 3; }
done
sha256sum -c "$REPLAY_ROOT/capture.sha256" > /dev/null
julia --startup-file=no --project=envs/ddp2026 "$REPLAY_ROOT/replay.jl"
