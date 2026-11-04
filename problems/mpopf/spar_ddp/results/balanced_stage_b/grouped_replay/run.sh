#!/usr/bin/env bash
set -euo pipefail
export REPLAY_ROOT="$SLURM_SUBMIT_DIR"
export REPLAY_REPO=/home/aryan.jha/kamiak_jobs/stage_b_replay_oct15/repo
export REPLAY_CAPTURE=/home/aryan.jha/kamiak_jobs/stage_b_replay_oct15/capture
export PATH="$HOME/bin:$PATH"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE
export FILTERDDP_TREE_KKT=1 REPLAY_MULTINODE=0
date -Is
lscpu
scontrol show job "$SLURM_JOB_ID"
node=$(hostname -s)
others=$(squeue -h -w "$node" -t RUNNING -o '%A' | sort -u | grep -v "^${SLURM_JOB_ID}$" || true)
test -z "$others" || { echo "Other active jobs on $node: $others"; exit 3; }
sha256sum -c source.sha256
sha256sum -c capture.sha256
sha256sum "$REPLAY_REPO"/ddp/examples/power_system/{ieee123c_filterddp,procs_pass}.jl
julia --version
julia --startup-file=no "$REPLAY_ROOT/check_schedule.jl"
cd "$REPLAY_REPO"
julia --startup-file=no --project=envs/ddp2026 "$REPLAY_ROOT/replay.jl"
cgroup=$(awk -F: '$1=="0"{print $3}' /proc/self/cgroup)
if test -n "$cgroup" && test -r "/sys/fs/cgroup${cgroup}/memory.peak"; then
    printf 'BATCH_CGROUP_PEAK_BYTES '
    cat "/sys/fs/cgroup${cgroup}/memory.peak"
fi
