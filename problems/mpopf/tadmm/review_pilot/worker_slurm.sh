#!/usr/bin/env bash
set -euo pipefail
# srun launches one process per period. First/last allocated CPUs are reserved
# for the coordinator and sampler on each host; workers occupy the middle.
CPU_RANGE=$(awk '/Cpus_allowed_list:/ {print $2}' /proc/self/status)
export TADMM_PIN_CPU=$(python3 -c 'import sys
a=[]
for p in sys.argv[1].split(","):
 b=p.split("-");a.extend(range(int(b[0]),int(b[-1])+1))
i=int(sys.argv[2])+1
assert i<len(a)-1, (i,a)
print(a[i])' "$CPU_RANGE" "$SLURM_LOCALID")
export TADMM_MONITOR_CPU=$(python3 -c 'import sys
p=sys.argv[1].split(",")[-1];print(p.split("-")[-1])' "$CPU_RANGE")
# The first task on each node also starts that node's memory sampler on the spare last CPU
# (Kamiak's srun has no --overlap, so the sampler cannot be a step of its own).
if [[ -n "${SC_NODE_SAMPLE_DIR:-}" && "${SLURM_LOCALID:-1}" == 0 ]]; then
  taskset -c "$TADMM_MONITOR_CPU" python3 "$(dirname "$0")/sample_memory.py" sample-node \
    "$SC_NODE_SAMPLE_DIR/memory_node_$(hostname -s).csv" "${SC_SAMPLE_S:-2}" "$SC_NODE_SAMPLE_DIR/STOP_SAMPLING" &
fi
exec taskset -c "$TADMM_PIN_CPU" "${TADMM_JULIA_EXE:?}" "$@"
