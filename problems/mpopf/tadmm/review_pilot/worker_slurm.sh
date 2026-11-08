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
exec taskset -c "$TADMM_PIN_CPU" "${TADMM_JULIA_EXE:?}" "$@"
