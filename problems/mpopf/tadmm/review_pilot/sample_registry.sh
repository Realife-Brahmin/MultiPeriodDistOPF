#!/usr/bin/env bash
set -euo pipefail
CPU_RANGE=$(awk '/Cpus_allowed_list:/ {print $2}' /proc/self/status)
PIN=$(python3 -c 'import sys
a=[]
for p in sys.argv[1].split(","):
 b=p.split("-");a.extend(range(int(b[0]),int(b[-1])+1))
print(a[-1])' "$CPU_RANGE")
exec taskset -c "$PIN" python3 "$(dirname "$0")/sample_registry.py" "$1"
