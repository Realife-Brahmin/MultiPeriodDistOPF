#!/usr/bin/env bash
set -euo pipefail
C=$PWD
OLD=/home/aryan.jha/kamiak_jobs/value_coupling_audit_oct22
BASE="$OLD/source"
test ! -e "$C/source"
mkdir -p "$C/source/spar_ddp/examples/power_system" "$C/source/spar_ddp/results/network_filterddp" "$C/logs" "$C/results"
cp "$BASE"/spar_ddp/examples/power_system/*.jl "$C/source/spar_ddp/examples/power_system/"
cp "$C/capture_helper.jl" "$C/source/spar_ddp/examples/power_system/frozen_value_capture.jl"
python3 - "$C/source/spar_ddp/examples/power_system/procs_pass.jl" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]);s=p.read_text()
needle='fc_before(M,nx,d.k,d.j,μ)'
assert s.count(needle)==1
p.write_text(s.replace(needle,'fc_before(M,nx,d.k,d.j,μ,d.reg_last)'))
PY
ln -s "$BASE/envs" "$C/source/envs"
ln -s "$BASE/rawData" "$C/source/rawData"
ln -s "$BASE/spar_ddp/DDP4OPF.jl" "$C/source/spar_ddp/DDP4OPF.jl"
ln -s "$BASE/spar_ddp/experiments" "$C/source/spar_ddp/experiments"
ln -s "$BASE/spar_ddp/results/network_filterddp/network_data_large10kC_1ph_T24_periodic.jls" "$C/source/spar_ddp/results/network_filterddp/network_data_large10kC_1ph_T24_periodic.jls"
find "$C/source/spar_ddp/examples/power_system" -maxdepth 1 -type f -name '*.jl' -print0 | sort -z | xargs -0 sha256sum > "$C/runtime.sha256"
sha256sum "$C"/*.jl "$C"/*.sh > "$C/source.sha256"
sha256sum "$C/source/spar_ddp/results/network_filterddp/"network_data* > "$C/input.sha256"
df -h "$C" > "$C/disk-before.txt"
[[ $(df --output=avail -B1 "$C" | tail -1) -ge 32212254720 ]] || { echo 'DISK_GUARD less than30GiB';exit 2; }
bash -n "$C/capture_job.sh"

