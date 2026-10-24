#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1")
BASE=/home/aryan.jha/kamiak_jobs/value_coupling_audit_oct22/source
test ! -e "$C/source"
mkdir -p "$C/source/spar_ddp/examples/power_system" "$C/source/spar_ddp/results/network_filterddp" "$C/logs" "$C/results"
# Only1.1MB of drivers copied. Everything large/reused remains immutable.
cp "$BASE"/spar_ddp/examples/power_system/*.jl "$C/source/spar_ddp/examples/power_system/"
cp "$C/campaign_capture_helper.jl" "$C/source/spar_ddp/examples/power_system/frozen_value_capture.jl"
ln -s "$BASE/envs" "$C/source/envs"
ln -s "$BASE/rawData" "$C/source/rawData"
ln -s "$BASE/spar_ddp/DDP4OPF.jl" "$C/source/spar_ddp/DDP4OPF.jl"
ln -s "$BASE/spar_ddp/experiments" "$C/source/spar_ddp/experiments"
for spec in ieee123C_1ph_T3 large10kC_1ph_T24;do
    ln -s "$BASE/spar_ddp/results/network_filterddp/network_data_${spec}_periodic.jls" "$C/source/spar_ddp/results/network_filterddp/network_data_${spec}_periodic.jls"
done
find "$C/source/spar_ddp/examples/power_system" -maxdepth 1 -type f -name '*.jl' -print0 | sort -z | xargs -0 sha256sum > "$C/runtime.sha256"
sha256sum "$C"/*.jl "$C"/*.sh > "$C/source.sha256"
sha256sum "$C/source/spar_ddp/results/network_filterddp/"network_data* > "$C/input.sha256"
df -h "$C" > "$C/disk-before.txt"
[[ $(df --output=avail -B1 "$C" | tail -1) -ge 2147483648 ]] || { echo 'DISK_GUARD less than2GiB';exit 2; }
