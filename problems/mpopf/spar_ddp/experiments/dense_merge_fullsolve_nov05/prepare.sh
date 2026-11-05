#!/usr/bin/env bash
set -euo pipefail
C=$PWD
BASE=/home/aryan.jha/kamiak_jobs/battery_count_oct21_campaign
FROM="$BASE/source"
S="$C/variants/10k_1020/source"
sha256sum -c "$BASE/original-input.sha256" > inherited-input-check.log
mkdir -p "$S/spar_ddp"
cp -a "$FROM/envs" "$FROM/rawData" "$S/"
for p in "$FROM"/spar_ddp/*; do
    [[ ${p##*/} == results ]] || cp -a "$p" "$S/spar_ddp/"
done
PS="$S/spar_ddp/examples/power_system"
cp right_kernel.jl "$PS/study_right_kernel.jl"
cp schedule.jl "$PS/study_schedule.jl"
cp control.jl "$PS/study_control.jl"
python3 patch.py "$S"
mkdir -p "$S/spar_ddp/results/network_filterddp"
for T in 3 96 288; do
    F="network_data_large10kC_1ph_T${T}_periodic.jls"
    ln -s "$BASE/inputs/$F" "$S/spar_ddp/results/network_filterddp/$F"
done
mkdir -p results jobs
(cd "$S"; find envs rawData spar_ddp \( -type f -o -type l \) -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256")
sha256sum *.jl *.py *.sh > harness.sha256
