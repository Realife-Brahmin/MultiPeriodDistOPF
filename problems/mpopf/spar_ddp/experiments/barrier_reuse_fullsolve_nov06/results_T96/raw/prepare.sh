#!/usr/bin/env bash
set -euo pipefail
C=$PWD
OLD=/home/aryan.jha/kamiak_jobs/dense_merge_teams96_nov05
BASE=/home/aryan.jha/kamiak_jobs/battery_count_oct21_campaign
test -s "$OLD/T96.finished"
(cd "$OLD/variants/10k_1020/source"; sha256sum -c "$OLD/execution-source.sha256") > inherited-source-check.log
sha256sum -c "$BASE/original-input.sha256" > inherited-input-check.log
cp -a "$OLD/variants" .
cp "$OLD/affinity.py" .
cp "$OLD/compare_solutions.jl" .
mkdir jobs results
S="$C/variants/10k_1020/source"
F=network_data_large10kC_1ph_T24_periodic.jls
if [[ ! -e "$S/spar_ddp/results/network_filterddp/$F" ]]; then
 ln -s "$BASE/inputs/$F" "$S/spar_ddp/results/network_filterddp/$F"
fi
cp cache_kernel.jl "$S/spar_ddp/examples/power_system/br_cache_kernel.jl"
cp integration.jl "$S/spar_ddp/examples/power_system/br_integration.jl"
python3 patch.py "$S"
(cd "$S"; find envs rawData spar_ddp \( -type f -o -type l \) -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256")
sha256sum *.sh *.py *.jl *.md > harness.sha256
