#!/usr/bin/env bash
set -euo pipefail
C=$PWD
OLD=/home/aryan.jha/kamiak_jobs/dense_merge_fullsolve_nov05
BASE=/home/aryan.jha/kamiak_jobs/battery_count_oct21_campaign
test -s "$OLD/T288.finished"
(cd "$OLD/variants/10k_1020/source"; sha256sum -c "$OLD/execution-source.sha256") > inherited-source-check.log
sha256sum -c "$BASE/original-input.sha256" > inherited-input-check.log
cp -a "$OLD/variants" "$C/"
for file in job.sh affinity.py compare_solutions.jl right_kernel.jl schedule.jl control.jl patch.py; do cp "$OLD/$file" "$C/$file"; done
mkdir posthoc_v2
cp "$OLD/posthoc_v2/compare_solutions.jl" posthoc_v2/
cp "$OLD/T96.validated" "$C/T96.validated"
printf '%s\n' c4 > selected_candidate.txt
printf '%s\n' 0a5f17114d058974004ea037dab207664ebe1136 > parent.commit
S="$C/variants/10k_1020/source"
F=network_data_large10kC_1ph_T144_periodic.jls
test ! -e "$S/spar_ddp/results/network_filterddp/$F"
ln -s "$BASE/inputs/$F" "$S/spar_ddp/results/network_filterddp/$F"
mkdir jobs results
(cd "$S"; find envs rawData spar_ddp \( -type f -o -type l \) -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256")
sha256sum *.jl *.py *.sh > harness.sha256
sha256sum posthoc_v2/compare_solutions.jl > posthoc_v2/source.sha256
