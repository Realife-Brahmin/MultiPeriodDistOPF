#!/usr/bin/env bash
set -euo pipefail
T=$1; [[ $T == 144 || $T == 288 ]] || exit 2
C=$PWD; BASE=/home/aryan.jha/kamiak_jobs/battery_count_oct21_campaign
# prepare.sh copies ONLY the completed original immutable fast-kernel snapshot.
# It does not copy or edit the currently running v4 source.
bash prepare.sh
S="$C/variants/10k_1020/source"
F=network_data_large10kC_1ph_T${T}_periodic.jls
test -s "$BASE/inputs/$F"
if [[ ! -e "$S/spar_ddp/results/network_filterddp/$F" ]]; then
 ln -s "$BASE/inputs/$F" "$S/spar_ddp/results/network_filterddp/$F"
fi
printf '%s\n' "$T" > case.txt
(cd "$S"; find envs rawData spar_ddp \( -type f -o -type l \) -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256")
python3 compare_source_manifests.py /home/aryan.jha/kamiak_jobs/barrier_reuse_fullsolve_nov06_v4/execution-source.sha256 execution-source.sha256 "$T" > common-source-match.json
sha256sum *.sh *.py *.jl *.md case.txt quality_reference.json reused_reference_rows.csv common-source-match.json > harness.sha256
