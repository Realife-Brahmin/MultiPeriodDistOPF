#!/usr/bin/env bash
set -euo pipefail
# Add an exact, private runtime without changing ~/bin/julia or other campaigns.
cd /home/aryan.jha/opt
if [[ ! -x julia-1.12.4/bin/julia ]]; then
 curl --fail --location --retry 2 -o julia-1.12.4-linux-x86_64.tar.gz https://julialang-s3.julialang.org/bin/linux/x64/1.12/julia-1.12.4-linux-x86_64.tar.gz
 curl --fail --location --retry 2 -o julia-1.12.4.sha256 https://julialang-s3.julialang.org/bin/checksums/julia-1.12.4.sha256
 grep 'julia-1.12.4-linux-x86_64.tar.gz$' julia-1.12.4.sha256 | sha256sum -c -
 tar --no-same-owner -xzf julia-1.12.4-linux-x86_64.tar.gz
fi
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export TREE_ABLATION_EXPORT_SMOKE=1
cd /home/aryan.jha/kamiak_jobs/tree_ablation_oct20_checkout
bash spar_ddp/experiments/tree_ablation/submit.sh /home/aryan.jha/kamiak_jobs/proper_logging_campaign_oct16/repo/ddp/results/network_filterddp /home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign
