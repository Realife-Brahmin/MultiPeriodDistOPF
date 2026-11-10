#!/usr/bin/env bash
set -euo pipefail
C=$PWD
OLD=/home/aryan.jha/kamiak_jobs/dense_merge_core128_nov05
test -s "$OLD/core128.finished"
(cd "$OLD/variants/10k_1020/source"; sha256sum -c "$OLD/execution-source.sha256") > inherited-source-check.log
cp -a "$OLD/variants" .
for f in affinity.py affinity_legacy.py compare_solutions.jl check_ipopt_logs.py test_affinity.py execution-source.sha256; do cp "$OLD/$f" .; done
mkdir jobs results
cp "$OLD/source.commit" inherited.commit
sha256sum *.sh *.py *.jl > harness.sha256
bash -n job.sh
python3 test_affinity.py affinity.py > affinity-test.log

