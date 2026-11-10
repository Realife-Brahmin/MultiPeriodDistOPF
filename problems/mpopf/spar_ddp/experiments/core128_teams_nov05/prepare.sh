#!/usr/bin/env bash
set -euo pipefail
C=$PWD
OLD=/home/aryan.jha/kamiak_jobs/dense_merge_core128_nov05
BASE=/home/aryan.jha/kamiak_jobs/battery_count_oct21_campaign
test -s "$OLD/core128.finished"
(cd "$OLD/variants/10k_1020/source"; sha256sum -c "$OLD/execution-source.sha256") > inherited-source-check.log
sha256sum -c "$BASE/original-input.sha256" > inherited-input-check.log
cp -a "$OLD/variants" .
for f in affinity.py compare_solutions.jl execution-source.sha256; do cp "$OLD/$f" .; done
cp "$OLD/source.commit" inherited.commit
mkdir jobs results
python3 test_affinity.py affinity.py > affinity-test.log
bash -n job.sh
sha256sum *.sh *.py *.jl > harness.sha256

