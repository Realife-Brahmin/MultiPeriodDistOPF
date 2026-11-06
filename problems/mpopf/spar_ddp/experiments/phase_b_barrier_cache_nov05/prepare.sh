#!/usr/bin/env bash
# Run ONLY in a new immutable snapshot containing this experiment's source.
set -euo pipefail
C=$PWD
BASE=/home/aryan.jha/kamiak_jobs/phase_b_influence_replay_nov05
HIST=/home/aryan.jha/kamiak_jobs/dense_merge_T24_nov05
[[ ! -e "$C/results" && ! -e "$C/lib" ]]
mkdir "$C/lib"
cp "$BASE/lib/"{replay,production_reference,hierarchy}.jl "$C/lib/"
cp /home/aryan.jha/kamiak_jobs/dense_merge_core128_nov05/affinity.py "$C/affinity.py"
mkdir "$C/frequency"
for rep in 1 2;do
 mkdir "$C/frequency/r$rep"
 D="$HIST/results/10k_1020_T24_ddp_roriginal_$rep"
 cp "$D/"{iteration_timing.csv,passes.csv,settings.txt,COMPLETE,timing_summary.json} "$C/frequency/r$rep/"
done
cp "$HIST/execution-source.sha256" "$C/frequency/"
bash -n "$C/job.sh"
sha256sum *.jl *.sh *.py lib/*.jl > source.sha256
