#!/usr/bin/env bash
set -euo pipefail
cd /home/aryan.jha/kamiak_jobs/stage_b_replay_oct15
capture=$1
sha256sum replay.jl capture.sh run_replay.sh repo/ddp/examples/power_system/procs_pass.jl > replay_source.sha256
one=$(sbatch --parsable --dependency=afterok:$capture --partition=kamiak --account=dubey --job-name=stageB_one --nodes=1 --ntasks=64 --cpus-per-task=1 --mem=600G --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --time=00:20:00 --output=logs/one_%j.out --error=logs/one_%j.err run_replay.sh 0)
two=$(sbatch --parsable --dependency=afterok:$one --partition=kamiak --account=dubey --job-name=stageB_two --nodes=2 --ntasks=128 --cpus-per-task=1 --mem=600G --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --time=00:20:00 --output=logs/two_%j.out --error=logs/two_%j.err run_replay.sh 1)
printf 'capture=%s\none=%s\ntwo=%s\n' "$capture" "$one" "$two" | tee jobs.txt
