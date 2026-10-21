#!/usr/bin/env bash
set -euo pipefail
cd /home/aryan.jha/kamiak_jobs/balanced_stage_b_oct15
bash -n run_replay.sh
sha256sum replay.jl run_replay.sh > source.sha256
sbatch --parsable --partition=kamiak --account=dubey --job-name=stageB_balanced --nodes=1 --ntasks=64 --cpus-per-task=1 --mem=600G --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --time=00:12:00 --output=logs/balanced_%j.out --error=logs/balanced_%j.err run_replay.sh 0
