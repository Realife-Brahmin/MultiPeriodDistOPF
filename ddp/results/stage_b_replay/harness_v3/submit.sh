#!/usr/bin/env bash
set -euo pipefail
cd /home/aryan.jha/kamiak_jobs/stage_b_replay_oct15_v3
chmod +x bin/srun
sha256sum replay.jl run_replay.sh bin/srun > replay_source.sha256
sbatch --parsable --partition=kamiak --account=dubey --job-name=stageB_two_v3 --nodes=2 --ntasks=128 --cpus-per-task=1 --mem=600G --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --time=00:10:00 --output=logs/two_%j.out --error=logs/two_%j.err run_replay.sh 1
