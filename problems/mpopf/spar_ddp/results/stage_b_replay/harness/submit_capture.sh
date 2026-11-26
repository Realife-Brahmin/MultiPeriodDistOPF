#!/usr/bin/env bash
set -euo pipefail
cd /home/aryan.jha/kamiak_jobs/stage_b_replay_oct15
bash -n capture.sh
bash -n run_replay.sh
/home/aryan.jha/bin/julia --startup-file=no -e 'Meta.parseall(read(ARGS[1],String)); println("PARSE_OK")' replay.jl
sbatch --parsable --partition=kamiak --account=dubey --job-name=stageB_capture --nodes=1 --ntasks=64 --cpus-per-task=1 --mem=600G --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --time=00:10:00 --output=logs/capture_%j.out --error=logs/capture_%j.err capture.sh
