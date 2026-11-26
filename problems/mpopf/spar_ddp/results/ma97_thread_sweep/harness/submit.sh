#!/usr/bin/env bash
set -euo pipefail
cd /home/aryan.jha/kamiak_jobs/ma97_thread_sweep_oct16
bash -n run_case.sh
sha256sum run_case.sh thread_entry.jl check_ipopt_logs.py > source.sha256
grep '  repo/' parent_source.sha256 >> source.sha256
previous=''
for threads in 1 4 16 64; do
 args=()
 [[ -z "$previous" ]] || args+=(--dependency=afterok:$previous)
 job=$(sbatch --parsable "${args[@]}" --partition=kamiak --account=dubey --job-name="ma97_omp${threads}" --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=300G --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --time=00:30:00 --export=ALL,OMP_SWEEP_THREADS=$threads --output=logs/job_%j.out --error=logs/job_%j.err run_case.sh large10kC_1ph 96 1 ma97 amd 1)
 printf '%s,%s\n' "$job" "$threads" | tee -a submitted_jobs.csv
 previous=$job
done
