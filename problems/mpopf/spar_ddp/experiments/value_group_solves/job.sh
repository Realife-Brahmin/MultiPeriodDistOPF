#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1");cd "$C/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=2
unset NO_COLOR FORCE_COLOR
{ hostname;julia --version;lscpu;taskset -pc $$;scontrol show job "$SLURM_JOB_ID"; } > "$C/hardware.txt"
grep -q 'Gold 6438Y+' "$C/hardware.txt"
[[ $(julia --version) == 'julia version 1.12.4' ]]
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER";exit 2; }
sha256sum -c "$C/preparation-source.sha256" > "$C/logs/source-check.log"
julia --startup-file=no --project=envs/spar_ddp -e 'using Pkg;Pkg.develop(path="spar_ddp/DDP4OPF.jl");Pkg.instantiate();Pkg.precompile()' > "$C/logs/setup-ddp.log" 2>&1
julia --startup-file=no --project=envs/tadmm -e 'using Pkg;Pkg.instantiate();using OpenDSSDirect' > "$C/logs/setup-opendss.log" 2>&1
julia --startup-file=no spar_ddp/experiments/value_group_solves/unit.jl > "$C/logs/unit.log" 2>&1
# All smoke arms must validate before full cases; never bypass this gate.
for MODE in dense root chunks16;do
 bash spar_ddp/experiments/value_group_solves/run.sh "$C" ieee123C_1ph 3 "$MODE"
done
python3 spar_ddp/experiments/value_group_solves/compare.py "$C" ieee123C_1ph 3 > "$C/logs/smoke-comparison.json"
date -u +%FT%TZ > "$C/smoke.complete"
printf 'system,T,mode,exit_code\n' > "$C/outcomes.csv"
for SYS in ieee123C_1ph ieee2522C_1ph large10kC_1ph;do
 for MODE in dense root chunks16;do
  rc=0;bash spar_ddp/experiments/value_group_solves/run.sh "$C" "$SYS" 24 "$MODE" || rc=$?
  printf '%s,24,%s,%s\n' "$SYS" "$MODE" "$rc" >> "$C/outcomes.csv"
  # A failed dense baseline is a setup/reference blocker, not an approximation result.
  if [[ $MODE == dense && $rc != 0 ]];then exit "$rc";fi
  # A failed approximate solve must not leave a dispatch blocking later arms.
  SOL="spar_ddp/results/network_filterddp/filterddp_solution_${SYS}_T24_periodic_CBsystem.jls"
  if [[ -f "$SOL" ]];then mv "$SOL" "$C/results/${SYS}_T24_${MODE}/unvalidated_solution.jls";fi
 done
done
python3 spar_ddp/experiments/value_group_solves/compare.py "$C" all 24 > "$C/comparison.json"
date -u +%FT%TZ > "$C/campaign.finished"
