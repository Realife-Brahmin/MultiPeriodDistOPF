#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1");cd "$C/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=2
unset NO_COLOR FORCE_COLOR
{ hostname;julia --version;lscpu;taskset -pc $$;scontrol show job "$SLURM_JOB_ID"; } > "$C/hardware.txt"
grep -q 'Gold 6438Y+' "$C/hardware.txt";[[ $(julia --version) == 'julia version 1.12.4' ]]
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER";exit 2; }
sha256sum -c "$C/preparation-source.sha256" > "$C/logs/source-check.log"
julia --startup-file=no --project=envs/spar_ddp -e 'using Pkg;Pkg.develop(path="spar_ddp/DDP4OPF.jl");Pkg.instantiate();Pkg.precompile()' > "$C/logs/setup-ddp.log" 2>&1
julia --startup-file=no --project=envs/tadmm -e 'using Pkg;Pkg.instantiate();using OpenDSSDirect' > "$C/logs/setup-opendss.log" 2>&1
julia --startup-file=no spar_ddp/experiments/value_coupling_audit/unit.jl > "$C/logs/unit.log" 2>&1
bash spar_ddp/experiments/value_coupling_audit/run_capture.sh "$C" ieee123C_1ph 3 2
date -u +%FT%TZ > "$C/smoke.complete"
bash spar_ddp/experiments/value_coupling_audit/run_capture.sh "$C" ieee123C_1ph 24 11
bash spar_ddp/experiments/value_coupling_audit/run_capture.sh "$C" ieee2522C_1ph 24 26
bash spar_ddp/experiments/value_coupling_audit/run_capture.sh "$C" large10kC_1ph 24 7
date -u +%FT%TZ > "$C/capture.finished"
