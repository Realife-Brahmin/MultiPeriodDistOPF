#!/usr/bin/env bash
set -euo pipefail
C=$SLURM_SUBMIT_DIR;cd "$C"
CAP=/home/aryan.jha/kamiak_jobs/phase_b_influence_capture_nov05
OLD=/home/aryan.jha/kamiak_jobs/value_coupling_audit_oct22
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$CAP/julia-depot:$OLD/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VALUE_REPLAY_NUMERICAL_ROOT="$CAP/source"
export FILTERDDP_TREE_KKT=1 FILTERDDP_FACTOR_BACKED_POLICY=1
exec > >(tee "$C/stdout.log") 2> >(tee "$C/stderr.log" >&2)
{ date -u +%FT%TZ;hostname;julia --version;lscpu;taskset -pc $$;scontrol show job "$SLURM_JOB_ID"; } > hardware.txt
grep -q 'Gold 6438Y+' hardware.txt
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER";exit 2; }
test -f "$CAP/capture.finished"
sha256sum -c source.sha256 > source-check.log
(cd "$CAP/source";sha256sum -c "$CAP/runtime.sha256";sha256sum -c "$CAP/input.sha256") > capture-runtime-input-check.log
julia --startup-file=no influence_unit.jl > unit.log 2>&1
julia --startup-file=no --project="$CAP/source/envs/spar_ddp" runtime_unit.jl > runtime-unit.log 2>&1
INPUT="$CAP/source/spar_ddp/results/network_filterddp/network_data_large10kC_1ph_T24_periodic.jls"
# Batch2h is the overall bound, including diagnostic baseline/SVD/reconstruction.
timeout --signal=TERM --kill-after=30s 110m julia --startup-file=no --project="$CAP/source/envs/spar_ddp" influence_replay.jl "$CAP/results/large10k_T24_influence_capture/captures" "$INPUT" "$C/results" > replay.log 2>&1
sha256sum -c source.sha256 > source-check-after.log
(cd "$CAP/source";sha256sum -c "$CAP/runtime.sha256";sha256sum -c "$CAP/input.sha256") > capture-runtime-input-check-after.log
date -u +%FT%TZ > replay.finished
