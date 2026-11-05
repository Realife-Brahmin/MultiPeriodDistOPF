#!/usr/bin/env bash
set -euo pipefail
C="$SLURM_SUBMIT_DIR"; cd "$C"
export REPLAY_REPO=/home/aryan.jha/kamiak_jobs/stage_b_replay_oct15/repo
export REPLAY_CAPTURE=/home/aryan.jha/kamiak_jobs/stage_b_replay_oct15/capture
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_DYNAMIC=FALSE OPENBLAS_MAIN_FREE=1
unset GOMP_CPU_AFFINITY KMP_AFFINITY OMP_PROC_BIND OMP_PLACES
export FILTERDDP_TREE_KKT=1 REPLAY_MULTINODE=0
export PATH="$HOME/bin:$PATH"
{ date -u +%FT%TZ; hostname; lscpu; scontrol show job "$SLURM_JOB_ID"; taskset -pc $$; command -v julia; julia --version; } > hardware.txt
grep -q 'Gold 6438Y+' hardware.txt
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER"; exit 2; }
sha256sum -c source.sha256 > source-check.log
sha256sum -c capture.sha256 > capture-check.log
sha256sum "$REPLAY_REPO"/ddp/examples/power_system/{ieee123c_filterddp,procs_pass,paired_pass}.jl > numerical-source.sha256
julia --startup-file=no --project="$REPLAY_REPO/envs/ddp2026" dispatch.jl > dispatch.log 2>&1
julia --startup-file=no --project="$REPLAY_REPO/envs/ddp2026" unit.jl > unit.log 2>&1
grep -q RIGHT_KERNEL_SYNTHETIC_PASS unit.log
export VALUE_FOLLOWUP_BASE=/home/aryan.jha/kamiak_jobs/value_coupling_replay_oct22_v3
CAP=/home/aryan.jha/kamiak_jobs/compact_campaign_capture_oct23_v3
export VALUE_REPLAY_NUMERICAL_ROOT="$CAP/source"
export JULIA_DEPOT_PATH="$C/julia-depot:$CAP/julia-depot:$VALUE_FOLLOWUP_BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_LEAN_EXACT=1 FILTERDDP_LEAN_VALUE=1 FILTERDDP_HALF_BLOCK=1 FILTERDDP_BLOCKED_SOLVE=16
sha256sum -c "$CAP/runtime.sha256" > late-runtime-check.log
export STUDY_BLAS_THREADS=1
timeout --signal=TERM --kill-after=20s 10m /home/aryan.jha/opt/julia-1.12.4/bin/julia --startup-file=no --project="$CAP/source/envs/spar_ddp" late.jl "$CAP/results/large10kC_1ph_T24_campaign_capture/captures/pass0023_k0014_j008" > late-1.log 2>&1
grep -q 'LATE_COMPLETE accurate=true' late-1.log
# The first capture uses the historical project; select its existing depot again.
unset JULIA_DEPOT_PATH
run_config() {
    local label=$1 workers=$2 threads=$3 kernels=$4 widths=$5
    local directory="$C/$label"
    mkdir "$directory"
    (cd "$directory"; python3 "$C/affinity.py" "$workers" $widths)
    ln -s "$C/right_kernel.jl" "$directory/right_kernel.jl"
    export REPLAY_ROOT="$directory" STUDY_WORKERS="$workers" STUDY_BLAS_THREADS="$threads" STUDY_KERNELS="$kernels"
    timeout --signal=TERM --kill-after=20s 12m julia --startup-file=no --project="$REPLAY_REPO/envs/ddp2026" "$C/replay.jl" > "$directory/replay.log" 2>&1
    grep -q BALANCED_COMPLETE "$directory/replay.log"
}
run_config kernel63 63 1 current,right "1"
run_config threads31 31 1 current_1,current_2 "1 2"
run_config threads15 15 1 current_1,current_4 "1 4"
python3 choose.py "$C" > choice.txt
read -r workers threads < choice.txt
if [[ $workers != 0 ]]; then
    run_config combination "$workers" "$threads" current,right "$threads"
    export STUDY_BLAS_THREADS="$threads"
    export JULIA_DEPOT_PATH="$C/julia-depot:$CAP/julia-depot:$VALUE_FOLLOWUP_BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
    timeout --signal=TERM --kill-after=20s 10m /home/aryan.jha/opt/julia-1.12.4/bin/julia --startup-file=no --project="$CAP/source/envs/spar_ddp" late.jl "$CAP/results/large10kC_1ph_T24_campaign_capture/captures/pass0023_k0014_j008" > "late-$threads.log" 2>&1
    grep -q 'LATE_COMPLETE accurate=true' "late-$threads.log"
fi
date -u +%FT%TZ > study.finished
