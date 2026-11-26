#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1"); cd "$C/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=2
hostname > "$C/hostname"; lscpu > "$C/lscpu.txt"
scontrol show job "$SLURM_JOB_ID" > "$C/allocation.txt"
taskset -pc $$ > "$C/affinity.txt"
[[ $(julia --version) == 'julia version 1.12.4' ]] || exit 2
grep -q '6438Y+' "$C/lscpu.txt" || exit 2
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "Occupied node: $OTHER"; exit 2; }
julia --startup-file=no --project=envs/spar_ddp -e 'using Pkg; Pkg.develop(path="spar_ddp/DDP4OPF.jl"); Pkg.instantiate(); Pkg.precompile()' > "$C/logs/setup-ddp.log" 2>&1
julia --startup-file=no --project=envs/tadmm -e 'using Pkg; Pkg.instantiate(); using OpenDSSDirect' > "$C/logs/setup-opendss.log" 2>&1
cp envs/spar_ddp/Manifest.toml "$C/Manifest-ddp.toml"
cp envs/tadmm/Manifest.toml "$C/Manifest-tadmm.toml"
find spar_ddp/DDP4OPF.jl spar_ddp/examples spar_ddp/experiments/screening_ablation envs rawData -type f -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256"
# All-eight-arm smoke gates the full factorial. Mask trajectories may differ.
for M in 7 0 1 2 3 4 5 6; do
 SCREEN_WORKERS=2 bash spar_ddp/experiments/screening_ablation/run_arm.sh ieee123C_1ph 3 0 "$M" "$C/smoke/m$M"
done
for SYS in ieee123C_1ph ieee2522C_1ph large10kC_1ph; do
 for REP in 1 2; do
  ORDER='7 0 1 2 3 4 5 6'; [[ $REP == 2 ]] && ORDER='6 5 4 3 2 1 0 7'
  for M in $ORDER; do
   OUT="$C/results/${SYS}_T96_r${REP}_m${M}"
   date -u +%FT%TZ | tee -a "$C/progress.log"
   if bash spar_ddp/experiments/screening_ablation/run_arm.sh "$SYS" 96 "$REP" "$M" "$OUT"; then
    echo "PASS $SYS $REP $M" | tee -a "$C/progress.log"
   else
    CODE=$?; echo "FAIL $SYS $REP $M exit=$CODE" | tee -a "$C/progress.log"
    if pgrep -u "$USER" -f "$C/source" > /dev/null; then echo 'Lingering workers: stop safely'; exit 1; fi
   fi
  done
 done
done
date -u +%FT%TZ > "$C/campaign.finished"
