#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1");MODE=$2;V=${3:-};T=${4:-}
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=2
META="$C/jobs/$SLURM_JOB_ID";mkdir -p "$META"
exec > >(tee "$META/stdout.log") 2> >(tee "$META/stderr.log" >&2)
trap 'c=/sys/fs/cgroup/memory/slurm/uid_$(id -u)/job_${SLURM_JOB_ID}; for f in memory.max_usage_in_bytes memory.limit_in_bytes; do [[ ! -f "$c/$f" ]] || echo "$f $(cat "$c/$f")"; done > "$META/node_memory.txt"' EXIT
{ hostname; julia --version; lscpu; scontrol show job "$SLURM_JOB_ID"; taskset -pc $$; } > "$META/hardware.txt"
grep -q '6438Y+' "$META/hardware.txt"; [[ $(julia --version) == 'julia version 1.12.4' ]]
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "Other node jobs: $OTHER";exit 2; }
if [[ $MODE == prepare ]]; then
 cd "$C/source"
 julia --startup-file=no --project=envs/spar_ddp -e 'using Pkg;Pkg.develop(path="spar_ddp/DDP4OPF.jl");Pkg.instantiate();Pkg.precompile()' > "$META/setup-ddp.log" 2>&1
 julia --startup-file=no --project=envs/tadmm -e 'using Pkg;Pkg.instantiate();using OpenDSSDirect' > "$META/setup-opendss.log" 2>&1
 mkdir -p spar_ddp/results/network_filterddp
 PROFILE_PERIODIC=1 julia --startup-file=no --project=envs/tadmm spar_ddp/examples/power_system/export_ieee123c_data.jl large10kC_1ph 3 > "$META/export-T3.log" 2>&1
 cp spar_ddp/results/network_filterddp/network_data_large10kC_1ph_T3_periodic.jls "$C/inputs/"
 julia --startup-file=no spar_ddp/experiments/battery_count/prepare.jl "$C" > "$META/construction.log" 2>&1
 for V in 10k_X 10k_Y; do
  cd "$C/variants/$V/source"
  julia --startup-file=no --project=envs/spar_ddp -e 'using Pkg;Pkg.develop(path="spar_ddp/DDP4OPF.jl");Pkg.instantiate();using DDP4OPF,JuMP,Ipopt' > "$META/$V-environment.log" 2>&1
  cp envs/spar_ddp/Manifest.toml "$C/variants/$V/Manifest-ddp.toml"
  cp envs/tadmm/Manifest.toml "$C/variants/$V/Manifest-tadmm.toml"
  find spar_ddp envs rawData -type f -not -path '*/results/*' -print0 | sort -z | xargs -0 sha256sum > "$C/variants/$V/execution-source.sha256"
  sha256sum spar_ddp/results/network_filterddp/network_data_*.jls > "$C/variants/$V/input.sha256"
  for M in ma97 ma57 ddp; do bash spar_ddp/experiments/battery_count/run.sh "$C" "$V" 3 "$M" 0; done
  python3 spar_ddp/experiments/battery_count/check_group.py "$C" "$V" 3 > "$META/$V-smoke-validation.json"
 done
 date -u +%FT%TZ > "$C/smoke.complete"
elif [[ $MODE == case ]]; then
 test -s "$C/smoke.complete"
 cd "$C/variants/$V/source"
 sha256sum -c "$C/variants/$V/execution-source.sha256" > "$META/source-check.txt"
 sha256sum -c "$C/variants/$V/input.sha256" > "$META/input-check.txt"
 for R in 1 2; do
  ORDER='ma97 ma57 ddp'; [[ $R == 2 ]] && ORDER='ddp ma57 ma97'
  for M in $ORDER; do bash spar_ddp/experiments/battery_count/run.sh "$C" "$V" "$T" "$M" "$R"; done
 done
 python3 spar_ddp/experiments/battery_count/check_group.py "$C" "$V" "$T" > "$META/group-validation.json"
 date -u +%FT%TZ > "$META/COMPLETE"
else exit 2; fi
