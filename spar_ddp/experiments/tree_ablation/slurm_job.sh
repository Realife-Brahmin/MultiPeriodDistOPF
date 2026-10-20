#!/usr/bin/env bash
set -euo pipefail
CAMPAIGN=$(realpath "$1");cd "$CAMPAIGN/source"
export TREE_ABLATION_SOURCE_SHA=$(cat "$CAMPAIGN/source.commit")
export JULIA_DEPOT_PATH="$CAMPAIGN/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=2
export REDUCED_PROFILE=periodic
unset NO_COLOR FORCE_COLOR
# Load Julia before submission, or supply an explicit module via the site's shell.
command -v julia >/dev/null || { echo 'Julia is not on PATH; load Julia 1.12.4 before submission' >&2; exit 2; }
julia --version | tee "$CAMPAIGN/julia.version"
hostname > "$CAMPAIGN/hostname";lscpu > "$CAMPAIGN/lscpu.txt"
scontrol show job "$SLURM_JOB_ID" > "$CAMPAIGN/allocation.txt"
# Separate projects/depot: no Pkg operations against MSOPF or shared environments.
julia --startup-file=no --project=envs/spar_ddp -e 'using Pkg; Pkg.develop(path="spar_ddp/DDP4OPF.jl"); Pkg.instantiate(); Pkg.precompile()' > "$CAMPAIGN/logs/setup-ddp.log" 2>&1
julia --startup-file=no --project=envs/tadmm -e 'using Pkg; Pkg.instantiate(); using OpenDSSDirect' > "$CAMPAIGN/logs/setup-opendss.log" 2>&1
cp envs/spar_ddp/Manifest.toml "$CAMPAIGN/Manifest-spar_ddp.toml"
cp envs/tadmm/Manifest.toml "$CAMPAIGN/Manifest-tadmm.toml"
julia --startup-file=no --project=envs/spar_ddp spar_ddp/examples/power_system/tree_ablation_smoke.jl > "$CAMPAIGN/logs/synthetic-smoke.log" 2>&1
export FILTERDDP_PROCS_WORKERS=2
bash spar_ddp/experiments/tree_ablation/run_pair.sh ieee123C_1ph 3 0 "$CAMPAIGN/results/smoke" tree
export FILTERDDP_PROCS_WORKERS=63
for SYS in ieee123C_1ph ieee2522C_1ph large10kC_1ph; do
 for REP in 1 2 3; do
  FIRST=tree; (( REP % 2 == 0 )) && FIRST=umfpack
  bash spar_ddp/experiments/tree_ablation/run_pair.sh "$SYS" 96 "$REP" "$CAMPAIGN/results/${SYS}_T96_r${REP}" "$FIRST"
 done
done
python3 spar_ddp/experiments/tree_ablation/summarize.py "$CAMPAIGN/results" > "$CAMPAIGN/summary.csv"
printf 'TREE_ABLATION_CAMPAIGN_COMPLETE %s\n' "$(date -u +%FT%TZ)" > "$CAMPAIGN/campaign.complete"
