#!/usr/bin/env bash
set -euo pipefail
ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
C=$(realpath -m "$1")
[[ ! -e "$C" && -z $(git -C "$ROOT" status --porcelain --untracked-files=no) ]] || exit 2
INPUT=/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/source/spar_ddp/results/network_filterddp
mkdir -p "$C/source" "$C/results" "$C/logs"
git -C "$ROOT" rev-parse HEAD > "$C/source.commit"
git -C "$ROOT" archive HEAD AGENTS.md spar_ddp/DDP4OPF.jl spar_ddp/examples/power_system spar_ddp/experiments/screening_ablation spar_ddp/experiments/tree_ablation envs/spar_ddp envs/tadmm rawData/LoadShapeDefault.dss rawData/LoadShapePVDefault.dss rawData/LoadShapePSubsCostDefault.dss rawData/ieee123C_1ph rawData/ieee2522C_1ph rawData/large10kC_1ph | tar --no-same-owner -x -C "$C/source"
touch "$C/source/.screening-ablation-private"
mkdir -p "$C/source/spar_ddp/results/network_filterddp"
for F in network_data_ieee123C_1ph_T3_periodic.jls network_data_ieee123C_1ph_T96_periodic.jls network_data_ieee2522C_1ph_T96_periodic.jls network_data_large10kC_1ph_T96_periodic.jls; do
 cp "$INPUT/$F" "$C/source/spar_ddp/results/network_filterddp/"
done
(cd "$C/source" && sha256sum spar_ddp/results/network_filterddp/network_data_*.jls) > "$C/input.sha256"
sbatch --parsable --job-name=ddp-screen-factorial --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=600G --time=08:00:00 --partition=kamiak --account=dubey --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --output="$C/logs/slurm-%j.out" --error="$C/logs/slurm-%j.err" "$C/source/spar_ddp/experiments/screening_ablation/job.sh" "$C" | tee "$C/job.id"
