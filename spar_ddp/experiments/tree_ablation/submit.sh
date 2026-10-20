#!/usr/bin/env bash
# Submit one serialized, exclusive-node diagnostic; never cancel other jobs.
set -euo pipefail
[[ $# == 2 ]] || { echo 'usage: submit.sh INPUT_EXPORT_DIRECTORY NEW_CAMPAIGN_DIRECTORY' >&2; exit 2; }
ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
INPUT=$(realpath "$1"); CAMPAIGN=$(realpath -m "$2")
command -v sbatch >/dev/null || { echo 'Run this command from a Kamiak-connected shell with Slurm available' >&2; exit 2; }
[[ ! -e "$CAMPAIGN" ]] || { echo 'Campaign directory already exists; refusing overwrite' >&2; exit 2; }
[[ -z $(git -C "$ROOT" status --porcelain --untracked-files=no) ]] || { echo 'Commit source changes before submitting' >&2; exit 2; }
# Validate data before reserving resources; exports are inputs, not regenerated.
for SYS in ieee123C_1ph ieee2522C_1ph large10kC_1ph; do
 [[ -s "$INPUT/network_data_${SYS}_T96_periodic.jls" ]] || { echo "Missing input: $SYS T96 periodic" >&2; exit 2; }
done
[[ -s "$INPUT/network_data_ieee123C_1ph_T3_periodic.jls" || ${TREE_ABLATION_EXPORT_SMOKE:-0} == 1 ]] || { echo 'Missing IEEE123 T3 periodic smoke input' >&2; exit 2; }
mkdir -p "$CAMPAIGN/source" "$CAMPAIGN/results" "$CAMPAIGN/logs"
git -C "$ROOT" rev-parse HEAD > "$CAMPAIGN/source.commit"
git -C "$ROOT" archive HEAD AGENTS.md spar_ddp/DDP4OPF.jl spar_ddp/examples/power_system \
 spar_ddp/experiments/tree_ablation envs/spar_ddp envs/tadmm \
 rawData/ieee123C_1ph rawData/ieee2522C_1ph rawData/large10kC_1ph | tar --no-same-owner -x -C "$CAMPAIGN/source"
touch "$CAMPAIGN/source/.tree-ablation-private"
mkdir -p "$CAMPAIGN/source/spar_ddp/results/network_filterddp"
for SYS in ieee123C_1ph ieee2522C_1ph large10kC_1ph; do
 cp "$INPUT/network_data_${SYS}_T96_periodic.jls" "$CAMPAIGN/source/spar_ddp/results/network_filterddp/"
done
if [[ -s "$INPUT/network_data_ieee123C_1ph_T3_periodic.jls" ]]; then
 cp "$INPUT/network_data_ieee123C_1ph_T3_periodic.jls" "$CAMPAIGN/source/spar_ddp/results/network_filterddp/"
else
 touch "$CAMPAIGN/export-periodic-smoke"
fi
(cd "$CAMPAIGN/source" && find spar_ddp envs rawData -type f -print0 | sort -z | xargs -0 sha256sum) > "$CAMPAIGN/source.sha256"
# A single exclusive allocation executes the two arms and all repeats in turn.
# Inherits the site's normal account/partition; caller may use SBATCH_* settings.
JOB=$(sbatch --parsable --job-name=ddp-tree-ablation --nodes=1 --ntasks=1 \
 --cpus-per-task=64 --mem=600G --time=04:00:00 \
 --partition=kamiak --account=dubey --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' \
 --output="$CAMPAIGN/logs/slurm-%j.out" --error="$CAMPAIGN/logs/slurm-%j.err" \
 "$CAMPAIGN/source/spar_ddp/experiments/tree_ablation/slurm_job.sh" "$CAMPAIGN")
printf '%s\n' "$JOB" | tee "$CAMPAIGN/job.id"
