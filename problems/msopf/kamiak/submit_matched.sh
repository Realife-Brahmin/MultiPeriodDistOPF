#!/usr/bin/env bash
set -euo pipefail
cd "$HOME/msopf/jobs/matched_secondary_oct25"
[[ ! -e submitted.tsv ]] || { echo 'Campaign already submitted';exit 2; }
squeue -u "$USER" -o '%.18i %.20j %.8T %.6C %.12P %.30R' > queue_before.txt
boot=$(sbatch --parsable --partition=vcea --account=dubey --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=6G --time=00:15:00 --job-name=msm_boot --output=logs/%x_%j.out --error=logs/%x_%j.err matched_job.sh bootstrap init)
boot=${boot%%;*};printf 'bootstrap\t%s\n' "$boot" >> submitted.tsv
for label in primary4 secondary4 sf4 primary3 secondary3 primary2 secondary2; do
 prep=$(sbatch --parsable --dependency="afterok:$boot" --partition=vcea --account=dubey --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=12G --time=00:20:00 --job-name="msm_${label}_prep" --output=logs/%x_%j.out --error=logs/%x_%j.err matched_job.sh "$label" prep)
 prep=${prep%%;*}
 count=$(wc -l < "tasks_${label}.txt")
 concurrency=6;[[ "$label" == *3 ]] && concurrency=2;[[ "$label" == *2 ]] && concurrency=1
 arr=$(sbatch --parsable --dependency="afterok:$prep" --array="0-$((count-1))%$concurrency" --partition=vcea --account=dubey --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=12G --time=00:45:00 --job-name="msm_${label}" --output=logs/%x_%A_%a.out --error=logs/%x_%A_%a.err matched_job.sh "$label" array)
 arr=${arr%%;*};printf '%s\t%s\t%s\t%s\n' "$label" "$prep" "$arr" "$count" >> submitted.tsv
done
cat submitted.tsv
