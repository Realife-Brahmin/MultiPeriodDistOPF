#!/usr/bin/env bash
set -euo pipefail
C=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_persistent_caps
test -s "$C/PERSISTENT_CAP_GATE_VALIDATED"
cd "$C/source"
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
squeue -u aryan.jha
for CAP in 30 10;do
 J=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=64 --hint=nomultithread --mem=128G --time=01:00:00 --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --job-name="tadmm-large24-always$CAP" --output="$C/logs/cap$CAP-%j.out" --error="$C/logs/cap$CAP-%j.err" "$C/campaign/run.sh" "$C" "$CAP")
 A=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:15:00 --dependency="afterany:$J" --job-name="tadmm-always$CAP-energy" --output="$C/logs/audit-%j.out" --error="$C/logs/audit-%j.err" "$C/campaign/audit.sh" "$C" 24 "$J")
 printf 'resubmitted_cap%s=%s\nresubmitted_cap%s_audit=%s\n' "$CAP" "$J" "$CAP" "$A" >> "$C/campaign/submissions.txt"
done
cat "$C/campaign/submissions.txt"
