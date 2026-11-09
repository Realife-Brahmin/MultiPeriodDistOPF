#!/usr/bin/env bash
set -euo pipefail
OLD=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_historical_anchors
C=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_med2522_single_node
test ! -e "$C"
mkdir -p "$C/campaign" "$C/logs"
cp -a "$OLD/source" "$C/source"
cp -a "$OLD/julia-depot" "$C/julia-depot"
cp "$OLD/input.sha256" "$OLD/FULL_PERIOD_READY" "$OLD/ANCHOR_GATE_VALIDATED" "$C/"
cp "$OLD/campaign/audit_sweep_energy.jl" "$OLD/campaign/audit.sh" "$C/campaign/"
cp /home/aryan.jha/kamiak_jobs/med2522-single-node-overlay.tar.gz "$C/campaign/overlay.tar.gz"
tar -xzf "$C/campaign/overlay.tar.gz" -C "$C/source"
printf '%s\n' "$1" > "$C/source/source.commit"
cd "$C/source";P=problems/mpopf/tadmm/review_pilot
sha256sum --check --quiet "$C/input.sha256"
find "$P" test_networks problems/mpopf/spar_ddp/examples/power_system -type f ! -path '*/results/*' -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256"
bash -n "$P/kamiak_full_period.sh"
cat > "$C/campaign/run.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1; SOLVER=$2
test -s "$C/ANCHOR_GATE_VALIDATED"
[[ $SLURM_JOB_NUM_NODES == 1 ]]
export TADMM_CPU_MODEL_MATCH=6760P
export TADMM_RHO_POLICY=bootstrap_fixed TADMM_RHO_FLOOR=48 TADMM_LOCAL_SOLVE_POLICY=exact
exec bash "$C/source/problems/mpopf/tadmm/review_pilot/kamiak_full_period.sh" "$C" 96 run primal_dual 20000 "$SOLVER" ieee2522C_1ph
BASH
bash -n "$C/campaign/run.sh"
squeue -u aryan.jha > "$C/campaign/initial-queue.txt"
DEP=''
for SOLVER in ma57 mumps; do
 OPTIONS=()
 if [[ -n $DEP ]]; then OPTIONS=(--dependency="afterok:$DEP" --kill-on-invalid-dep=yes); fi
 J=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=128 --hint=nomultithread --mem=128G --time=03:00:00 --constraint='graniterapids&6760P-2.20ghz&cores-128' "${OPTIONS[@]}" --job-name="tadmm-med2522-single-$SOLVER" --output="$C/logs/$SOLVER-%j.out" --error="$C/logs/$SOLVER-%j.err" "$C/campaign/run.sh" "$C" "$SOLVER")
 A=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:15:00 --dependency="afterany:$J" --job-name=tadmm-single-energy --output="$C/logs/audit-%j.out" --error="$C/logs/audit-%j.err" "$C/campaign/audit.sh" "$C" 96 "$J")
 printf '%s=%s\n%s_audit=%s\n' "$SOLVER" "$J" "$SOLVER" "$A" >> "$C/campaign/submissions.txt"
 DEP=$J
done
cat "$C/campaign/submissions.txt"
