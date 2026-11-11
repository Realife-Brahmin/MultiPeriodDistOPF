#!/usr/bin/env bash
set -euo pipefail
JOBS=/home/aryan.jha/kamiak_jobs
BASE="$JOBS/tadmm_review_nov09_med2522_caps_sequence"
POLICY=${1:-cap10}
case "$POLICY" in
 cap10|exact) ;;
 *) echo "Unknown local solve policy: $POLICY" >&2; exit 2 ;;
esac
C="$JOBS/tadmm_review_nov10_large10k_t48_$POLICY"
LAUNCH_FILES="$JOBS/large10k-t48-cap10-nov10-launch"
test -s "$BASE/MED_SEQUENCE_GATE_VALIDATED"
test ! -e "$C"
mkdir -p "$C/campaign" "$C/logs"
cp -a "$BASE/source" "$BASE/julia-depot" "$C/"
cp "$BASE/input.sha256" "$BASE/execution-source.sha256" "$BASE/FULL_PERIOD_READY" "$BASE/MED_SEQUENCE_GATE_VALIDATED" "$C/"
cp "$BASE/campaign/audit_sweep_energy.jl" "$BASE/campaign/audit.sh" "$C/campaign/"
cp "$LAUNCH_FILES/audit-socp-storage-nov09.py" "$LAUNCH_FILES/summarize-tadmm-round-timing-nov10.py" "$C/campaign/"
cp "$LAUNCH_FILES/posthoc-med2522-caps-nov09.sh" "$C/campaign/posthoc.sh"
cd "$C/source"
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
cat > "$C/campaign/run.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1
POLICY=$2
[[ ${SLURM_JOB_NUM_NODES:?} == 1 ]]
test -s "$C/MED_SEQUENCE_GATE_VALIDATED"
export TADMM_RHO_POLICY=bootstrap_fixed TADMM_RHO_FLOOR=48
if [[ $POLICY == cap10 ]]; then
 export TADMM_LOCAL_SOLVE_POLICY=persistent_cap TADMM_CAP_ITERATIONS=10 TADMM_CAP_SCREEN=finite
else
 export TADMM_LOCAL_SOLVE_POLICY=exact
fi
exec bash "$C/source/problems/mpopf/tadmm/review_pilot/kamiak_full_period.sh" "$C" 48 run primal_dual 120000 ma57 large10kC_1ph
BASH
for script in "$C/campaign/"*.sh; do bash -n "$script"; done
J=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=64 --hint=nomultithread --mem=128G --time=01:00:00 --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --job-name="tadmm-large48-$POLICY" --output="$C/logs/solve-%j.out" --error="$C/logs/solve-%j.err" "$C/campaign/run.sh" "$C" "$POLICY")
printf 'solve=%s\n' "$J" > "$C/campaign/submissions.txt"
A=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:15:00 --dependency="afterany:$J" --job-name=large48-energy --output="$C/logs/energy-%j.out" --error="$C/logs/energy-%j.err" "$C/campaign/audit.sh" "$C" 48 "$J")
printf 'energy_audit=%s\n' "$A" >> "$C/campaign/submissions.txt"
E=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:15:00 --dependency="afterany:$J" --job-name=large48-explicit-validation --output="$C/logs/explicit-%j.out" --error="$C/logs/explicit-%j.err" "$C/campaign/posthoc.sh" "$J" "$C" 48 120000 large10kC_1ph)
printf 'explicit_audit=%s\n' "$E" >> "$C/campaign/submissions.txt"
cat "$C/campaign/submissions.txt"
