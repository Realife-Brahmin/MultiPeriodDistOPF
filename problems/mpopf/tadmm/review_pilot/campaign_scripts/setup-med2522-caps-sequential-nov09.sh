#!/usr/bin/env bash
set -euo pipefail
C=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_med2522_caps_sequence
OLD=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_persistent_caps
test ! -e "$C"
mkdir -p "$C/campaign" "$C/logs"
cp -a "$OLD/source" "$C/source"
cp -a "$OLD/julia-depot" "$C/julia-depot"
cp "$OLD/input.sha256" "$OLD/FULL_PERIOD_READY" "$C/"
cp "$OLD/campaign/audit_sweep_energy.jl" "$OLD/campaign/audit.sh" "$C/campaign/"
cp /home/aryan.jha/kamiak_jobs/med2522-caps-sequence-overlay.tar.gz "$C/campaign/overlay.tar.gz"
tar -xzf "$C/campaign/overlay.tar.gz" -C "$C/source"
printf '%s\n' "$1" > "$C/source/source.commit"
cd "$C/source";P=problems/mpopf/tadmm/review_pilot
sha256sum --check --quiet "$C/input.sha256"
find "$P" test_networks problems/mpopf/spar_ddp/examples/power_system -type f ! -path '*/results/*' ! -path '*/__pycache__/*' -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256"
cp "$OLD/campaign/gate.sh" "$C/campaign/cap-gate.sh"
cat > "$C/campaign/gate.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1;P=problems/mpopf/tadmm/review_pilot
bash "$C/campaign/cap-gate.sh" "$C"
source /home/aryan.jha/kamiak_jobs/tadmm_review_nov09_sequential_dss_v2/campaign/env.sh
cd "$C/source"
sha256sum --check --quiet "$C/execution-source.sha256"
julia --startup-file=no --project="$P" "$P/validate_sequential_dss.jl"
# Existing two-node endpoints are used only to diagnose chronological physics.
# Their timing is not a replacement for the requested single-node benchmark.
CASE=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_historical_anchors/results/run_T96_primal_dual_31660059
export SEQUENTIAL_DSS_GZIP=1
for NAME in central endpoint;do
 if [[ $NAME == central ]];then SOL="$CASE/central_dispatch.jls";else SOL="$CASE/rho20000/endpoint.jls";fi
 julia --startup-file=no --project="$P" "$P/run_sequential_dss.jl" "$CASE/instance.jls" "$SOL" "$C/saved-T96-$NAME-model" model
done
date -u +%FT%TZ > "$C/MED_SEQUENCE_GATE_VALIDATED"
BASH
cat > "$C/campaign/run.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1;CAP=$2;test -s "$C/MED_SEQUENCE_GATE_VALIDATED"
[[ $SLURM_JOB_NUM_NODES == 1 ]]
export TADMM_CPU_MODEL_MATCH=6760P
export TADMM_RHO_POLICY=bootstrap_fixed TADMM_RHO_FLOOR=48
export TADMM_LOCAL_SOLVE_POLICY=persistent_cap TADMM_CAP_ITERATIONS=$CAP TADMM_CAP_SCREEN=finite
# make_model(local_window=false) forces exact policy for centralized.
# Centralized retains max_iter=3000 and must return solved, never iteration limit.
exec bash "$C/source/problems/mpopf/tadmm/review_pilot/kamiak_full_period.sh" "$C" 96 run primal_dual 20000 ma57 ieee2522C_1ph
BASH
for F in "$C/campaign/"*.sh;do bash -n "$F";done
squeue -u aryan.jha > "$C/campaign/initial-queue.txt"
G=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:15:00 --job-name=med2522-sequence-gate --output="$C/logs/gate-%j.out" --error="$C/logs/gate-%j.err" "$C/campaign/gate.sh" "$C")
printf 'gate=%s\n' "$G" > "$C/campaign/submissions.txt"
DEP=$G
for CAP in 10 30;do
 J=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=128 --hint=nomultithread --mem=128G --time=03:00:00 --constraint='graniterapids&6760P-2.20ghz&cores-128' --dependency="afterok:$DEP" --kill-on-invalid-dep=yes --job-name="tadmm-med96-always$CAP" --output="$C/logs/cap$CAP-%j.out" --error="$C/logs/cap$CAP-%j.err" "$C/campaign/run.sh" "$C" "$CAP")
 A=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:15:00 --dependency="afterany:$J" --job-name="med96-cap$CAP-energy" --output="$C/logs/audit-%j.out" --error="$C/logs/audit-%j.err" "$C/campaign/audit.sh" "$C" 96 "$J")
 printf 'cap%s=%s\ncap%s_audit=%s\n' "$CAP" "$J" "$CAP" "$A" >> "$C/campaign/submissions.txt"
 DEP=$J
done
cat "$C/campaign/submissions.txt"
