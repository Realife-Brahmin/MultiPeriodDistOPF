#!/usr/bin/env bash
set -euo pipefail
OLD=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_targeted
C=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_mumps_fixed
test ! -e "$C"
mkdir -p "$C/campaign" "$C/logs"
cp -a "$OLD/source" "$C/source"
cp -a "$OLD/julia-depot" "$C/julia-depot"
cp "$OLD/input.sha256" "$OLD/execution-source.sha256" "$OLD/FULL_PERIOD_READY" "$C/"
cp "$OLD/campaign/audit_sweep_energy.jl" "$C/campaign/"
cd "$C/source"
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
cat > "$C/campaign/gate.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1; P=problems/mpopf/tadmm/review_pilot
cd "$C/source"
export JULIA_DEPOT_PATH="$C/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov08_horizon_sweep/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_accuracy_fix/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export TADMM_LINEAR_SOLVER=mumps TADMM_IPOPT_OPTIONS=linear_system_scaling=none
/home/aryan.jha/opt/julia-1.12.4/bin/julia --startup-file=no --project="$P" "$P/validate_stage_reuse.jl"
date -u +%FT%TZ > "$C/MUMPS_FIXED_GATE_VALIDATED"
BASH
cat > "$C/campaign/run.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1; T=$2; RHO=$3
test -s "$C/MUMPS_FIXED_GATE_VALIDATED"
if [[ $T == 96 ]]; then
 test -s "$C/T48_VALIDATED"
fi
export TADMM_RHO_POLICY=bootstrap_fixed TADMM_RHO_FLOOR=48 TADMM_LOCAL_SOLVE_POLICY=exact
bash "$C/source/problems/mpopf/tadmm/review_pilot/kamiak_full_period.sh" "$C" "$T" run primal_dual "$RHO" mumps
if [[ $T == 48 ]]; then
 test -s "$C/results/run_T48_primal_dual_$SLURM_JOB_ID/SWEEP_VALIDATED"
 date -u +%FT%TZ > "$C/T48_VALIDATED"
fi
BASH
cat > "$C/campaign/audit.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1; T=$2; J=$3; CASE="$C/results/run_T${T}_primal_dual_$J"
if [[ ! -s "$CASE/SWEEP_VALIDATED" ]]; then echo 'No validated endpoint; no completed audit claimed.'; exit 0; fi
export TADMM_AUDIT_LAYOUT="$C/source/problems/mpopf/spar_ddp/examples/power_system/control_layout.jl"
export TADMM_AUDIT_VALIDATED_ONLY=1 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/home/aryan.jha/opt/julia-1.12.4/bin/julia --startup-file=no "$C/campaign/audit_sweep_energy.jl" "$CASE" "$C/campaign/T$T-$J-energy.toml"
BASH
for F in "$C/campaign/"*.sh; do bash -n "$F"; done
squeue -u aryan.jha > "$C/campaign/initial-queue.txt"
G=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:10:00 --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --job-name=tadmm-mumps-fixed-gate --output="$C/logs/gate-%j.out" --error="$C/logs/gate-%j.err" "$C/campaign/gate.sh" "$C")
printf 'gate=%s\n' "$G" > "$C/campaign/submissions.txt"
DEP=$G
for T in 48 96; do
 N=$((T/48)); RHO=$((T*2500))
 J=$(sbatch --parsable --account=dubey --partition=kamiak --nodes="$N" --ntasks="$N" --ntasks-per-node=1 --cpus-per-task=64 --hint=nomultithread --mem=128G --time=03:00:00 --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --dependency="afterok:$DEP" --kill-on-invalid-dep=yes --job-name="tadmm-mumps-fixed-T$T" --output="$C/logs/T$T-%j.out" --error="$C/logs/T$T-%j.err" "$C/campaign/run.sh" "$C" "$T" "$RHO")
 A=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:15:00 --dependency="afterany:$J" --job-name="tadmm-mumps-energy-T$T" --output="$C/logs/audit-%j.out" --error="$C/logs/audit-%j.err" "$C/campaign/audit.sh" "$C" "$T" "$J")
 printf 'T%s=%s\nT%s_audit=%s\n' "$T" "$J" "$T" "$A" >> "$C/campaign/submissions.txt"
 DEP=$J
done
cat "$C/campaign/submissions.txt"
