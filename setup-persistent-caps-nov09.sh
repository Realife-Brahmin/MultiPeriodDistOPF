#!/usr/bin/env bash
set -euo pipefail
OLD=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_early_cap_finite
C=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_persistent_caps
test ! -e "$C"
mkdir -p "$C/campaign" "$C/logs"
cp -a "$OLD/source" "$C/source"
cp -a "$OLD/julia-depot" "$C/julia-depot"
cp "$OLD/input.sha256" "$OLD/FULL_PERIOD_READY" "$C/"
cp "$OLD/campaign/audit_sweep_energy.jl" "$OLD/campaign/audit.sh" "$C/campaign/"
cp /home/aryan.jha/kamiak_jobs/persistent-caps-nov09-overlay.tar.gz "$C/campaign/overlay.tar.gz"
tar -xzf "$C/campaign/overlay.tar.gz" -C "$C/source"
printf '%s\n' "$1" > "$C/source/source.commit"
cd "$C/source";P=problems/mpopf/tadmm/review_pilot
sha256sum --check --quiet "$C/input.sha256"
find "$P" test_networks problems/mpopf/spar_ddp/examples/power_system -type f ! -path '*/results/*' -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256"
bash -n "$P/kamiak_full_period.sh"
python3 -m py_compile "$P/validate_horizon_sweep.py"
cat > "$C/campaign/gate.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1;P=problems/mpopf/tadmm/review_pilot;cd "$C/source"
export JULIA_DEPOT_PATH="$C/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov08_horizon_sweep/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_accuracy_fix/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export TADMM_LINEAR_SOLVER=ma57
export TADMM_IPOPT_OPTIONS="hsllib=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so;linear_system_scaling=none;ma57_pivot_order=2"
/home/aryan.jha/opt/julia-1.12.4/bin/julia --startup-file=no --project="$P" "$P/validate_early_cap.jl"
date -u +%FT%TZ > "$C/PERSISTENT_CAP_GATE_VALIDATED"
BASH
cat > "$C/campaign/run.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1;CAP=$2;test -s "$C/PERSISTENT_CAP_GATE_VALIDATED"
export TADMM_RHO_POLICY=bootstrap_fixed TADMM_RHO_FLOOR=48 TADMM_LOCAL_SOLVE_POLICY=persistent_cap
export TADMM_CAP_ITERATIONS=$CAP TADMM_CAP_SCREEN=finite
exec bash "$C/source/problems/mpopf/tadmm/review_pilot/kamiak_full_period.sh" "$C" 24 run primal_dual 15000 ma57 large10kC_1ph
BASH
for F in "$C/campaign/"*.sh;do bash -n "$F";done
squeue -u aryan.jha > "$C/campaign/initial-queue.txt"
G=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:10:00 --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --job-name=tadmm-persistent-gate --output="$C/logs/gate-%j.out" --error="$C/logs/gate-%j.err" "$C/campaign/gate.sh" "$C")
printf 'gate=%s\n' "$G" > "$C/campaign/submissions.txt"
for CAP in 30 10;do
 J=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=64 --hint=nomultithread --mem=128G --time=01:00:00 --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --dependency="afterok:$G" --kill-on-invalid-dep=yes --job-name="tadmm-large24-always$CAP" --output="$C/logs/cap$CAP-%j.out" --error="$C/logs/cap$CAP-%j.err" "$C/campaign/run.sh" "$C" "$CAP")
 A=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:15:00 --dependency="afterany:$J" --job-name="tadmm-always$CAP-energy" --output="$C/logs/audit-%j.out" --error="$C/logs/audit-%j.err" "$C/campaign/audit.sh" "$C" 24 "$J")
 printf 'cap%s=%s\ncap%s_audit=%s\n' "$CAP" "$J" "$CAP" "$A" >> "$C/campaign/submissions.txt"
done
cat "$C/campaign/submissions.txt"
