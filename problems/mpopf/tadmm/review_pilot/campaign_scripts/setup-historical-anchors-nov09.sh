#!/usr/bin/env bash
set -euo pipefail
OLD=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_targeted
C=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_historical_anchors
test ! -e "$C"
mkdir -p "$C/campaign" "$C/logs"
cp -a "$OLD/source" "$C/source"
cp -a "$OLD/julia-depot" "$C/julia-depot"
cp "$OLD/FULL_PERIOD_READY" "$C/"
cp "$OLD/campaign/audit_sweep_energy.jl" "$C/campaign/"
cp /home/aryan.jha/kamiak_jobs/historical-anchors-nov09-overlay.tar.gz "$C/campaign/overlay.tar.gz"
tar -xzf "$C/campaign/overlay.tar.gz" -C "$C/source"
printf '%s\n' "$1" > "$C/source/source.commit"
cd "$C/source"; P=problems/mpopf/tadmm/review_pilot
find "$P" test_networks problems/mpopf/spar_ddp/examples/power_system -type f ! -path '*/results/*' -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256"
bash -n "$P/kamiak_full_period.sh"
python3 -m py_compile "$P/validate_horizon_sweep.py"
cat > "$C/campaign/gate.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1; P=problems/mpopf/tadmm/review_pilot
cd "$C/source"
export JULIA_DEPOT_PATH="$C/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov08_horizon_sweep/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_accuracy_fix/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PROFILE_PERIODIC=1 OPENDSS_ZERO_LINE_C=1 TERMINAL_SOC_SOFT=0 DROP_ROOT_DERS=0
julia --startup-file=no --project="$P" problems/mpopf/spar_ddp/examples/power_system/export_ieee123c_data.jl ieee2522C_1ph 96 > "$C/campaign/export-med2522-T96.log" 2>&1
sha256sum problems/mpopf/spar_ddp/results/network_filterddp/network_data_large10kC_1ph_T24_periodic.jls problems/mpopf/spar_ddp/results/network_filterddp/network_data_ieee2522C_1ph_T96_periodic.jls > "$C/input.sha256"
for SOLVER in ma57 mumps; do
 export TADMM_LINEAR_SOLVER=$SOLVER
 if [[ $SOLVER == ma57 ]]; then
  export TADMM_IPOPT_OPTIONS="hsllib=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so;linear_system_scaling=none;ma57_pivot_order=2"
 else export TADMM_IPOPT_OPTIONS=linear_system_scaling=none; fi
 julia --startup-file=no --project="$P" "$P/validate_stage_reuse.jl" > "$C/campaign/gate-$SOLVER.log" 2>&1
done
sha256sum --check --quiet "$C/execution-source.sha256"
date -u +%FT%TZ > "$C/ANCHOR_GATE_VALIDATED"
BASH
cat > "$C/campaign/run.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1; SYSTEM=$2; T=$3; RHO=$4; SOLVER=$5
test -s "$C/ANCHOR_GATE_VALIDATED"
export TADMM_RHO_POLICY=bootstrap_fixed TADMM_RHO_FLOOR=48 TADMM_LOCAL_SOLVE_POLICY=exact
exec bash "$C/source/problems/mpopf/tadmm/review_pilot/kamiak_full_period.sh" "$C" "$T" run primal_dual "$RHO" "$SOLVER" "$SYSTEM"
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
G=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:15:00 --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --job-name=tadmm-anchor-gate --output="$C/logs/gate-%j.out" --error="$C/logs/gate-%j.err" "$C/campaign/gate.sh" "$C")
printf 'gate=%s\n' "$G" > "$C/campaign/submissions.txt"
DEP=$G
for SYSTEM in large10kC_1ph ieee2522C_1ph; do
 if [[ $SYSTEM == large10kC_1ph ]]; then T=24; N=1; RHO=15000; else T=96; N=2; RHO=20000; fi
 for SOLVER in ma57 mumps; do
  J=$(sbatch --parsable --account=dubey --partition=kamiak --nodes="$N" --ntasks="$N" --ntasks-per-node=1 --cpus-per-task=64 --hint=nomultithread --mem=128G --time=03:00:00 --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --dependency="afterok:$DEP" --kill-on-invalid-dep=yes --job-name="tadmm-$SYSTEM-T$T-$SOLVER" --output="$C/logs/$SYSTEM-T$T-$SOLVER-%j.out" --error="$C/logs/$SYSTEM-T$T-$SOLVER-%j.err" "$C/campaign/run.sh" "$C" "$SYSTEM" "$T" "$RHO" "$SOLVER")
  A=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:15:00 --dependency="afterany:$J" --job-name="tadmm-anchor-energy-T$T" --output="$C/logs/audit-%j.out" --error="$C/logs/audit-%j.err" "$C/campaign/audit.sh" "$C" "$T" "$J")
  printf '%s_T%s_%s=%s\naudit_%s=%s\n' "$SYSTEM" "$T" "$SOLVER" "$J" "$J" "$A" >> "$C/campaign/submissions.txt"
  DEP=$J
 done
done
cat "$C/campaign/submissions.txt"
