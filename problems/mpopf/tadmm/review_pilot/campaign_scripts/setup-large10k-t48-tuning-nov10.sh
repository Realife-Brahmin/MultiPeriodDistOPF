#!/usr/bin/env bash
set -euo pipefail
JOBS=/home/aryan.jha/kamiak_jobs
BASE="$JOBS/tadmm_review_nov09_med2522_caps_sequence"
C="$JOBS/tadmm_review_nov10_large10k_t48_tuning_v3"
LAUNCH="$JOBS/large10k-t48-tuning-nov10-launch"
P=problems/mpopf/tadmm/review_pilot
test -s "$BASE/MED_SEQUENCE_GATE_VALIDATED"
test ! -e "$C"
mkdir -p "$C/campaign" "$C/logs"
cp -a "$BASE/source" "$BASE/julia-depot" "$C/"
cp "$BASE/input.sha256" "$BASE/FULL_PERIOD_READY" "$BASE/MED_SEQUENCE_GATE_VALIDATED" "$C/"
cp "$BASE/campaign/audit_sweep_energy.jl" "$BASE/campaign/audit.sh" "$C/campaign/"
cp "$LAUNCH/audit-socp-storage-nov09.py" "$LAUNCH/summarize-tadmm-round-timing-nov10.py" "$C/campaign/"
cp "$LAUNCH/posthoc-med2522-caps-nov09.sh" "$C/campaign/posthoc.sh"
cp "$LAUNCH/candidates.tsv" "$C/campaign/"
cp "$LAUNCH/setup-large10k-t48-tuning-nov10.sh" "$C/campaign/"
cd "$C/source"
sha256sum --check --quiet "$BASE/execution-source.sha256"
tar -xzf "$LAUNCH/controller-overlay.tgz"
printf '%s\n' "$1" > source.commit
sha256sum --check --quiet "$C/input.sha256"
find "$P" test_networks problems/mpopf/spar_ddp/examples/power_system -type f ! -path '*/results/*' ! -path '*/__pycache__/*' -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256"
cat > "$C/campaign/gate.jl" <<'JULIA'
include(abspath(ARGS[1]))
let rho = 120000.0, path = Float64[]
    for k in 1:26
        push!(path, rho)
        rho = bootstrap_fixed_rho(rho,k)
    end
    @assert path[1:5:26] == [120000,24000,4800,960,192,48]
end
@assert bootstrap_fixed_rho(120000,5;decrease=2,residual=2e-4,gate=1e-4)==120000
@assert bootstrap_fixed_rho(120000,5;decrease=2,residual=1e-4,gate=1e-4)==60000
@assert bootstrap_fixed_rho(60,5;floor=48,decrease=10,residual=0,gate=1e-4)==48
@assert bootstrap_fixed_rho(120000,4;decrease=2,residual=0,gate=1e-4)==120000
for path in ARGS[2:end]
    expr = Meta.parseall(read(path,String))
    @assert !any(x -> x isa Expr && x.head in (:error,:incomplete),expr.args)
end
println("Bootstrap schedule, defaults, residual gate and Julia parse checks passed.")
JULIA
/home/aryan.jha/opt/julia-1.12.4/bin/julia --startup-file=no "$C/campaign/gate.jl" "$P/legacy_rho.jl" "$P/run.jl" "$P/run_horizon_sweep.jl" > "$C/logs/controller-gate.log" 2>&1
cat > "$C/campaign/run.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1
CAP=$2
RHO=$3
[[ ${SLURM_JOB_NUM_NODES:?} == 1 ]]
export TADMM_RHO_POLICY=bootstrap_fixed TADMM_BOOTSTRAP_DECREASE=$4 TADMM_RHO_FLOOR=$5 TADMM_BOOTSTRAP_GATE=$6
export TADMM_RHO_PROVENANCE="Predeclared exploratory T48 tuning grid"
if [[ $CAP == exact ]]; then
 export TADMM_LOCAL_SOLVE_POLICY=exact
else
 export TADMM_LOCAL_SOLVE_POLICY=persistent_cap TADMM_CAP_ITERATIONS=$CAP TADMM_CAP_SCREEN=finite
fi
exec bash "$C/source/problems/mpopf/tadmm/review_pilot/kamiak_full_period.sh" "$C" 48 run primal_dual "$RHO" ma57 large10kC_1ph
BASH
for script in "$C/campaign/"*.sh; do bash -n "$script"; done
printf 'label\tcap\trhos\tdecrease\tfloor\tgate\tsolve\taudit\n' > "$C/campaign/submissions.tsv"
while IFS=$'\t' read -r LABEL CAP RHO DECREASE FLOOR GATE; do
 [[ $LABEL != label ]] || continue
 J=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=64 --hint=nomultithread --mem=128G --time=02:00:00 --constraint='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64' --job-name="l48-$LABEL" --output="$C/logs/solve-%j.out" --error="$C/logs/solve-%j.err" "$C/campaign/run.sh" "$C" "$CAP" "$RHO" "$DECREASE" "$FLOOR" "$GATE")
 A=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:20:00 --dependency="afterany:$J" --job-name="l48-audit-$LABEL" --output="$C/logs/audit-%j.out" --error="$C/logs/audit-%j.err" "$C/campaign/posthoc.sh" "$J" "$C" 48 "$RHO" large10kC_1ph)
 printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$LABEL" "$CAP" "$RHO" "$DECREASE" "$FLOOR" "$GATE" "$J" "$A" >> "$C/campaign/submissions.tsv"
done < "$C/campaign/candidates.tsv"
cat "$C/campaign/submissions.tsv"
