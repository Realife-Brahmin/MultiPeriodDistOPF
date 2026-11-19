#!/usr/bin/env bash
set -euo pipefail
C=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_sequential_current_gate
OLD=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_sequential_dss_v2
test ! -e "$C";mkdir -p "$C/logs"
cp -a "$OLD/source" "$C/source"
tar -xzf /home/aryan.jha/kamiak_jobs/sequential-current-gate-overlay.tar.gz -C "$C/source"
cat > "$C/run.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_sequential_current_gate
source /home/aryan.jha/kamiak_jobs/tadmm_review_nov09_sequential_dss_v2/campaign/env.sh
cd "$C/source";P=problems/mpopf/tadmm/review_pilot
julia --startup-file=no --project="$P" "$P/validate_sequential_dss.jl"
cat > "$C/parse.jl" <<'JULIA'
for file in ARGS
 Meta.parseall(read(file,String));println("PARSE_OK ",file)
end
JULIA
julia --startup-file=no "$C/parse.jl" "$P/run_horizon_sweep.jl" "$P/run_sequential_dss.jl"
CASE=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_persistent_caps/results/run_T24_primal_dual_31665889
export SEQUENTIAL_DSS_GZIP=1
julia --startup-file=no --project="$P" "$P/run_sequential_dss.jl" "$CASE/instance.jls" "$CASE/central_dispatch.jls" "$C/central-model" model
test -s "$C/central-model/battery.csv.gz"
test ! -e "$C/central-model/battery.csv"
date -u +%FT%TZ > "$C/CURRENT_SEQUENCE_GATE_VALIDATED"
BASH
bash -n "$C/run.sh"
sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:15:00 --job-name=dss-current-unit-cli --output="$C/logs/gate-%j.out" --error="$C/logs/gate-%j.err" "$C/run.sh"
