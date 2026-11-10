#!/usr/bin/env bash
set -euo pipefail
OLD=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_persistent_caps
C=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_sequential_dss
test ! -e "$C"
mkdir -p "$C/campaign" "$C/logs"
cp -a "$OLD/source" "$C/source"
cp /home/aryan.jha/kamiak_jobs/sequential-dss-nov09-overlay.tar.gz "$C/campaign/overlay.tar.gz"
tar -xzf "$C/campaign/overlay.tar.gz" -C "$C/source"
printf '%s\n' "$1" > "$C/source/source.commit"
cd "$C/source";P=problems/mpopf/tadmm/review_pilot
find "$P" test_networks problems/mpopf/spar_ddp/examples/power_system -type f ! -path '*/results/*' ! -path '*/__pycache__/*' -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256"
cat > "$C/campaign/env.sh" <<'BASH'
export JULIA_DEPOT_PATH="/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_persistent_caps/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov08_horizon_sweep/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
BASH
cat > "$C/campaign/gate.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1;P=problems/mpopf/tadmm/review_pilot;cd "$C/source";source "$C/campaign/env.sh"
sha256sum --check --quiet "$C/execution-source.sha256"
julia --startup-file=no --project="$P" "$P/validate_sequential_dss.jl"
date -u +%FT%TZ > "$C/SEQUENTIAL_GATE_VALIDATED"
BASH
cat > "$C/campaign/run.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1;P=problems/mpopf/tadmm/review_pilot;cd "$C/source";source "$C/campaign/env.sh"
test -s "$C/SEQUENTIAL_GATE_VALIDATED"
sha256sum --check --quiet "$C/execution-source.sha256"
OLD=/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_persistent_caps
for NAME in central cap30 cap10;do
 if [[ $NAME == central ]]; then
  CASE="$OLD/results/run_T24_primal_dual_31665889";SOL="$CASE/central_dispatch.jls"
 elif [[ $NAME == cap30 ]]; then
  CASE="$OLD/results/run_T24_primal_dual_31665887";SOL="$CASE/rho15000/endpoint.jls"
 else
  CASE="$OLD/results/run_T24_primal_dual_31665889";SOL="$CASE/rho15000/endpoint.jls"
 fi
 for EFF in model dss;do
  OUT="$C/results/$NAME-$EFF"
  julia --startup-file=no --project="$P" "$P/run_sequential_dss.jl" "$CASE/instance.jls" "$SOL" "$OUT" "$EFF" > "$C/logs/$NAME-$EFF.log" 2> "$C/logs/$NAME-$EFF.err"
  for F in "$OUT/"*.csv;do gzip "$F";done
 done
done
date -u +%FT%TZ > "$C/SEQUENTIAL_REPLAYS_RECORDED"
BASH
for F in "$C/campaign/"*.sh;do bash -n "$F";done
squeue -u aryan.jha > "$C/campaign/initial-queue.txt"
G=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=3G --time=00:10:00 --job-name=dss-sequence-gate --output="$C/logs/gate-%j.out" --error="$C/logs/gate-%j.err" "$C/campaign/gate.sh" "$C")
J=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=01:00:00 --dependency="afterok:$G" --kill-on-invalid-dep=yes --job-name=dss-sequence-replays --output="$C/logs/run-%j.out" --error="$C/logs/run-%j.err" "$C/campaign/run.sh" "$C")
printf 'gate=%s\nreplays=%s\n' "$G" "$J" > "$C/campaign/submissions.txt"
cat "$C/campaign/submissions.txt"
