#!/usr/bin/env bash
set -euo pipefail
JOBS=/home/aryan.jha/kamiak_jobs
BASE="$JOBS/tadmm_review_nov10_large10k_t48_tuning_v3"
C="$JOBS/tadmm_review_nov11_large10k_interim"
LAUNCH="$JOBS/large10k-interim-nov11-launch"
P=problems/mpopf/tadmm/review_pilot
CASES=(60:20:64 72:20:76 48:20:76 60:exact:64 72:exact:76)
REUSE_PROFILES=${2:-}
if [[ $REUSE_PROFILES == ice ]]; then
 BASE=$C
 C="${C}_ice"
 CASES=(72:20:76 48:20:76 72:exact:76)
 test -s "$BASE/INTERIM_PROFILES_VALIDATED"
fi
test ! -e "$C"
mkdir -p "$C/campaign" "$C/logs"
cp -a "$BASE/source" "$BASE/julia-depot" "$C/"
cp "$BASE/input.sha256" "$C/campaign/base-input.sha256"
cp "$BASE/FULL_PERIOD_READY" "$C/"
cp "$BASE/campaign/audit_sweep_energy.jl" "$BASE/campaign/audit.sh" "$C/campaign/"
cp "$LAUNCH/audit-socp-storage-nov09.py" "$LAUNCH/summarize-tadmm-round-timing-nov10.py" "$C/campaign/"
cp "$LAUNCH/posthoc-med2522-caps-nov09.sh" "$C/campaign/posthoc.sh"
cp "$LAUNCH/setup-large10k-interim-horizons-nov11.sh" "$C/campaign/"
cd "$C/source"
sha256sum --check --quiet "$BASE/execution-source.sha256"
tar -xzf "$LAUNCH/horizons-overlay.tgz"
printf '%s\n' "$1" > source.commit
sha256sum --check --quiet "$C/campaign/base-input.sha256"
find "$P" test_networks problems/mpopf/spar_ddp/examples/power_system -type f ! -path '*/results/*' ! -path '*/__pycache__/*' -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256"
cat > "$C/campaign/export.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1
cd "$C/source"
P=problems/mpopf/tadmm/review_pilot
source /home/aryan.jha/kamiak_jobs/tadmm_review_nov09_sequential_dss_v2/campaign/env.sh
export JULIA_DEPOT_PATH="$C/julia-depot:$JULIA_DEPOT_PATH"
export PROFILE_PERIODIC=1 TERMINAL_SOC_SOFT=0 OPENDSS_ZERO_LINE_C=1
unset DROP_ROOT_DERS FILTERDDP_EXPORT_DIR
sha256sum --check --quiet "$C/execution-source.sha256"
INPUTS=problems/mpopf/spar_ddp/results/network_filterddp
EXPORT=problems/mpopf/spar_ddp/examples/power_system/export_ieee123c_data.jl
for T in 60 72; do
 test ! -e "$INPUTS/network_data_large10kC_1ph_T${T}_periodic.jls"
 julia --startup-file=no --project="$P" "$EXPORT" large10kC_1ph "$T" > "$C/logs/export-T$T.log" 2>&1
done
cat > "$C/campaign/check-profiles.jl" <<'JULIA'
using Serialization
folder = ARGS[1]
read_profile(t) = deserialize(joinpath(folder, "network_data_large10kC_1ph_T$(t)_periodic.jls"))
reference = read_profile(48)
for t in (60,72)
    d = read_profile(t)
    @assert d[:T] == t && d[:profile_periodic] && isapprox(t*d[:delta_t_h],24.0)
    for key in (:Nset,:Lset,:Bset,:Dset,:NLset,:B0_pu,:B_R_pu,:kVA_B,:kV_B)
        @assert d[key] == reference[key] key
    end
    common = gcd(t,48)
    for key in (:LoadShapeLoad,:LoadShapeCost)
        @assert length(d[key]) == t && all(isfinite,d[key])
        @assert all(isapprox.(d[key][1:t÷common:end],reference[key][1:48÷common:end];atol=1e-12))
    end
    @assert length(d[:LoadShapePV]) == t && all(x->0<=x<=1,d[:LoadShapePV]) && maximum(d[:LoadShapePV])>0.99
    println("T=",t," dt_h=",d[:delta_t_h]," common profile samples=",common," battery count=",length(d[:Bset]))
end
JULIA
julia --startup-file=no "$C/campaign/check-profiles.jl" "$INPUTS"
cat "$C/campaign/base-input.sha256" > "$C/input.sha256"
sha256sum "$INPUTS/network_data_large10kC_1ph_T60_periodic.jls" "$INPUTS/network_data_large10kC_1ph_T72_periodic.jls" >> "$C/input.sha256"
sha256sum --check --quiet "$C/input.sha256"
date -u +%FT%TZ > "$C/INTERIM_PROFILES_VALIDATED"
BASH
cat > "$C/campaign/run.sh" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
C=$1
T=$2
CAP=$3
MODEL=$4
[[ ${SLURM_JOB_NUM_NODES:?} == 1 ]]
test -s "$C/INTERIM_PROFILES_VALIDATED"
CPU_INFO=$(lscpu)
grep -Fq "$MODEL" <<< "$CPU_INFO"
export TADMM_CPU_MODEL_MATCH=$MODEL
export TADMM_RHO_POLICY=bootstrap_fixed TADMM_BOOTSTRAP_DECREASE=5 TADMM_BOOTSTRAP_GATE=Inf TADMM_RHO_FLOOR=24
export TADMM_RHO_PROVENANCE="Transfer from fastest large10k T48 tuning configuration"
if [[ $CAP == exact ]]; then
 export TADMM_LOCAL_SOLVE_POLICY=exact
else
 export TADMM_LOCAL_SOLVE_POLICY=persistent_cap TADMM_CAP_ITERATIONS=$CAP TADMM_CAP_SCREEN=finite
fi
exec bash "$C/source/problems/mpopf/tadmm/review_pilot/kamiak_full_period.sh" "$C" "$T" run primal_dual 15000 ma57 large10kC_1ph
BASH
for script in "$C/campaign/"*.sh; do bash -n "$script"; done
DEPENDENCY=()
if [[ $REUSE_PROFILES == ice ]]; then
 cp "$BASE/input.sha256" "$BASE/INTERIM_PROFILES_VALIDATED" "$C/"
 cp "$BASE/campaign/profile-job.txt" "$C/campaign/"
 sha256sum --check --quiet "$C/input.sha256"
else
 G=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:20:00 --job-name=large-interim-profiles --output="$C/logs/profiles-%j.out" --error="$C/logs/profiles-%j.err" "$C/campaign/export.sh" "$C")
 printf 'profiles=%s\n' "$G" > "$C/campaign/profile-job.txt"
 DEPENDENCY=(--dependency="afterok:$G" --kill-on-invalid-dep=yes)
fi
printf 'T\tcap\tcores\tmodel\tsolve\taudit\n' > "$C/campaign/submissions.tsv"
for CASE in "${CASES[@]}"; do
 IFS=: read -r T CAP CORES <<< "$CASE"
 if [[ $CORES == 64 ]]; then
  MODEL=6438Y+
  FEATURES='sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64'
 else
  MODEL=8368
  FEATURES='icelake&platinum-8368-2.40ghz&mem-1TB&cores-76'
 fi
 J=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task="$CORES" --hint=nomultithread --mem=128G --time=02:00:00 --constraint="$FEATURES" "${DEPENDENCY[@]}" --job-name="large$T-$CAP-$CORES" --output="$C/logs/solve-%j.out" --error="$C/logs/solve-%j.err" "$C/campaign/run.sh" "$C" "$T" "$CAP" "$MODEL")
 A=$(sbatch --parsable --account=dubey --partition=kamiak --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:20:00 --dependency="afterany:$J" --job-name="large$T-audit" --output="$C/logs/audit-%j.out" --error="$C/logs/audit-%j.err" "$C/campaign/posthoc.sh" "$J" "$C" "$T" 15000 large10kC_1ph)
 printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$T" "$CAP" "$CORES" "$MODEL" "$J" "$A" >> "$C/campaign/submissions.tsv"
done
cat "$C/campaign/profile-job.txt" "$C/campaign/submissions.tsv"
