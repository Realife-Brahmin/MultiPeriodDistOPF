#!/usr/bin/env bash
set -euo pipefail
C=$SLURM_SUBMIT_DIR; cd "$C"
BASE=/home/aryan.jha/kamiak_jobs/battery_count_oct21_campaign
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$BASE/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_MAIN_FREE=1 STUDY_CORE_BUDGET=128
M="$C/jobs/$SLURM_JOB_ID"; mkdir "$M"
exec > >(tee "$M/stdout.log") 2> >(tee "$M/stderr.log" >&2)
{ date -u +%FT%TZ; hostname; julia --version; lscpu; taskset -pc $$; scontrol show job "$SLURM_JOB_ID"; } > "$M/hardware.txt"
grep -q '6760P' "$M/hardware.txt"
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER"; exit 2; }
sha256sum -c harness.sha256 > "$M/harness-check.log"
S="$C/variants/10k_1020/source"
(cd "$S"; sha256sum -c "$C/execution-source.sha256") > "$M/source-check.log"
export VALUE_REPLAY_NUMERICAL_ROOT="$S"
julia --startup-file=no unit.jl > "$M/unit.log" 2>&1
julia --startup-file=no parse.jl "$S/spar_ddp/examples/power_system/procs_pass.jl" "$S/spar_ddp/examples/power_system/br_integration.jl" > "$M/parse.log" 2>&1
run_one() {
 local T=$1 arm=$2 rep=$3 debug=$4
 export STUDY_SCAN=grouped STUDY_KERNEL=right STUDY_BLAS_THREADS=2
 if [[ $T == 24 ]]; then export STUDY_WORKERS=24; else export STUDY_WORKERS=48; fi
 if [[ $arm == base ]]; then export STUDY_BARRIER_REUSE=0; else export STUDY_BARRIER_REUSE=1; fi
 export STUDY_BARRIER_DEBUG=$debug
 local label="${arm}_${rep}" mapdir="$M/${T}_${arm}_${rep}"
 mkdir "$mapdir"
 (cd "$mapdir"; python3 "$C/affinity.py" "$STUDY_WORKERS" "$STUDY_BLAS_THREADS")
 export STUDY_CPU_MAP="$mapdir/cpu_groups.csv"
 export STUDY_OUT="$C/results/10k_1020_T${T}_ddp_r${label}"
 echo "BR_RUN_BEGIN T=$T arm=$arm repeat=$rep debug=$debug"
 bash "$S/spar_ddp/experiments/battery_count/run.sh" "$C" 10k_1020 "$T" ddp "$label"
 if [[ $debug == 1 ]]; then python3 gate.py "$STUDY_OUT" debug; else python3 gate.py "$STUDY_OUT"; fi
 sha256sum "$STUDY_OUT/solution.jls" >> "$M/solution.sha256"
 echo "BR_RUN_COMPLETE T=$T arm=$arm repeat=$rep debug=$debug"
}
compare() {
 local T=$1 ref=$2 cand=$3 name=$4
 mkdir -p "$C/comparisons_T${T}"
 julia --startup-file=no "$C/compare_solutions.jl" "$C/results/10k_1020_T${T}_ddp_r${ref}/solution.jls" "$C/results/10k_1020_T${T}_ddp_r${cand}/solution.jls" "$C/comparisons_T${T}/$name.toml" > "$C/comparisons_T${T}/$name.log" 2>&1
 grep -q 'strict_existing_match_gate = true' "$C/comparisons_T${T}/$name.toml"
}
# Expensive fingerprints/shadows only here, outside every clean clock.
run_one 24 debug 1 1
for T in 24 96; do
 run_one "$T" base 1 0
 run_one "$T" cache 1 0
 compare "$T" base_1 cache_1 cache_r1
 if [[ $T == 24 ]]; then compare 24 base_1 debug_1 debug; fi
 run_one "$T" cache 2 0
 run_one "$T" base 2 0
 compare "$T" base_2 cache_2 cache_r2
 compare "$T" base_1 base_2 reference_repeats
 date -u +%FT%TZ > "$C/T${T}.finished"
done
(cd "$S"; sha256sum -c "$C/execution-source.sha256") > "$M/source-check-after.log"
sha256sum -c harness.sha256 > "$M/harness-check-after.log"
date -u +%FT%TZ > "$C/fullsolve.finished"
