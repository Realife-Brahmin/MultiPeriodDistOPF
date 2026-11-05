#!/usr/bin/env bash
set -euo pipefail
C=$SLURM_SUBMIT_DIR; cd "$C/source"
OLD=/home/aryan.jha/kamiak_jobs/value_coupling_audit_oct22
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$OLD/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
exec > >(tee "$C/logs/stdout.log") 2> >(tee "$C/logs/stderr.log" >&2)
{ date -u +%FT%TZ; hostname;julia --version;lscpu;taskset -pc $$;scontrol show job "$SLURM_JOB_ID"; } > "$C/hardware.txt"
grep -q 'Gold 6438Y+' "$C/hardware.txt"
OTHER=$(squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%A' | sort -u | grep -vx "$SLURM_JOB_ID" || true)
[[ -z "$OTHER" ]] || { echo "OTHER_JOBS=$OTHER";exit 2; }
sha256sum -c "$C/source.sha256" > "$C/logs/source-check.log"
sha256sum -c "$C/runtime.sha256" > "$C/logs/runtime-check.log"
sha256sum -c "$C/input.sha256" > "$C/logs/input-check.log"
(cd "$OLD/source";sha256sum -c "$OLD/preparation-source.sha256") > "$C/logs/inherited-source-check.log"
julia --startup-file=no "$C/capture_unit.jl" > "$C/logs/unit.log" 2>&1
[[ $(df --output=avail -B1 "$C" | tail -1) -ge 32212254720 ]] || exit 2
OUT="$C/results/large10k_T24_influence_capture";test ! -e "$OUT";mkdir "$OUT"
while IFS= read -r K;do unset "$K";done < <(compgen -v | grep '^FILTERDDP_' || true)
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
export FILTERDDP_PROCS_WORKERS=24 FILTERDDP_WARMUP=0 FILTERDDP_TREE_KKT=1 FILTERDDP_NETWORK_BACKEND=tree
export FILTERDDP_STRUCTURED_DYNAMICS=1 FILTERDDP_TYPED_EQUATIONS=1 FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1
export FILTERDDP_CACHE_KKT_PATTERN=1 FILTERDDP_DIRECT_DIAG_HESSIAN=1 FILTERDDP_FACTOR_BACKED_POLICY=1
export FILTERDDP_LEAN_VALUE=1 FILTERDDP_LEAN_EXACT=1 FILTERDDP_HALF_BLOCK=1
export FILTERDDP_SCREEN=substation,vupper,ell,psubs FILTERDDP_LOADFLOW_START=1 FILTERDDP_AFFINE_LINESEARCH=1
export FILTERDDP_PROCS_GC=auto FILTERDDP_PROCS_BALANCED=1 FILTERDDP_BLOCKED_SOLVE=16
export FILTERDDP_MAX_ITERATIONS=150 FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1
export FILTERDDP_SCREEN_ABLATION_ORIGINAL_BOUNDS=1 FILTERDDP_OPTIMALITY_TOLERANCE=1e-5
export FILTERDDP_FROZEN_CAPTURE_DIR="$OUT/captures"
SOL=spar_ddp/results/network_filterddp/filterddp_solution_large10kC_1ph_T24_periodic_CBsystem.jls
test ! -e "$SOL"
timeout --signal=TERM --kill-after=30s 40m julia --startup-file=no --project=envs/spar_ddp spar_ddp/examples/power_system/procs_main.jl large10kC_1ph 24 solve > "$OUT/solve.log" 2>&1
mv "$SOL" "$OUT/solution.jls"
python3 spar_ddp/experiments/screening_ablation/validate.py "$OUT/solve.log" large10kC_1ph > "$OUT/numerical-check.log"
julia --startup-file=no --project=envs/tadmm spar_ddp/examples/power_system/opendss_check.jl large10kC_1ph 24 "$OUT/solution.jls" influence_capture > "$OUT/opendss.log" 2>&1
python3 spar_ddp/experiments/tree_ablation/validate_opendss.py "$OUT/opendss.log" > "$OUT/opendss-check.log"
julia --startup-file=no spar_ddp/experiments/value_group_solves/validate_physical.jl "$PWD" 24 "$OUT/solution.jls" large10kC_1ph > "$OUT/physical.log" 2>&1
julia --startup-file=no "$C/verify_capture.jl" "$OUT" > "$OUT/coverage-validation.log" 2>&1
env | grep -E '^(FILTERDDP_|OMP_NUM_THREADS|OPENBLAS_NUM_THREADS|JULIA_NUM_THREADS|REDUCED_|TERMINAL_SOC_)' | sort > "$OUT/settings.txt"
sha256sum "$OUT/solution.jls" > "$OUT/solution.sha256"
sha256sum -c "$C/runtime.sha256" > "$C/logs/runtime-check-after.log"
sha256sum -c "$C/source.sha256" > "$C/logs/source-check-after.log"
sha256sum -c "$C/input.sha256" > "$C/logs/input-check-after.log"
date -u +%FT%TZ > "$OUT/COMPLETE"
date -u +%FT%TZ > "$C/capture.finished"
df -h "$C" > "$C/disk-after.txt"

