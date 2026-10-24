#!/usr/bin/env bash
set -euo pipefail
C=$(realpath "$1"); V=$2; T=$3; METHOD=$4; REP=$5
cd "$C/variants/$V/source"
OUT="$C/results/${V}_T${T}_${METHOD}_r${REP}"
[[ ! -e "$OUT" ]] || exit 2
mkdir -p "$OUT"
while IFS= read -r NAME; do unset "$NAME"; done < <(compgen -v | grep -E '^(FILTERDDP_|IPOPT_)')
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_DYNAMIC=FALSE
unset OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD NO_COLOR FORCE_COLOR
START=$(date +%s)
printf 'BATTERY_COUNT_RUN variant=%s T=%s method=%s repeat=%s\n' "$V" "$T" "$METHOD" "$REP" > "$OUT/solve.log"
if [[ $METHOD == ddp ]]; then
 export FILTERDDP_PROCS_WORKERS=$((T<63?T:63)) FILTERDDP_WARMUP=full
 export FILTERDDP_TREE_KKT=1 FILTERDDP_NETWORK_BACKEND=tree FILTERDDP_STRUCTURED_DYNAMICS=1 FILTERDDP_TYPED_EQUATIONS=1
 export FILTERDDP_TRIPLET_SECOND_DERIVATIVES=1 FILTERDDP_CACHE_KKT_PATTERN=1 FILTERDDP_DIRECT_DIAG_HESSIAN=1
 export FILTERDDP_FACTOR_BACKED_POLICY=1 FILTERDDP_LEAN_VALUE=1 FILTERDDP_LEAN_EXACT=1 FILTERDDP_HALF_BLOCK=1
 export FILTERDDP_SCREEN=substation,vupper,ell,psubs FILTERDDP_LOADFLOW_START=1 FILTERDDP_AFFINE_LINESEARCH=1
 export FILTERDDP_PROCS_GC=pass FILTERDDP_PROCS_BALANCED=1 FILTERDDP_BLOCKED_SOLVE=16
 export FILTERDDP_MAX_ITERATIONS=400 FILTERDDP_OPTIMALITY_TOLERANCE=1e-5
 export FILTERDDP_TIMING_DIAGNOSTIC=1 FILTERDDP_FEASIBILITY_DIAGNOSTIC=1 FILTERDDP_SCREEN_ABLATION_ORIGINAL_BOUNDS=1
 SOL=spar_ddp/results/network_filterddp/filterddp_solution_large10kC_1ph_T${T}_periodic_CBsystem.jls
 [[ ! -e "$SOL" ]] || { echo 'Unclaimed dispatch exists';exit 2; }
 timeout --signal=TERM --kill-after=30s 50m julia --startup-file=no --project=envs/spar_ddp spar_ddp/examples/power_system/procs_main.jl large10kC_1ph "$T" solve >> "$OUT/solve.log" 2>&1
 mv "$SOL" "$OUT/solution.jls"
 python3 spar_ddp/experiments/screening_ablation/validate.py "$OUT/solve.log" large10kC_1ph
else
 export IPOPT_SCREEN=substation,vupper,ell,psubs IPOPT_LOADFLOW_START=1 IPOPT_MATCH_DDP_START=1 IPOPT_TRACE=1 IPOPT_WARMUP=full
 if [[ $METHOD == ma97 ]]; then
  LIB=/home/aryan.jha/opt/hsl-private/install/hsl_ma97-2.8.1-metis-openmp/lib/libhsl_ma97.so; EXTRA=ma97_order=amd
 elif [[ $METHOD == ma57 ]]; then
  LIB=/home/aryan.jha/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so; EXTRA=ma57_pivot_order=2
 else exit 2; fi
 sha256sum "$LIB" > "$OUT/private-library.sha256"
 export IPOPT_EXTRA_OPTIONS="linear_solver=$METHOD;hsllib=$LIB;linear_system_scaling=none;tol=1e-6;mu_strategy=monotone;inf_pr_output=original;$EXTRA"
 export IPOPT_SAVE_DISPATCH="$OUT/solution.jls"
 timeout --signal=TERM --kill-after=30s 50m julia --startup-file=no --project=envs/spar_ddp spar_ddp/examples/power_system/centralized_ipopt_matched.jl large10kC_1ph "$T" "$OUT/ipopt_internal.log" >> "$OUT/solve.log" 2>&1
 python3 "$C/check_ipopt_logs.py" "$OUT/solve.log" "$OUT/ipopt_internal.log" > "$OUT/logging_check.json"
 grep -q 'IPOPT_WARMUP_CHECK' "$OUT/solve.log"
fi
julia --startup-file=no --project=envs/tadmm spar_ddp/examples/power_system/opendss_check.jl large10kC_1ph "$T" "$OUT/solution.jls" "${V}_${METHOD}" > "$OUT/opendss.log" 2>&1
python3 spar_ddp/experiments/tree_ablation/validate_opendss.py "$OUT/opendss.log"
julia --startup-file=no spar_ddp/experiments/battery_count/validate_physical.jl "$PWD" "$T" "$OUT/solution.jls" > "$OUT/physical.log" 2>&1
python3 spar_ddp/experiments/battery_count/extract.py "$OUT" "$METHOD"
env | grep -E '^(FILTERDDP_|IPOPT_|JULIA_NUM_THREADS=|OPENBLAS_NUM_THREADS=|OMP_NUM_THREADS=|REDUCED_|TERMINAL_SOC_)' | LC_ALL=C sort > "$OUT/settings.txt"
echo "total_process_wall_s=$(($(date +%s)-START))" > "$OUT/wall.txt"
date -u +%FT%TZ > "$OUT/COMPLETE"
