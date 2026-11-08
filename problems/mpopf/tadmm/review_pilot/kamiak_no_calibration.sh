#!/usr/bin/env bash
# One bounded candidate. No full warm-up, repeat, or clean timing claim.
set -euo pipefail
C=$(realpath "$1"); REF=$(realpath "$2"); RHO=${3:-384}; POLICY=${4:-relative}
[[ "$POLICY" == relative || "$POLICY" == fixed ]] || exit 2
cd "$C/source"
PILOT=problems/mpopf/tadmm/review_pilot
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
export JULIA_DEPOT_PATH="$C/julia-depot:$REF/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_parallel/julia-depot:/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/julia-depot:/home/aryan.jha/kamiak_jobs/tree_ablation_oct20_campaign_v2/julia-depot"
export PYTHONPATH=/home/aryan.jha/kamiak_jobs/tadmm_review_nov07_v2/python-deps${PYTHONPATH:+:$PYTHONPATH}
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export TADMM_LINEAR_SOLVER=ma57 TADMM_NWORKERS=24 TADMM_RHO_INIT=$RHO TADMM_RHO_POLICY=$POLICY
export TADMM_SCREEN=substation,vupper,ell
export TADMM_IPOPT_OPTIONS="hsllib=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so;linear_system_scaling=none;ma57_pivot_order=2"
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=0 OPENDSS_ZERO_LINE_C=1
unset NO_COLOR FORCE_COLOR OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
test -s "$C/DDP_TRANSFER_VALIDATED"
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
OUT="$C/results/large10k_T24_rho${RHO}_${POLICY}_NO_$SLURM_JOB_ID"
[[ ! -e "$OUT" ]] || exit 2
mkdir -p "$OUT"
CPU_RANGE=$(awk '/Cpus_allowed_list:/ {print $2}' /proc/self/status)
export TADMM_CPU_IDS=$(python3 -c 'import sys; a=[]
for p in sys.argv[1].split(","):
 b=p.split("-"); a.extend(range(int(b[0]),int(b[-1])+1))
assert len(a)>=25; print(",".join(map(str,a[:25])))' "$CPU_RANGE")
taskset -apc "$TADMM_CPU_IDS" "$$" > "$OUT/affinity.txt"
{ hostname; lscpu; scontrol show job "$SLURM_JOB_ID"; env | grep -E '^(TADMM_|OMP_|OPENBLAS_|JULIA_)'; } > "$OUT/allocation.txt"
grep -q '6438Y+' "$OUT/allocation.txt"
grep -q 'NumCPUs=64' "$OUT/allocation.txt"
python3 "$PILOT/capture.py" "$OUT" large10kC_1ph
sha256sum "$0" > "$OUT/orchestration.sha256"
REFERENCE="$REF/results/large10kC_1ph_T24_ma57_firstpass_31645790_a0"
taskset -c "$TADMM_CPU_IDS" julia --startup-file=no --project="$PILOT" "$PILOT/calibrate_no.jl" large10kC_1ph 24 "$OUT" "$REFERENCE" > "$OUT/solve.log" 2> "$OUT/solve.err" &
PID=$!
python3 "$PILOT/sample_linux.py" "$PID" "$OUT/process_samples.csv" > "$OUT/sampler.log" 2>&1 &
SAMPLER=$!
wait "$PID"
wait "$SAMPLER"
if [[ -f "$OUT/calibration_no_dispatch.jls" ]]; then
 julia --startup-file=no --project="$PILOT" problems/mpopf/spar_ddp/examples/power_system/opendss_check.jl large10kC_1ph 24 "$OUT/calibration_no_dispatch.jls" tadmm_no > "$OUT/no_opendss.log" 2>&1
fi
python3 "$PILOT/validate_no_calibration.py" "$OUT" > "$OUT/validation.log"
date -u +%FT%TZ > "$OUT/CALIBRATION_RECORDED"
