#!/usr/bin/env bash
# Kamiak: how a period solve slows when k copies run at once on one node, Gurobi vs Ipopt.
# usage: gurobi_contention.sh C run_dir   (BENCH_SET=3, SC_CONIC_CONE as exported)
set -euo pipefail
C=$(realpath "$1"); R=$2; P=problems/mpopf/tadmm/review_pilot
cd "$C/source"
export PATH=/home/aryan.jha/opt/julia-1.12.4/bin:$PATH
D=/home/aryan.jha/kamiak_jobs
export JULIA_DEPOT_PATH="$C/julia-depot:$D/tadmm_review_nov11_large10k_interim_ice/julia-depot:$D/tadmm_review_nov08_horizon_sweep/julia-depot:$D/tadmm_review_nov08_watchdog_rho/julia-depot:$D/tadmm_review_nov07_accuracy_fix/julia-depot:$D/tadmm_review_nov07_parallel/julia-depot:$D/tadmm_review_nov07_v2/julia-depot:$D/tree_ablation_oct20_campaign_v2/julia-depot"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
unset NO_COLOR FORCE_COLOR LD_PRELOAD
LIB=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so
export LD_LIBRARY_PATH=$HOME/opt/hsl-private/deps${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export TADMM_LINEAR_SOLVER=ma57 TADMM_IPOPT_OPTIONS="hsllib=$LIB;linear_system_scaling=none;ma57_pivot_order=2"
export TADMM_CB_PROFILE=july REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=0 SC_ETA=0.95
export GRB_LICENSE_FILE=/opt/apps/gurobi/13.0.1/licenses/gurobi.lic SC_STAGE_SOLVER=gurobi BENCH_SET=3 TADMM_SCREEN=substation,vupper,ell
grep -m1 "Model name" /proc/cpuinfo; nproc
CPUS=($(python3 -c 'import os; print(*sorted(os.sched_getaffinity(0)))'))
for k in 1 8 24 48; do
 [ $k -le ${#CPUS[@]} ] || continue
 tmp=$(mktemp -d)
 for i in $(seq 0 $((k-1))); do
  t=$(( (i % 48) + 1 ))
  taskset -c ${CPUS[$i]} julia --startup-file=no --project="$P/conic_env" "$P/gurobi_stage_bench.jl" "$R/instance.jls" "$R/lean_trace_timed.csv" $t > $tmp/$i.out 2>&1 &
 done
 wait
 for s in tol1e-8 ipopt_warm; do
  grep -h "setting=$s " $tmp/*.out | sed -E 's/.*wall_mean=([0-9.]+).*/\1/' | awk -v k=$k -v s=$s '{n++; x+=$1; if($1>m)m=$1} END{printf "CONTENTION k=%d setting=%s n=%d mean_wall=%.3f max_wall=%.3f\n",k,s,n,x/n,m}'
 done
 rm -rf $tmp
done
