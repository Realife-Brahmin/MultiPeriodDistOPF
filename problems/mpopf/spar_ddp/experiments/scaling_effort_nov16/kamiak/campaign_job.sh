#!/usr/bin/env bash
# Two-repeat matched runs on the October 16 campaign's protocol (run_case.sh of
# proper_logging_campaign_oct16) with the current DDP: chunked scan where it is
# shallower, separate primal/dual step lengths (DDP_SPLIT, default 1), balanced
# combinations on the rounds path, collector once a pass. Ipopt as in October:
# HSL MA57/MA97 on one thread, tol 1e-6, monotone, original-residual output.
# Both: periodic profile, per-system C_B and soft terminal SOC, screening
# substation,vupper,ell,psubs, shared load-flow start, full warm-up; original
# primal residual logged for DDP; OpenDSS check after every solve.
# Repeat 1 runs the arms in the given order, repeat 2 in reverse.
#   sbatch ... campaign_job.sh <system> "<T list>" "<arms: ddp ma57 ma97>"
set -uo pipefail
C=$SLURM_SUBMIT_DIR; cd "$C/repo" || exit 2
SYS=$1; TS=$2; ARMS=$3
export PATH=$HOME/opt/julia-1.12.7/bin:$PATH
export JULIA_DEPOT_PATH=/home/aryan.jha/kamiak_jobs/block_value_nov17/julia-depot
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE
unset NO_COLOR FORCE_COLOR OMP_PROC_BIND OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY LD_PRELOAD
M="$C/jobs/${SLURM_JOB_ID}_$SYS"; mkdir -p "$M"
{ date -Is; hostname; julia --version; lscpu; scontrol show job "$SLURM_JOB_ID"; } > "$M/hardware.txt"
grep -q '6438Y+' "$M/hardware.txt" || { echo "WRONG_NODE"; exit 2; }
squeue -h -w "$SLURMD_NODENAME" -t RUNNING,COMPLETING,SUSPENDED -o '%i %u %C' > "$M/node_jobs.txt"
P=problems/mpopf/spar_ddp; R=$P/results
mkdir -p $R/network_filterddp $R/parallel_in_time/logs $R/opendss_check $R/kkt_ordering/captures/solutions
export REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1
TOL=1e-5; [ "$SYS" = ieee123C_1ph ] && TOL=1e-6
REV=$(echo $ARMS | tr ' ' '\n' | tac | tr '\n' ' ')
for T in $TS; do
  D=$R/network_filterddp/network_data_${SYS}_T${T}_periodic.jls
  if [ ! -s "$D" ]; then
    PROFILE_PERIODIC=1 julia --startup-file=no --project=$P $P/examples/power_system/export_ieee123c_data.jl "$SYS" "$T" > "$M/export_T$T.log" 2>&1 || { echo "EXPORT_FAILED T=$T" >> "$M/summary.txt"; continue; }
  fi
  sha256sum "$D" >> "$M/inputs.sha256"
  W=$(( T < 63 ? T : 63 ))
  for REP in 1 2; do
    LIST=$ARMS; [ $REP = 2 ] && LIST=$REV
    for A in $LIST; do
      tag=camp_${SLURM_JOB_ID}_${A}_r$REP
      echo "BEGIN T=$T arm=$A repeat=$REP $(date -Is)" >> "$M/summary.txt"; t0=$(date +%s)
      if [ "$A" = ddp ]; then
        FILTERDDP_PROCS_WORKERS=$W FILTERDDP_PROCS_GC=pass FILTERDDP_PROCS_CHUNKED=1 FILTERDDP_PROCS_BALANCED=1 \
        FILTERDDP_SPLIT_STEP=${DDP_SPLIT:-1} FILTERDDP_FEASIBILITY_DIAGNOSTIC=1 FILTERDDP_SCREEN_ABLATION_ORIGINAL_BOUNDS=1 \
          PAIRED_WARMUP=full PAIRED_TOL=$TOL bash $P/examples/power_system/run_paired.sh "$SYS" "$T" procs "$tag" > "$M/run_T${T}_${A}_r$REP.out" 2>&1
        L=$R/parallel_in_time/logs/fddp_${SYS}_T${T}_paired_${tag}.log; O=$R/opendss_check/parsim_${SYS}_T${T}_paired_${tag}.log
        RES="$(grep -E 'solve complete' "$L" | tail -1) | $(grep -oE 'FilterDDP objective=[-0-9.eE+]+ max_equality_residual=[0-9.eE+-]+' "$L")"
      else
        if [ "$A" = ma97 ]; then lib=$HOME/opt/hsl-private/install/hsl_ma97-2.8.1-metis-openmp/lib/libhsl_ma97.so; extra=ma97_order=amd
        else lib=$HOME/opt/hsl-private/install/ma57-3.11.3-metis/lib/libma57.so; extra=ma57_pivot_order=2; fi
        L=$M/ipopt_T${T}_${A}_r$REP.log; SOL=$M/sol_T${T}_${A}_r$REP.jls
        IPOPT_SCREEN=substation,vupper,ell,psubs IPOPT_LOADFLOW_START=1 IPOPT_MATCH_DDP_START=1 IPOPT_TRACE=1 IPOPT_WARMUP=full \
          IPOPT_EXTRA_OPTIONS="linear_solver=$A;hsllib=$lib;linear_system_scaling=none;tol=1e-6;mu_strategy=monotone;inf_pr_output=original;$extra" \
          IPOPT_SAVE_DISPATCH=$SOL julia --startup-file=no --project=$P $P/examples/power_system/centralized_ipopt_matched.jl "$SYS" "$T" "$M/ipopt_T${T}_${A}_r${REP}_internal.log" > "$L" 2>&1
        bash $P/examples/power_system/run_opendss_check.sh "$SYS" "$T" "$SOL" "$tag" > "$M/opendss_T${T}_${A}_r$REP.out" 2>&1
        rm -f "$SOL"; O=$R/opendss_check/$tag.log
        RES="$(grep -oE 'status=[A-Z_]+ iterations=[0-9]+ objective=[-0-9.eE+]+ solve_time_s=[0-9.]+' "$L" | head -1)"
      fi
      cp "$L" "$M/" 2>/dev/null; cp "$O" "$M/" 2>/dev/null
      echo "END T=$T arm=$A repeat=$REP wall=$(( $(date +%s) - t0 ))s | $RES | $(grep -oE 'below limit: [0-9]+ of [0-9]+|above limit: [0-9]+' "$O" 2>/dev/null | tr '\n' ' ')" >> "$M/summary.txt"
    done
  done
done
date -Is > "$M/finished"
