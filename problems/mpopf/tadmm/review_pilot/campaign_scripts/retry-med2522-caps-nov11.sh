#!/usr/bin/env bash
set -euo pipefail
root=/home/aryan.jha/kamiak_jobs
caps="$root/tadmm_review_nov09_med2522_caps_sequence"
exact="$root/tadmm_review_nov09_med2522_single_node"
tools="$root/tadmm_review_nov11_large10k_interim/campaign"
p=problems/mpopf/tadmm/review_pilot
mode=$1
if [[ $mode == solve ]]; then
  cap=$2
  cd "$caps/source"
  sha256sum --check --quiet "$caps/execution-source.sha256"
  test -s "$caps/MED_SEQUENCE_GATE_VALIDATED"
  [[ $SLURM_JOB_NUM_NODES == 1 && $SLURM_CPUS_PER_TASK == 128 ]]
  lscpu | grep '6760P'
  export TADMM_RHO_POLICY=bootstrap_fixed TADMM_RHO_FLOOR=48
  export TADMM_LOCAL_SOLVE_POLICY=persistent_cap TADMM_CAP_ITERATIONS=$cap TADMM_CAP_SCREEN=finite
  # Replace only the obsolete CPU guard in a transient copy of the verified launcher.
  launcher=$(mktemp)
  trap 'rm -f -- "$launcher"' EXIT
  python3 - "$p/kamiak_full_period.sh" "$launcher" <<'PY'
import sys
from pathlib import Path
original = Path(sys.argv[1]).read_text()
old = "grep -q '6438Y+'"
assert original.count(old) == 1
Path(sys.argv[2]).write_text(original.replace(old, "grep -q '6760P'"))
PY
  sha256sum "$launcher"
  bash "$launcher" "$caps" 96 run primal_dual 20000 ma57 ieee2522C_1ph
elif [[ $mode == audit ]]; then
  c=$2
  j=$3
  run="$c/results/run_T96_primal_dual_$j"
  report="$c/results/chronological-$j-$SLURM_JOB_ID"
  test -s "$run/CASE_RECORDED"
  mkdir "$report"
  cd "$caps/source"
  sha256sum --check --quiet "$caps/execution-source.sha256"
  sha256sum --check --quiet "$caps/input.sha256"
  cmp "$caps/input.sha256" "$c/input.sha256"
  cp "$c/input.sha256" "$c/execution-source.sha256" "$report/"
  sha256sum "$run/instance.jls" "$run/central_dispatch.jls" "$run/rho20000/endpoint.jls" > "$report/replay-inputs.sha256"
  source "$root/tadmm_review_nov09_sequential_dss_v2/campaign/env.sh"
  export SEQUENTIAL_DSS_GZIP=1
  export PYTHONPATH="$root/tadmm_review_nov07_v2/python-deps${PYTHONPATH:+:$PYTHONPATH}"
  for name in central endpoint; do
    if [[ $name == central ]]; then dispatch="$run/central_dispatch.jls"; else dispatch="$run/rho20000/endpoint.jls"; fi
    julia --startup-file=no --project="$p" "$p/run_sequential_dss.jl" "$run/instance.jls" "$dispatch" "$report/$name" model
    python3 "$tools/audit-socp-storage-nov09.py" "$report/$name" test_networks/ieee2522C_1ph/Storage.dss
  done
  mkdir "$report/rounds"
  cp "$run/rho20000/trace.csv" "$run/rho20000/trace_subproblems.csv" "$report/rounds/"
  python3 "$tools/summarize-tadmm-round-timing-nov10.py" "$report/rounds"
  sacct -j "$j,$SLURM_JOB_ID" --format=JobID,State,Elapsed,ExitCode,MaxRSS,MaxVMSize,ReqMem,AllocCPUS,NodeList -P > "$report/slurm.txt"
  sha256sum --check --quiet "$report/replay-inputs.sha256"
  date -u +%FT%TZ > "$report/REPLAY_RECORDED"
else
  exit 2
fi
