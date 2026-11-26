#!/usr/bin/env bash
set -euo pipefail
ROOT=/home/aryan.jha/kamiak_jobs/stage_b_replay_oct15
cd "$ROOT/repo"
export PATH="$HOME/bin:$PATH"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export STAGE_B_CAPTURE_DIR="$ROOT/capture"
export FILTERDDP_PROCS_WORKERS=16 FILTERDDP_PROCS_GC=pass FILTERDDP_PROCS_REMOTE_PAIR=1
export PAIRED_WARMUP=0 PAIRED_ROUNDS=1
bash ddp/results/paired_pass/run_paired.sh large10kC_1ph 96 procs capture_stage_b > "$ROOT/logs/capture_launcher.log" 2>&1
grep STAGE_B_CAPTURE_COMPLETE ddp/results/parallel_in_time/logs/fddp_large10kC_1ph_T96_paired_capture_stage_b.log
test "$(find "$ROOT/capture" -name 'stage_*.jls' | wc -l)" -eq 96
sha256sum "$ROOT"/capture/stage_*.jls > "$ROOT/capture.sha256"
