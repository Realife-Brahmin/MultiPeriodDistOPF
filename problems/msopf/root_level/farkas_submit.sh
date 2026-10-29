#!/usr/bin/env bash
# Submit Farkas bracketing cases on Kamiak as chains A (one task per direction) -> B -> C (one task per direction).
# Usage (from ~/msopf/jobs/<study>, beside farkas_job.sh): bash farkas_submit.sh [-d <jobid>] <case[:dirs]> ...
#   case   a farkas_bracketing.jl system (ieee123_s12, med3980_s134, ...); pairs have 2 directions, triples 16
#   :dirs  an array spec for phases A and C instead (e.g. med3980_two:1 reruns direction 1)
#   -d     start every chain after this job ends (afterany), e.g. to queue behind running work
# Memory: 24 GB (ieee123, med3980), 32 GB (sf6937). Prints one line per case with its three job ids.
set -euo pipefail
dep=""
if [ "${1:-}" = "-d" ]; then dep="-d afterany:$2"; shift 2; fi
for arg in "$@"; do
  sys=${arg%%:*}; dirs=${arg#*:}
  if [ "$dirs" = "$arg" ]; then
    digits=${sys##*_s}
    case $sys in *_two) n=2;; *_three) n=16;; *) [ ${#digits} -eq 2 ] && n=2 || n=16;; esac
    dirs=1-$n
  fi
  case $sys in sfp1u*) mem=32G;; *) mem=24G;; esac
  P="-p vcea -A dubey -c 1 -t 10:00:00 --mem=$mem"
  a=$(sbatch --parsable $P $dep --array=$dirs -J fbA_$sys -o logs_fbA_${sys}_%a.out farkas_job.sh $sys A)
  b=$(sbatch --parsable $P -d afterany:$a -J fbB_$sys -o logs_fbB_$sys.out farkas_job.sh $sys B)
  c=$(sbatch --parsable $P -d afterany:$b --array=$dirs -J fbC_$sys -o logs_fbC_${sys}_%a.out farkas_job.sh $sys C)
  echo "$sys dirs $dirs: A $a B $b C $c"
done
