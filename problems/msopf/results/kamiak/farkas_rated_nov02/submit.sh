#!/usr/bin/env bash
# Rated campaign (2026-11-02): every line at its rating. Per case: A (idle edges) -> C (idle certificates) and
# A -> AD (edges with batteries over the day) -> CD (their certificates), chained direction by direction (aftercorr).
set -euo pipefail
P="-p vcea -A dubey"
for spec in med3980_s123:16:24G:32G sfp1u_four:78:32G:48G; do
  IFS=: read sys n memA memAD <<< "$spec"; th=$([ $n -gt 16 ] && echo 12 || echo 16)
  a=$(sbatch --parsable $P -c 1 -t 12:00:00 --mem=$memA --array=1-$n%$th -J rA_$sys -o logs_rA_${sys}_%a.out farkas_job.sh $sys A)
  c=$(sbatch --parsable $P -c 1 -t 12:00:00 --mem=$memA --array=1-$n%$th -d aftercorr:$a -J rC_$sys -o logs_rC_${sys}_%a.out farkas_job.sh $sys C)
  ad=$(FB_SCRIPT=farkas_dispatched.jl sbatch --parsable $P -c 1 -t 24:00:00 --mem=$memAD --array=1-$n%8 -d aftercorr:$a -J rAD_$sys -o logs_rAD_${sys}_%a.out farkas_job.sh $sys AD)
  cd=$(FB_SCRIPT=farkas_dispatched.jl FB_ENV=sdp_pardiso FB_LDL=mkl NTH=8 FB_CUTS=kcl0,jbound,lsec sbatch --parsable $P -c 8 -t 24:00:00 --mem=80G --array=1-$n%4 -d aftercorr:$ad -J rCD_$sys -o logs_rCD_${sys}_%a.out farkas_job.sh $sys CD)
  echo "$sys: A $a C $c AD $ad CD $cd"
done
