#!/usr/bin/env bash
# med3980 rated, no-battery certificates with derived bounds + lsec, directions given, 6 at a time (309 PC, 2026-11-03)
cd "$(dirname "$0")"
printf '%s\n' "$@" | xargs -P 6 -I{} env FB_DERIVED=1 FB_CUTS=kcl0,jbound,lsec SYS=med3980_s123 bash run.sh med_derived 1 farkas_bracketing.jl:C:{}
