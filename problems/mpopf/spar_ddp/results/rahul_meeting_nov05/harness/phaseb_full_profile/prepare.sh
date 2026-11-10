#!/usr/bin/env bash
set -euo pipefail
C=$PWD
BASE=/home/aryan.jha/kamiak_jobs/battery_count_oct21_campaign
(cd "$BASE"; sha256sum -c original-input.sha256)
(cd "$BASE/variants/10k_Y/source"; sha256sum -c ../execution-source.sha256; sha256sum -c ../input.sha256)
sha256sum job.sh profile_main.jl prepare.sh > harness.sha256
for V in 10k_Y 10k_1020; do
  FROM="$BASE/source"
  [[ $V != 10k_Y ]] || FROM="$BASE/variants/10k_Y/source"
  S="$C/variants/$V/source"
  mkdir -p "$S"
  cp -a "$FROM/envs" "$FROM/rawData" "$S/"
  mkdir "$S/spar_ddp"
  for p in "$FROM"/spar_ddp/*; do
    [[ ${p##*/} == results ]] || cp -a "$p" "$S/spar_ddp/"
  done
  R="$S/spar_ddp/results/network_filterddp"; mkdir -p "$R"
  for T in 96 144; do
    file="network_data_large10kC_1ph_T${T}_periodic.jls"
    input="$BASE/inputs/$file"
    [[ $V != 10k_Y ]] || input="$FROM/spar_ddp/results/network_filterddp/$file"
    ln -s "$input" "$R/$file"
  done
  driver="$S/spar_ddp/examples/power_system/procs_main.jl"
  python3 - "$driver" "$C/profile_main.jl" <<'PY'
from pathlib import Path
import sys
p, fragment = map(Path, sys.argv[1:])
source = p.read_text()
needle = "\nmain()\nrmprocs(workers())"
assert source.count(needle) == 1
p.write_text(source.replace(needle, '\nusing Profile\nif isempty(get(ENV, "PROFILE_RUN_DIR", ""))\n    main()\nelse\n' + fragment.read_text() + '\nend\nrmprocs(workers())'))
PY
  (cd "$S"; find envs rawData spar_ddp \( -type f -o -type l \) -print0 | sort -z | xargs -0 sha256sum > "$C/$V.sha256")
done
