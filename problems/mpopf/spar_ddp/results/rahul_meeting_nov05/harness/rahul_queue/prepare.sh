#!/usr/bin/env bash
set -euo pipefail
C=$PWD
BASE=/home/aryan.jha/kamiak_jobs/battery_count_oct21_campaign
FROM="$BASE/source"
S="$C/variants/10k_1020/source"
mkdir -p "$S/spar_ddp"
cp -a "$FROM/envs" "$FROM/rawData" "$S/"
for p in "$FROM"/spar_ddp/*; do
  [[ ${p##*/} == results ]] || cp -a "$p" "$S/spar_ddp/"
done
R="$S/spar_ddp/results/network_filterddp"; mkdir -p "$R"
F=network_data_large10kC_1ph_T96_periodic.jls
ln -s "$BASE/inputs/$F" "$R/$F"
cp schedule.jl "$S/spar_ddp/examples/power_system/meeting_grouped_schedule.jl"
python3 - "$S/spar_ddp/examples/power_system/procs_pass.jl" <<'PY'
from pathlib import Path
import sys
p = Path(sys.argv[1])
s = p.read_text()
a = s.index("    n = N + 1; h = 1\n")
b = s.index("        res = Vector{Any}", a)
start = s[a:b]
assert "for t in i:mid-1" in start
s = s[:a] + """    for pairs in meeting_scan_rounds(N, get(ENV, "SPAR_MEETING_GROUPED", "0") == "1")
""" + s[b:]
s = s.replace("        h *= 2\n    end\n    t_2 = time_ns()", "    end\n    t_2 = time_ns()", 1)
p.write_text('include(joinpath(@__DIR__, "meeting_grouped_schedule.jl"))\n' + '''
const MEETING_ROUND_CACHE = Dict{Tuple{Int,Bool},Vector{Vector{Tuple{Int,Int}}}}()
function meeting_scan_rounds(N, grouped)
    get!(MEETING_ROUND_CACHE, (N, grouped)) do
        replay_rounds(N, grouped ? 2 : 1)
    end
end
''' + s)
PY
(cd "$S"; find envs rawData spar_ddp \( -type f -o -type l \) -print0 | sort -z | xargs -0 sha256sum > "$C/source.sha256")
sha256sum job.sh prepare.sh schedule.jl > harness.sha256
