#!/usr/bin/env bash
set -euo pipefail
C=${1:?July campaign directory required}
cd "$C/source"
P=problems/mpopf/tadmm/review_pilot
sha256sum --check --quiet "$C/execution-source.sha256"
sha256sum --check --quiet "$C/input.sha256"
cp -n "$C/execution-source.sha256" "$C/campaign/pre-flat-start-source.sha256"
cp -n "$P/model.jl" "$C/campaign/model-before-flat-start.jl"
python3 - "$P/model.jl" <<'PY'
import sys
from pathlib import Path

model = Path(sys.argv[1])
original = model.read_text()
loadflow = "        vv,pp,qq,ll,_=loadflow_initial_values(d,t)"
flat = """        # The archived 4.1569 kV input defeated the idle-battery sweep.
        # Use a flat network start for both solvers on this diagnostic.
        vv = fill(1.05^2, length(ns))
        pp = zeros(length(es))
        qq = zeros(length(es))
        ll = fill(0.01 / 1.05^2, length(es))"""
assert original.count(loadflow) == 1
model.write_text(original.replace(loadflow, flat))
PY
find "$P" test_networks problems/mpopf/spar_ddp/examples/power_system -type f ! -path '*/results/*' ! -path '*/__pycache__/*' -print0 | sort -z | xargs -0 sha256sum > "$C/execution-source.sha256"
sha256sum --check --quiet "$C/execution-source.sha256"
