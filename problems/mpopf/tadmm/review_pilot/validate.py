"""Validate completion, original feasibility, matching warm-up and OpenDSS logs."""
import csv
import hashlib
import json
from pathlib import Path
import re
import sys
import tomllib

out = Path(sys.argv[1]).resolve()
summary = tomllib.loads((out / 'summary.toml').read_text())
assert summary['converged']
assert summary['relative_objective_gap'] <= 1e-5
assert summary['tadmm_original_primal'] <= 5e-7
warm = list(csv.DictReader((out / 'warmup.csv').open()))
timed = list(csv.DictReader((out / 'trace.csv').open()))
assert len(warm) == len(timed) == summary['iterations']
for a, b in zip(warm, timed):
    for key in ('objective', 'original_primal', 'r_max', 's_max', 'rho'):
        assert abs(float(a[key]) - float(b[key])) <= 1e-10 * max(1, abs(float(a[key]))), (key, a, b)
checks = {}
for method in ('central', 'tadmm'):
    text = (out / f'{method}_opendss.log').read_text(encoding='utf-8-sig')
    assert 'converged=true' in text
    assert re.search(r'below limit: 0 of .*above limit: 0,', text)
    dv = float(re.search(r'voltage: largest difference ([\deE+.-]+) pu', text)[1])
    assert dv <= 1e-4, dv
    opf, dss = map(float, re.search(r'objective: optimizer ([\deE+.-]+), with OpenDSS import ([\deE+.-]+)', text).groups())
    assert abs(opf - dss) / max(1, abs(opf)) <= 1e-4
    checks[method] = dict(max_voltage_difference_pu=dv, objective_opf=opf, objective_opendss=dss)
samples = list(csv.DictReader((out / 'process_samples.csv').open()))
checks['sampled_process_peak_gib'] = max(int(s['working_set_bytes']) for s in samples) / 2**30
checks['peak_note'] = 'Sampled whole process, including warm-up and setup; not aggregate parallel-worker memory.'
checks['optimization_validation'] = 'pass'
checks['warmup_trace_validation'] = 'pass'
checks['opendss_validation'] = 'pass'
(out / 'validation.json').write_text(json.dumps(checks, indent=2) + '\n')
for row in (out / 'source.sha256').read_text().splitlines():
    expected, rel = row.split('  ', 1)
    assert hashlib.sha256((out / 'source' / rel).read_bytes()).hexdigest() == expected, rel
checks['source_snapshot_validation'] = 'pass'
(out / 'validation.json').write_text(json.dumps(checks, indent=2) + '\n')
print('CASE_VALIDATED', out)
