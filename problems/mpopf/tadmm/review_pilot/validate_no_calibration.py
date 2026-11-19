"""Validate a capped NO calibration without promoting it to a timed benchmark."""
import csv
import hashlib
import json
from pathlib import Path
import re
import sys
try:
    import tomllib
except ImportError:
    import tomli as tomllib

out = Path(sys.argv[1]).resolve()
meta = tomllib.loads((out / 'calibration.toml').read_text())
trace = list(csv.DictReader((out / 'calibration.csv').open()))
assert len(trace) == meta['iterations_executed'] <= 100
near = [r for r in trace if float(r['original_primal']) <= meta['no_original_tolerance']
    and abs(float(r['objective'])-meta['reference_objective'])/max(1,abs(meta['reference_objective'])) <= .001]
assert meta['no_iteration'] == (int(near[0]['iteration']) if near else 0)
assert meta['target_under_100'] == (0 < meta['no_iteration'] < 100)
sub = list(csv.DictReader((out / 'calibration_subproblems.csv').open()))
assert len(sub) == len(trace)*meta['T']
assert {(int(r['iteration']),int(r['period'])) for r in sub} == {(k,t) for k in range(1,len(trace)+1) for t in range(1,meta['T']+1)}
assert all(r['status'] in ('LOCALLY_SOLVED','OPTIMAL') for r in sub)
for line in (out / 'source.sha256').read_text().splitlines():
    digest, rel = line.split('  ',1)
    assert hashlib.sha256((out / 'source' / rel).read_bytes()).hexdigest() == digest
result = dict(kind='100-step NO calibration', no_iteration=meta['no_iteration'],
    target_under_100=meta['target_under_100'], strict_converged=meta['strict_converged'],
    benchmark_timing_claim=False, stage_status_validation='pass', source_validation='pass')
if near:
    text = (out / 'no_opendss.log').read_text(encoding='utf-8-sig')
    assert 'converged=true' in text
    assert re.search(r'below limit: 0 of .*above limit: 0,',text)
    dv = float(re.search(r'voltage: largest difference ([\deE+.-]+) pu',text)[1])
    assert dv <= 1e-4
    opf,dss = map(float,re.search(r'objective: optimizer ([\deE+.-]+), with OpenDSS import ([\deE+.-]+)',text).groups())
    assert abs(opf-dss)/max(1,abs(opf)) <= 1e-4
    assert meta['no_original_primal'] <= meta['no_original_tolerance']
    assert meta['no_objective_gap'] <= .001
    result.update(no_dispatch_opendss_validation='pass',no_max_voltage_difference_pu=dv)
    (out / 'NO_POINT_VALIDATED').write_text('Actual first NO dispatch independently validated; not a full-convergence benchmark.\n')
else:
    result['no_dispatch_opendss_validation'] = 'not applicable: NO not reached'
(out / 'calibration_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
