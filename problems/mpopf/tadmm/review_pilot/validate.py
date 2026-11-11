"""Validate completion, original feasibility, matching warm-up and OpenDSS logs."""
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
summary = tomllib.loads((out / 'summary.toml').read_text())
meta = tomllib.loads((out / 'metadata.toml').read_text())
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
if 'no_iteration' in summary:
    near = [r for r in timed if float(r['original_primal']) <= summary['no_original_primal_tolerance']
        and abs(float(r['objective'])-summary['central_objective'])/max(1,abs(summary['central_objective'])) <= .001]
    assert near and int(near[0]['iteration']) == summary['no_iteration']
    assert summary['no_target_under_100'] == (summary['no_iteration'] < 100)
for name, rows in (('warmup', warm), ('trace', timed)):
    sub = list(csv.DictReader((out / f'{name}_subproblems.csv').open()))
    assert len(sub) == len(rows) * summary['T']
    pairs = {(int(s['iteration']), int(s['period'])) for s in sub}
    assert pairs == {(k, t) for k in range(1, len(rows) + 1) for t in range(1, summary['T'] + 1)}
    assert all(s['status'] in ('LOCALLY_SOLVED', 'OPTIMAL') for s in sub)
methods = ['central', 'tadmm']
if 'no_iteration' in summary:
    methods.append('tadmm_no')
for method in methods:
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
if meta.get('execution_mode') == 'distributed':
    layout = meta['worker_layout']
    assert len(layout) == meta['period_workers']
    assert len({r['pid'] for r in layout}) == len(layout)
    assert len({r['affinity'] for r in layout}) == len(layout)
    assert sorted(t for r in layout for t in r['periods']) == list(range(1, summary['T'] + 1))
    assert all(r['julia_threads'] == r['blas_threads'] == 1 and r['openmp_threads'] == '1' for r in layout)
    assert all(r['affinity'] != meta['coordinator_affinity'] for r in layout)
    checks['worker_layout_validation'] = 'pass'
    checks['peak_note'] = 'Sampled aggregate RSS of coordinator and workers, including setup/warm-up. RSS double-counts shared pages; consult PSS samples and Slurm cgroup peak.'
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
