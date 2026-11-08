"""Validate July-weight assessments while keeping truncated NLPs explicit."""
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys
try:
    import tomllib
except ImportError:
    import tomli as tomllib

out = Path(sys.argv[1]).resolve()
stage = sys.argv[2]
assert stage in ('reference', 'exact', 'cap30', 'complete')
meta = tomllib.loads((out / 'comparison.toml').read_text())
assert meta['cb_profile'] == 'july' and meta['terminal_soc_penalty'] == 0
assert meta['linear_solver'] == 'ma57' and meta['T'] == 24
assert math.isclose(meta['C_B'], 8.013987384856764e-8, rel_tol=1e-14)
assert meta['no_objective_tolerance'] == .001
assert meta['no_original_tolerance'] == (1e-6 if meta['system'] == 'ieee123C_1ph' else 1e-5)
assert meta['reference_original_primal'] <= 1e-6
for line in (out / 'source.sha256').read_text().splitlines():
    digest, rel = line.split('  ', 1)
    assert hashlib.sha256((out / 'source' / rel).read_bytes()).hexdigest() == digest


def dss(log):
    text = log.read_text(encoding='utf-8-sig')
    assert 'converged=true' in text
    assert re.search(r'below limit: 0 of .*above limit: 0,', text)
    dv = float(re.search(r'voltage: largest difference ([\deE+.-]+) pu', text)[1])
    assert dv <= 1e-4
    opf, actual = map(float, re.search(r'objective: optimizer ([\deE+.-]+), with OpenDSS import ([\deE+.-]+)', text).groups())
    assert abs(opf-actual)/max(1, abs(opf)) <= 1e-4
    return dict(max_voltage_difference_pu=dv, objective_opf=opf, objective_opendss=actual)


reference_dss = dss(out / 'central_opendss.log')
assert math.isclose(reference_dss['objective_opf'], meta['reference_objective'], rel_tol=1e-10, abs_tol=1e-6)
result = dict(kind='July-weight exact/cap30 assessment', benchmark_timing_claim=False,
    reference_opendss=reference_dss, C_B=meta['C_B'], source_validation='pass', arms={})
if stage == 'reference':
    (out / 'REFERENCE_VALIDATED').write_text('Exact July-weight reference and independent OpenDSS checks passed.\n')
else:
    policies = tuple(meta.get('policies',('exact','cap30'))) if stage == 'complete' else (stage,)
    for policy in policies:
        arm = out / policy
        row = tomllib.loads((arm / 'assessment.toml').read_text())
        assert row['policy'] == policy
        trace = list(csv.DictReader((arm / 'trace.csv').open()))
        assert len(trace) == row['iterations'] <= meta['maximum_iterations']
        limited = 0
        for phase in ('warmup', 'trace'):
            outer = list(csv.DictReader((arm / f'{phase}.csv').open()))
            if phase == 'warmup':
                assert 0 < len(outer) <= meta['warmup_steps_per_arm']
            sub = list(csv.DictReader((arm / f'{phase}_subproblems.csv').open()))
            assert len(sub) == len(outer)*meta['T']
            assert {(int(s['iteration']), int(s['period'])) for s in sub} == {
                (k,t) for k in range(1,len(outer)+1) for t in range(1,meta['T']+1)}
            for s in sub:
                cutoff = s['status'] == 'ITERATION_LIMIT'
                assert s['iteration_limit_hit'] == str(cutoff).lower()
                assert s['status'] in ('LOCALLY_SOLVED','OPTIMAL') or (policy == 'cap30' and cutoff)
                if cutoff:
                    checked_local = meta.get('local_primal_report',True)
                    assert s.get('local_primal_evaluated',str(checked_local).lower()) == str(checked_local).lower()
                    if checked_local:
                        assert math.isfinite(float(s['local_primal'])) and float(s['local_primal']) >= 0
                    else:
                        assert math.isnan(float(s['local_primal']))
                    limited += (phase == 'trace')
                if policy == 'cap30':
                    assert int(s['ipopt_iterations']) <= 30
        qualified = [r for r in trace if float(r['original_primal']) <= meta['no_original_tolerance'] and
            abs(float(r['objective'])-meta['reference_objective'])/max(1,abs(meta['reference_objective'])) <= .001]
        first = int(qualified[0]['iteration']) if qualified else 0
        assert row['no_iteration'] == first
        assert row['target_under_100'] == (0 < first < 100)
        last = trace[-1]
        for key in ('objective','original_primal','r_max','s_max'):
            assert math.isclose(row[key], float(last[key]), rel_tol=1e-12, abs_tol=1e-10)
        layout = row['worker_layout']
        assert len(layout) == meta['T']
        assert sorted(t for worker in layout for t in worker['periods']) == list(range(1,meta['T']+1))
        assert len({worker['affinity'] for worker in layout}) == meta['T']
        assert all(w['local_solve_policy'] == policy and w['julia_threads'] == w['blas_threads'] == 1 and
            w['openmp_threads'] == '1' for w in layout)
        assert all(w.get('local_primal_report',True)==meta.get('local_primal_report',True) for w in layout)
        checked = dict(no_iteration=first, target_under_100=row['target_under_100'],
            truncated_local_solves=limited, stage_status_validation='pass', original_primal_at_end=row['original_primal'])
        if first:
            assert math.isclose(row['no_original_primal'], float(qualified[0]['original_primal']), abs_tol=1e-10, rel_tol=1e-12)
            assert math.isclose(row['no_objective'], float(qualified[0]['objective']), abs_tol=1e-9, rel_tol=1e-12)
            checked['no_opendss'] = dss(arm / 'no_opendss.log')
            assert math.isclose(checked['no_opendss']['objective_opf'],row['no_objective'],rel_tol=1e-10,abs_tol=1e-6)
            (arm / 'NO_POINT_VALIDATED').write_text('Actual first NO dispatch passed original feasibility and independent OpenDSS.\n')
        if row['global_strict_criteria_met']:
            assert row['original_primal'] <= 5e-7 and row['r_max'] <= 2.5e-7 and row['s_max'] <= 1e-4
            assert row['objective_gap'] <= 1e-5
            checked['final_opendss'] = dss(arm / 'final_opendss.log')
        result['arms'][policy] = checked
    if stage == 'complete':
        starts = []
        for policy in policies:
            row = tomllib.loads((out / policy / 'assessment.toml').read_text())
            starts.append({t:digest for worker in row['worker_layout'] for t,digest in worker['initial_starts_sha256'].items()})
        if len(starts)==2:
            assert starts[0] == starts[1] and len(starts[0]) == meta['T']
            result['identical_initial_primal_starts'] = 'pass'
        else:
            prior = tomllib.loads((out/'comparison_reference'/'assessment.toml').read_text())
            prior_starts={t:digest for worker in prior['worker_layout'] for t,digest in worker['initial_starts_sha256'].items()}
            assert starts[0] == prior_starts and len(starts[0]) == meta['T']
            previous=list(csv.DictReader((out/'comparison_reference'/'trace.csv').open()))
            current=list(csv.DictReader((out/'cap30'/'trace.csv').open()))
            assert len(previous)==len(current)
            equivalent=all(math.isclose(float(a[k]),float(b[k]),rel_tol=1e-10,abs_tol=1e-10)
                for a,b in zip(previous,current) for k in ('objective','original_primal','r_max','s_max','rho'))
            result['identical_initial_primal_starts']='pass'
            result['trajectory_matches_instrumented_cap']=equivalent
        (out / 'ASSESSMENT_VALIDATED').write_text('Paired assessment checked; inspect NO/strict flags separately.\n')
(out / f'validation_{stage}.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
