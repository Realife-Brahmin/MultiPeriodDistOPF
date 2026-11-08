"""Check historical stopping separately from physical feasibility and DSS error."""
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

out=Path(sys.argv[1])
m=tomllib.loads((out/'historical.toml').read_text())
assert m['system']=='large10kC_1ph' and m['T']==24 and m['linear_solver']=='ma57'
assert m['stopping_rule']=='historical_owned_rms'
assert math.isclose(m['C_B'],8.013987384856764e-8,rel_tol=1e-14)
assert m['rho_policy'] in ('fixed','legacy_phase') and m['rho_initial']>0
assert m.get('rho_tuning',False) or m['rho_initial']==48
assert (m['historical_primal_tolerance'],m['historical_dual_tolerance'],m['stagnation_window'],m['stagnation_range_tolerance'])==(1e-4,1e-2,5,1e-3)
for line in (out/'source.sha256').read_text().splitlines():
    digest,rel=line.split('  ',1)
    assert hashlib.sha256((out/'source'/rel).read_bytes()).hexdigest()==digest
trace=list(csv.DictReader((out/'trace.csv').open()))
assert len(trace)==m['iterations']<=200
objectives=[]; first=0; no=0
for k,r in enumerate(trace,1):
    assert int(r['iteration'])==k
    assert all(math.isfinite(float(r[key])) for key in ('objective','original_primal','owned_r_rms','owned_s_rms'))
    objectives.append(float(r['objective']))
    recent=objectives[-5:]
    stagnant=len(recent)==5 and max(recent)>0 and (max(recent)-min(recent))/abs(max(recent))<1e-3
    met=float(r['owned_r_rms'])<=1e-4 and (float(r['owned_s_rms'])<=1e-2 or stagnant)
    assert r['historical_stagnant']==str(stagnant).lower()
    assert r['historical_stop_met']==str(met).lower()
    if met and not first: first=k
    if not no and float(r['original_primal'])<=1e-5 and abs(objectives[-1]-m['reference_objective'])/max(1,abs(m['reference_objective']))<=.001: no=k
assert m['historical_stopped']==bool(first)
assert (first==len(trace)) if first else len(trace)==200
assert m['physical_NO_iteration']==no
if no:
    assert m['no_opendss_process_success']
    no_log=(out/'no_opendss.log').read_text()
    assert 'converged=true' in no_log
    assert re.search(r'below limit: 0 of .*above limit: 0,',no_log)
    assert float(re.search(r'voltage: largest difference ([\deE+.-]+) pu',no_log)[1])<=1e-4
prefix=None
if m['rho_initial']==48 and m['rho_policy']=='fixed':
    baseline=list(csv.DictReader((out/'baseline_trace.csv').open()))
    assert len(baseline)>=len(trace)
    prefix=all(math.isclose(float(a[key]),float(b[key]),rel_tol=1e-10,abs_tol=1e-10)
        for a,b in zip(trace,baseline) for key in ('objective','original_primal','r_max','s_max','rho'))
    assert prefix
rho=m['rho_initial']; seen=False
for r in trace:
    assert float(r['rho'])==rho
    if m['rho_policy']=='legacy_phase' and int(r['iteration'])%5==0:
        rr=float(r['owned_r_rms']); ss=float(r['owned_s_rms'])
        seen=seen or rr<=1e-4
        if not seen:
            if rr>1e-4:rho=min(1e6,2*rho)
        elif rr>5*ss:rho=min(1e6,2*rho)
        elif ss>5*rr and rr<=1e-4:rho=max(1,rho/m.get('rho_tau_decrease',5))
layout=m['worker_layout']
assert len(layout)==24 and len({w['affinity'] for w in layout})==24
assert all(w['local_solve_policy']=='exact' and w['julia_threads']==w['blas_threads']==1 and w['openmp_threads']=='1' for w in layout)
starts=lambda rows:{t:h for w in rows for t,h in w['initial_starts_sha256'].items()}
assert starts(layout)==starts(tomllib.loads((out/'baseline_assessment.toml').read_text())['worker_layout'])
for key in ('objective','original_primal','owned_r_rms','owned_s_rms'):
    assert math.isclose(m[key],float(trace[-1][key]),rel_tol=1e-12,abs_tol=1e-12)
for phase in ('warmup','trace'):
    outer=list(csv.DictReader((out/f'{phase}.csv').open()))
    sub=list(csv.DictReader((out/f'{phase}_subproblems.csv').open()))
    assert len(sub)==24*len(outer)
    assert {(int(r['iteration']),int(r['period'])) for r in sub}=={(k,t) for k in range(1,len(outer)+1) for t in range(1,25)}
    assert all(r['status'] in ('LOCALLY_SOLVED','OPTIMAL') and r['iteration_limit_hit']=='false' for r in sub)
log=(out/'endpoint_opendss.log').read_text()
converged='converged=true' in log
match=re.search(r'voltage: largest difference ([\deE+.-]+) pu',log)
dv=float(match[1]) if match else None
assert converged and dv is not None # Successful load flow is distinct from an accepted OPF point.
physical=(m['original_primal']<=1e-5 and abs(m['objective']-m['reference_objective'])/max(1,abs(m['reference_objective']))<=.001)
result=dict(historical_iteration=first,physical_NO_iteration=no,endpoint_physical_NO=physical,
    numerical_prefix_matches_prior_exact=prefix,initial_starts_match=True,rho=m['rho_initial'],
    original_primal=m['original_primal'],opendss_converged=converged,max_voltage_difference_pu=dv,
    algorithm_s=m['algorithm_s'],loop_wall_s=m['loop_wall_s'],benchmark_timing_claim=False,
    physical_NO_opendss_validated=bool(no),rho_policy=m['rho_policy'])
(out/'historical_validation.json').write_text(json.dumps(result,indent=2)+'\n')
(out/'ASSESSMENT_VALIDATED').write_text('Stopping/source/status checks passed; physical and OpenDSS metrics are separate.\n')
print(json.dumps(result,indent=2))
