"""Independent first-crossing, full-copy, warm-start and source audit."""
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

out=Path(sys.argv[1]); stage=sys.argv[2]
m=tomllib.loads((out/'sweep.toml').read_text())
assert (m['system']=='large10kC_1ph' and m['T'] in (24,48,96)) or (m['system']=='ieee2522C_1ph' and m['T']==96)
expected_workers=m['T'] if m.get('execution_policy')=='one_per_period' else 24
assert m['workers']==expected_workers and m['compute_cpu_budget']==expected_workers+1 and m['linear_solver'] in ('ma57','mumps')
actual_solver=re.search(r'running with linear solver ([\w-]+)',(out/'central_timed.log').read_text(),re.I)
assert actual_solver and actual_solver[1].casefold()==m['linear_solver'].casefold()
assert m['terminal_soc_penalty']==0 and m['central_original_primal']<=1e-6
for line in (out/'source.sha256').read_text().splitlines():
    digest,rel=line.split('  ',1)
    assert hashlib.sha256((out/'source'/rel).read_bytes()).hexdigest()==digest

def flags(r):
    owned=float(r['owned_r_rms']); allr=float(r['r_rms']); dual=float(r['owned_s_rms'])
    gap=abs(float(r['objective'])-m['central_objective'])/abs(m['central_objective'])
    result={}
    for key,eps in [('1e-3',1e-3),('1e-4',1e-4),('1e-5',1e-5)]:
        result['owned_raw_'+key]=owned<=eps
        result['all_raw_'+key]=allr<=eps
        result['owned_dual_'+key]=owned<=eps and dual<=.01
        result['all_dual_'+key]=owned<=eps and allr<=eps and dual<=.01
        result['paper_cost_'+key]=owned<=eps and gap<=.005
        result['current_cost_'+key]=owned<=eps and allr<=eps and gap<=.001
    result['physical_NO']=float(r['original_primal'])<=1e-5 and gap<=.001
    result['final_rms']=owned<=1e-5 and allr<=1e-5 and dual<=.01
    if m.get('local_solve_policy')=='early_cap30' and int(r['iteration'])<=m['cap_end_iteration']:
        result['final_rms']=False
    return result

def dss(path):
    text=path.read_text(errors='replace')
    dv=re.search(r'voltage: largest difference ([\deE+.-]+) pu',text)
    return dict(scope='network_snapshot_only',converged='converged=true' in text and 'converged=false' not in text,
        max_voltage_difference_pu=float(dv[1]) if dv else None,
        zero_voltage_violations=bool(re.search(r'below limit: 0 of .*above limit: 0,',text)))

central_dss=dss(out/'central_opendss.log')
assert central_dss['converged'] and central_dss['zero_voltage_violations'] and central_dss['max_voltage_difference_pu']<=1e-4
names=[f'rho{int(r)}' for r in m['rhos']] if stage=='complete' else [stage]
result=dict(central_opendss=central_dss,arms={},final_benchmark_claim=False)
central_sequence=None
if m.get('sequential_validation_version')==1:
    central_sequence=tomllib.loads((out/m['central_sequential_summary']).read_text())
    assert central_sequence['sequential_schedule_pass']==m['central_sequential_schedule_pass']
    assert central_sequence['efficiency_mode']=='model' and central_sequence['T']==m['T']
result['central_sequential']=central_sequence
starts=[]
for name in names:
    arm=out/name
    a=tomllib.loads((arm/'assessment.toml').read_text())
    marks=tomllib.loads((arm/'milestones.toml').read_text()) if (arm/'milestones.toml').exists() else {}
    trace=list(csv.DictReader((arm/'trace.csv').open()))
    warm=list(csv.DictReader((arm/'warmup/trace.csv').open()))
    assert len(trace)==a['iterations']<=m['maximum_iterations'] and len(warm)==a['warmup_iterations']
    assert len(warm)==len(trace)
    assert all(math.isclose(float(x[k]),float(y[k]),abs_tol=1e-10,rel_tol=1e-10)
        for x,y in zip(warm,trace) for k in ('objective','original_primal','owned_r_rms','r_rms','owned_s_rms','rho'))
    first={};rho=a['rho_initial'];seen=False;high=0
    for k,r in enumerate(trace,1):
        assert int(r['iteration'])==k and float(r['rho'])==rho
        for key,yes in flags(r).items():
            if yes:first.setdefault(key,k)
        rr=float(r['owned_r_rms']);ss=float(r['owned_s_rms'])
        if m.get('rho_policy','legacy_phase')=='bootstrap_fixed':
            if k%5==0 and float(r['r_rms'])<=m.get('bootstrap_gate',math.inf):
                rho=max(m['rho_floor'],rho/m.get('bootstrap_decrease',5))
        else:
            if rr>2e-4:high+=1
            elif rr<.5e-4:high=0
            if high>=15:rho=min(1e6,2.5*rho);high=0
            if k%5==0:
                seen=seen or rr<=1e-4
                if not seen:
                    if rr>1e-4:rho=min(1e6,2*rho)
                elif rr>5*ss:rho=min(1e6,2*rho)
                elif ss>5*rr and rr<=1e-4:rho=max(1,rho/5)
    assert set(first)==set(marks)
    for key,k in first.items():
        assert marks[key]['iteration']==k
        assert len(marks[key]['snapshot_sha256'])==64
        row=trace[k-1]
        for src,dst in [('algorithm_s','core_s'),('owned_r_rms','owned_r_rms'),('r_rms','all_r_rms'),('objective','objective'),('original_primal','original_primal')]:
            assert math.isclose(float(row[src]),marks[key][dst],abs_tol=1e-10,rel_tol=1e-10)
    assert a['final_rms_met']==('final_rms' in first)
    assert first.get('final_rms',len(trace))==len(trace)
    assert a['final_rms_met'] or len(trace)==m['maximum_iterations']
    layout=a['worker_layout'];assert len(layout)==expected_workers
    assert len({(w['hostname'],w['affinity']) for w in layout})==expected_workers
    assert sorted(t for w in layout for t in w['periods'])==list(range(1,m['T']+1))
    assert all(len(w['periods'])==m['T']//expected_workers and w['julia_threads']==w['blas_threads']==1 and w['openmp_threads']=='1' for w in layout)
    if m.get('execution_policy')=='one_per_period':
        assert len({w['hostname'] for w in layout})==m['nodes']
        assert all(w['affinity']==w['expected_pin'] and w['warm_start']==m['warm_start'] and w['objective_update']==m['objective_update'] for w in layout)
    starts.append({t:h for w in layout for t,h in w['initial_starts_sha256'].items()})
    for folder in [arm,arm/'warmup']:
        rows=list(csv.DictReader((folder/'trace.csv').open()));sub=list(csv.DictReader((folder/'trace_subproblems.csv').open()))
        assert len(sub)==len(rows)*m['T']
        assert {(int(r['iteration']),int(r['period'])) for r in sub}=={(k,t) for k in range(1,len(rows)+1) for t in range(1,m['T']+1)}
        if m.get('local_solve_policy')=='persistent_cap':
            assert m['cap_screen']=='finite'
            for r in sub:
                assert int(r['ipopt_iterations'])<=m['cap_iterations']
                assert r['cap_recovery']=='false' and r['cold_retry']=='false'
                assert r['local_primal_evaluated']=='false'
                limited=r['status']=='ITERATION_LIMIT'
                assert r['status'] in ('LOCALLY_SOLVED','OPTIMAL','ITERATION_LIMIT')
                assert (r['iteration_limit_hit']=='true')==limited
                assert (r['inexact_accepted']=='true')==limited
                assert int(r['cap_attempt_iterations'])==(int(r['ipopt_iterations']) if limited else 0)
        elif m.get('local_solve_policy')=='early_cap30':
            assert len(rows)>m['cap_end_iteration'] or not a['final_rms_met']
            for r in sub:
                k=int(r['iteration']);accepted=r['inexact_accepted']=='true';recovery=r['cap_recovery']=='true'
                if accepted:
                    assert 1<k<=m['cap_end_iteration'] and not recovery
                    assert r['status']=='ITERATION_LIMIT' and r['iteration_limit_hit']=='true'
                    assert r['local_primal_evaluated']=='true'
                    assert math.isfinite(float(r['local_primal']))
                    if m.get('cap_screen','feasible')=='feasible':
                        assert float(r['local_primal'])<=m['cap_primal_tolerance']
                    assert int(r['cap_attempt_iterations'])==m['cap_iterations']
                else:
                    assert r['status'] in ('LOCALLY_SOLVED','OPTIMAL') and r['iteration_limit_hit']=='false'
                if recovery:
                    assert 1<k<=m['cap_end_iteration'] and int(r['ipopt_iterations'])>=int(r['cap_attempt_iterations'])
                if k>m['cap_end_iteration'] or k==1:
                    assert not accepted and not recovery and int(r['cap_attempt_iterations'])==0
        else:
            assert all(r['status'] in ('LOCALLY_SOLVED','OPTIMAL') and r['iteration_limit_hit']=='false' for r in sub)
        if m.get('local_solve_policy')=='warm_retry':
            for r in sub:
                if r['cold_retry']=='true':
                    assert int(r['warm_attempt_iterations'])>=m['warm_retry_limit']
                    assert int(r['ipopt_iterations'])>=int(r['warm_attempt_iterations'])
    checks={p:dss(arm/(Path(p).stem+'_opendss.log')) for p in set(v['snapshot'] for v in marks.values())}
    checks['endpoint.jls']=dss(arm/'endpoint_opendss.log')
    if 'physical_NO' in marks:
        q=checks[marks['physical_NO']['snapshot']]
        assert q['converged'] and q['zero_voltage_violations'] and q['max_voltage_difference_pu']<=1e-4
    result['arms'][name]=dict(first_crossings=first,final_rms_met=a['final_rms_met'],core_s=a['core_s'],
        loop_wall_s=a['loop_wall_s'],original_primal=a['original_primal'],objective_gap=a['objective_gap'],opendss=checks)
    seq=None
    if m.get('sequential_validation_version')==1:
        seq=tomllib.loads((out/a['sequential_summary']).read_text())
        assert seq['sequential_schedule_pass']==a['sequential_schedule_pass']
        assert seq['efficiency_mode']=='model' and seq['T']==m['T']
        if (arm/'endpoint.jls').exists():assert hashlib.sha256((arm/'endpoint.jls').read_bytes()).hexdigest()==seq['solution_sha256']
    result['arms'][name]['sequential']=seq
    result['arms'][name]['full_mpopf_validation_pass']=(bool(central_sequence['sequential_schedule_pass'] and seq['sequential_schedule_pass'] and a['original_primal']<=1e-5 and a['objective_gap']<=.001) if seq is not None else None)
if len(starts)>1:assert all(s==starts[0] for s in starts)
(out/f'validation_{stage}.json').write_text(json.dumps(result,indent=2)+'\n')
if stage=='complete':(out/'SWEEP_VALIDATED').write_text('Collection/source/timing audit complete. Snapshot and chronological OpenDSS flags and original-model feasibility are separate; this marker alone does not certify full MPOPF feasibility.\n')
print(json.dumps(result,indent=2))
