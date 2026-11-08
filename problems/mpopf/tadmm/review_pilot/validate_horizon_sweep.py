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
assert m['T'] in (24,48,96) and m['system']=='large10kC_1ph'
expected_workers=m['T'] if m.get('execution_policy')=='one_per_period' else 24
assert m['workers']==expected_workers and m['compute_cpu_budget']==expected_workers+1 and m['linear_solver']=='ma57'
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
    return result

def dss(path):
    text=path.read_text(errors='replace')
    dv=re.search(r'voltage: largest difference ([\deE+.-]+) pu',text)
    return dict(converged='converged=true' in text and 'converged=false' not in text,
        max_voltage_difference_pu=float(dv[1]) if dv else None,
        zero_voltage_violations=bool(re.search(r'below limit: 0 of .*above limit: 0,',text)))

central_dss=dss(out/'central_opendss.log')
assert central_dss['converged'] and central_dss['zero_voltage_violations'] and central_dss['max_voltage_difference_pu']<=1e-4
names=[f'rho{int(r)}' for r in m['rhos']] if stage=='complete' else [stage]
result=dict(central_opendss=central_dss,arms={},final_benchmark_claim=False)
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
        assert all(r['status'] in ('LOCALLY_SOLVED','OPTIMAL') and r['iteration_limit_hit']=='false' for r in sub)
    checks={p:dss(arm/(Path(p).stem+'_opendss.log')) for p in set(v['snapshot'] for v in marks.values())}
    checks['endpoint.jls']=dss(arm/'endpoint_opendss.log')
    if 'physical_NO' in marks:
        q=checks[marks['physical_NO']['snapshot']]
        assert q['converged'] and q['zero_voltage_violations'] and q['max_voltage_difference_pu']<=1e-4
    result['arms'][name]=dict(first_crossings=first,final_rms_met=a['final_rms_met'],core_s=a['core_s'],
        loop_wall_s=a['loop_wall_s'],original_primal=a['original_primal'],objective_gap=a['objective_gap'],opendss=checks)
if len(starts)>1:assert all(s==starts[0] for s in starts)
(out/f'validation_{stage}.json').write_text(json.dumps(result,indent=2)+'\n')
if stage=='complete':(out/'SWEEP_VALIDATED').write_text('Sources, full warmups, milestones and physical/OpenDSS flags audited separately.\n')
print(json.dumps(result,indent=2))
