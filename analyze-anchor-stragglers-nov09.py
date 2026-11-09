"""Passive local-solve and critical-path analysis; clipping is hypothetical."""
import collections
import csv
import json
import math
from pathlib import Path
import statistics
import tomllib

root=Path(__file__).parent
base=root/'problems/mpopf/tadmm/review_pilot/results/historical_anchors_nov09'
cases=sorted(base.glob('tadmm_review_*/results/run_*'))
cases.append(root/'problems/mpopf/tadmm/review_pilot/results/targeted_nov09/results/run_T48_primal_dual_31657407')
def rows(p):
    with p.open() as f:return list(csv.DictReader(f))
def stats(x):
    x=sorted(x)
    return dict(median=statistics.median(x),p95=x[math.ceil(.95*len(x))-1],maximum=max(x))
result=[]
for case in cases:
    meta=tomllib.loads((case/'sweep.toml').read_text())
    for rho in meta['rhos']:
        arm=case/f'rho{int(rho)}'
        raw=rows(arm/'trace_subproblems.csv'); trace=rows(arm/'trace.csv')
        grouped=collections.defaultdict(list)
        for r in raw:
            if int(r['iteration'])>1:grouped[int(r['iteration'])].append(r)
        phases={}
        for phase,rounds in [('all_after_first',list(grouped)),('first40_after_first',[k for k in grouped if k<=40]),('after40',[k for k in grouped if k>40])]:
            if not rounds:continue
            selected=[r for k in rounds for r in grouped[k]]
            maxima=[];ratios=[];top=collections.Counter();period_rows=collections.defaultdict(list);worst=[]
            clip_s={2:0.,3:0.};critical_solver=critical_stage=median_solver_sum=0.
            for k in rounds:
                group=grouped[k]
                med=statistics.median(float(r['solver_s']) for r in group)
                high=max(group,key=lambda r:float(r['solver_s']))
                mx=float(high['solver_s']);maxima.append(mx);ratios.append(mx/med)
                top[int(high['period'])]+=1
                critical_solver+=mx
                critical_stage+=max(float(r['stage_elapsed_s']) for r in group)
                median_solver_sum+=med
                for multiple in clip_s:clip_s[multiple]+=max(0,mx-multiple*med)
                worst.append(dict(iteration=k,period=int(high['period']),solver_s=mx,
                    round_median_s=med,ratio=mx/med,ipopt_iterations=int(high['ipopt_iterations']),rho=float(trace[k-1]['rho'])))
                for r in group:period_rows[int(r['period'])].append(float(r['solver_s']))
            phases[phase]=dict(rounds=len(rounds),local_solver_s=stats([float(r['solver_s']) for r in selected]),
                local_ipopt_iterations=stats([int(r['ipopt_iterations']) for r in selected]),round_max_solver_s=stats(maxima),
                round_max_to_median_ratio=stats(ratios),rounds_over_2x=sum(x>2 for x in ratios),rounds_over_5x=sum(x>5 for x in ratios),
                top_slowest_periods=[dict(period=p,rounds=n,share=n/len(rounds),median_s=statistics.median(period_rows[p]),mean_s=statistics.mean(period_rows[p])) for p,n in top.most_common(6)],
                worst_rounds=sorted(worst,key=lambda r:r['solver_s'],reverse=True)[:8],
                critical_solver_sum_s=critical_solver,critical_stage_sum_s=critical_stage,round_median_solver_sum_s=median_solver_sum,
                idealized_critical_solver_seconds_removed_by_2x_or_3x_median=clip_s)
        result.append(dict(case=case.name,system=meta['system'],T=meta['T'],solver=meta['linear_solver'],phases=phases))
payload=dict(note='One period per worker; iteration1 excluded. Local times are solver-only. Round ratios compare max to median across time-step solves. Ideal clipping assumes no overhead or convergence/quality change and is not a measured speedup. Phase after40 still excludes first iteration. Raw arrays remain unchanged.',cases=result)
(base/'straggler-distributions.json').write_text(json.dumps(payload,indent=2)+'\n')
for r in result:
    p=r['phases']['all_after_first']
    print(r['system'],r['T'],r['solver'],json.dumps({k:p[k] for k in ('local_solver_s','local_ipopt_iterations','round_max_solver_s','round_max_to_median_ratio','rounds_over_2x','rounds_over_5x','top_slowest_periods','worst_rounds','idealized_critical_solver_seconds_removed_by_2x_or_3x_median')}))
