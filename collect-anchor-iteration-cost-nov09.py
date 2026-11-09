"""Observed iteration costs and hypothetical count-only break-even."""
import csv
import json
import math
from pathlib import Path
import statistics
import tomllib

base=Path(__file__).parent/'problems/mpopf/tadmm/review_pilot/results/historical_anchors_nov09'
def rows(p):
    with p.open() as f:return list(csv.DictReader(f))
def dist(x):
    x=sorted(x)
    return dict(median=statistics.median(x),p95=x[math.ceil(.95*len(x))-1],maximum=max(x))
result=[]
for c in sorted(base.glob('tadmm_review_*/results/run_*')):
    m=tomllib.loads((c/'sweep.toml').read_text())
    for rho in m['rhos']:
        arm=c/f'rho{int(rho)}'
        trace=rows(arm/'trace.csv'); local=rows(arm/'trace_subproblems.csv')
        core=[float(r['algorithm_s']) for r in trace]
        delta=[x-y for x,y in zip(core,[0]+core[:-1])]
        warmed=[r for r in local if int(r['iteration'])>1]
        maxima={}
        for r in local:
            k=int(r['iteration'])
            maxima[k]=max(maxima.get(k,0),float(r['stage_elapsed_s']))
        avg=core[-1]/len(trace)
        result.append(dict(system=m['system'],T=m['T'],solver=m['linear_solver'],case=c.name,
            iterations=len(trace),central_s=m['central_wall_s'],mean_core_round_s=avg,
            core_round_s=dist(delta),warmed_local_solver_s=dist([float(r['solver_s']) for r in warmed]),
            warmed_local_stage_s=dist([float(r['stage_elapsed_s']) for r in warmed]),
            warmed_local_ipopt_iterations=dist([int(r['ipopt_iterations']) for r in warmed]),
            warmed_round_max_stage_s=dist([s for k,s in maxima.items() if k>1]),
            prefix40_core_s=core[39] if len(core)>=40 else None,
            count_only_40_round_projection_s=40*avg,
            count_only_strictly_faster_max_iterations=math.ceil(m['central_wall_s']/avg)-1,
            observed_last_prefix_below_central=max((k for k,s in enumerate(core,1) if s<m['central_wall_s']),default=0),
            note='Projection holds measured average iteration cost fixed. Observed prefixes are unfinished states, not converged solutions. Existing core timer excludes build/warmup/diagnostics.'))
(base/'iteration-costs.json').write_text(json.dumps(result,indent=2)+'\n')
for r in result:
    print(r['system'],r['T'],r['solver'],json.dumps({k:r[k] for k in ('mean_core_round_s','warmed_local_solver_s','warmed_round_max_stage_s','prefix40_core_s','central_s','count_only_strictly_faster_max_iterations','observed_last_prefix_below_central')}))
