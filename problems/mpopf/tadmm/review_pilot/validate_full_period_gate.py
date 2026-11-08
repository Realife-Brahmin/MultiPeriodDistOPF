"""Independent audit of the saved full-period execution gate."""
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys
import tomllib

folder=Path(sys.argv[1]);out=Path(sys.argv[2])
assert (folder/'FULL_PERIOD_GATE_VALIDATED').is_file()
for line in (folder/'source.sha256').read_text().splitlines():
    digest,rel=line.split('  ',1)
    assert hashlib.sha256((folder/'source'/rel).read_bytes()).hexdigest()==digest
def rows(path):
    with path.open() as f:return list(csv.DictReader(f))
summary={};trajectories={}
for mode,warm,update in [('legacy','primal','legacy'),('batched','primal','batched'),('dual','primal_dual','batched')]:
    meta=tomllib.loads((folder/f'{mode}.toml').read_text());layout=meta['workers']
    assert len(layout)==96 and len({(w['hostname'],w['affinity']) for w in layout})==96
    assert sorted(t for w in layout for t in w['periods'])==list(range(1,97))
    hosts={w['hostname'] for w in layout};assert len(hosts)==2
    assert all(sum(w['hostname']==host for w in layout)==48 for host in hosts)
    assert all(len(w['periods'])==1 and w['affinity']==w['expected_pin'] and
               w['julia_threads']==w['blas_threads']==1 and w['openmp_threads']=='1' and
               w['warm_start']==warm and w['objective_update']==update for w in layout)
    trace=rows(folder/f'{mode}.csv');sub=rows(folder/f'{mode}_subproblems.csv')
    assert len(trace)==8 and len(sub)==8*96
    assert {(int(r['iteration']),int(r['period'])) for r in sub}=={(k,t) for k in range(1,9) for t in range(1,97)}
    assert all(r['status'] in ('LOCALLY_SOLVED','OPTIMAL') and r['iteration_limit_hit']=='false' for r in sub)
    log=(folder/f'{mode}-opendss.log').read_text(errors='replace')
    dv=float(re.search(r'voltage: largest difference ([\deE+.-]+) pu',log)[1])
    assert 'converged=true' in log and 'converged=false' not in log and dv<=1e-4
    assert re.search(r'below limit: 0 of .*above limit: 0,',log)
    steady=[r for r in sub if int(r['iteration'])>1]
    summary[mode]=dict(period_workers=96,workers_per_host={h:48 for h in hosts},local_solves=len(sub),
        median_inner_iterations=statistics.median(int(r['ipopt_iterations']) for r in steady),
        maximum_inner_iterations=max(int(r['ipopt_iterations']) for r in sub),
        **{f'median_{key}':statistics.median(float(r[key]) for r in steady)
           for key in ('solver_s','stage_elapsed_s','warm_start_s','objective_update_s')},
        max_opendss_voltage_difference_pu=dv,final_original_primal=float(trace[-1]['original_primal']),
        final_all_rms=float(trace[-1]['r_rms']),final_unique_dual_rms=float(trace[-1]['owned_s_rms']))
    trajectories[mode]=trace
keys=('objective','original_primal','owned_r_rms','r_rms','rho')
assert all(a[key]==b[key] for a,b in zip(trajectories['legacy'],trajectories['batched']) for key in keys)
assert all(math.isclose(float(a[key]),float(b[key]),abs_tol=1e-9,rel_tol=1e-10)
           for a,b in zip(trajectories['batched'],trajectories['dual']) for key in keys)
result=dict(short_gate_not_convergence_benchmark=True,iteration_one_excluded_from_local_medians=True,
            batched_numerical_path_exactly_matches=True,full_dual_path_agrees_within_audit_tolerances=True,
            source_checks_passed=True,modes=summary,
            warning='Modes run sequentially with different JIT warmup histories; their total gate times are not fair speedup measurements.')
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
