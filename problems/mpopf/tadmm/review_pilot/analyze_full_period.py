"""Matched implementation comparison, retaining every supplied completed arm."""
import csv
from datetime import datetime
import hashlib
import json
from pathlib import Path
import statistics
import sys
import tomllib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import configure_paper_style,TADMM,THRESHOLD,SYSTEM_COLORS

def rows(path):
    with path.open() as f:return list(csv.DictReader(f))

def stamp(value):return datetime.fromisoformat(value.replace('Z','+00:00')).timestamp()

def memory(case):
    files=sorted(case.glob('process_samples.*.csv'))
    if not files and (case/'process_samples.csv').exists():files=[case/'process_samples.csv']
    events=rows(case/'phase_events.csv');result={}
    for start in events:
        if not start['event'].endswith('_measured_begin'):continue
        phase=start['event'].removesuffix('_measured_begin')
        end=next((e for e in events if e['event']==phase+'_measured_end'),None)
        if end is None:continue
        nodes={}
        for file in files:
            active=[r for r in rows(file) if stamp(start['utc'])<=stamp(r['utc'])<=stamp(end['utc'])]
            if active:nodes[file.name]=dict(peak_pss_gib=max(int(r['pss_bytes']) for r in active)/2**30,
                                          maximum_processes=max(int(r['process_count']) for r in active))
        result[phase]=dict(nodes=nodes,sum_of_node_peak_pss_gib_upper_bound=sum(n['peak_pss_gib'] for n in nodes.values()),
                           note='Sum of independent node peaks is an upper bound, not a synchronized aggregate peak.')
    return result

def critical(details,layout,key):
    by_iteration={}
    for row in details:by_iteration.setdefault(int(row['iteration']),{})[int(row['period'])]=float(row.get(key,0))
    return sum(max(sum(periods[t] for t in w['periods']) for w in layout) for periods in by_iteration.values())

out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);cases=[Path(p) for p in sys.argv[2:]]
configure_paper_style()
horizons=sorted({tomllib.loads((p/'sweep.toml').read_text())['T'] for p in cases})
fig,axes=plt.subplots(3,len(horizons),figsize=(5*len(horizons),8),layout='constrained',squeeze=False,sharex='col')
summary=[];hashes=[];references={};curves={}
for i,case in enumerate(cases):
    assert (case/'SWEEP_VALIDATED').is_file()
    m=tomllib.loads((case/'sweep.toml').read_text());validation=json.loads((case/'validation_complete.json').read_text())
    identity=(m['T'],m['input_sha256'],m['C_B'],m['central_initial_starts_sha256'])
    references.setdefault(m['T'],identity);assert references[m['T']]==identity
    for rho in m['rhos']:
        arm=case/f'rho{int(rho)}';a=tomllib.loads((arm/'assessment.toml').read_text())
        trace=rows(arm/'trace.csv');details=rows(arm/'trace_subproblems.csv')
        label=f"T={m['T']}, {m['workers']} workers, rho={rho:g}, {m.get('warm_start','primal')}, {m.get('objective_update','legacy')}"
        layout=a['worker_layout'];native=critical(details,layout,'solver_s');stage=critical(details,layout,'stage_elapsed_s')
        inner=[int(r['ipopt_iterations']) for r in details]
        a.pop('worker_layout')
        a.update(label=label,case=case.name,T=m['T'],rho=rho,workers=m['workers'],
                 warm_start=m.get('warm_start','primal'),objective_update=m.get('objective_update','legacy'),
                 native_group_critical_s=native,stage_group_critical_s=stage,
                 core_outside_stage_critical_s=a['core_s']-stage,
                 warm_transfer_group_critical_s=critical(details,layout,'warm_start_s'),
                 objective_update_group_critical_s=critical(details,layout,'objective_update_s'),
                 inner_iterations_median=statistics.median(inner),inner_iterations_max=max(inner),
                 inner_iterations_p95=sorted(inner)[int(.95*(len(inner)-1))],
                 central_full_wall_s=m['central_wall_s'],memory=memory(case),
                 validation=validation['arms'][arm.name])
        summary.append(a)
        x=[float(r['algorithm_s']) for r in trace]
        series=[[max(1e-14,abs(float(r['objective'])-m['central_objective'])/abs(m['central_objective'])) for r in trace],
                [float(r['r_rms']) for r in trace],[max(1e-14,float(r['original_primal'])) for r in trace]]
        col=horizons.index(m['T']);curve=curves.get(m['T'],0);curves[m['T']]=curve+1
        short=f"{m.get('warm_start','primal').replace('_','–')} start, $\\rho_0={rho:g}$"
        for ax,y in zip(axes[:,col],series):ax.semilogy(x,y,color=TADMM,ls=['-','--',':','-.'][curve%4],label=short)
        for f in (arm/'assessment.toml',arm/'trace.csv',arm/'trace_subproblems.csv',case/'sweep.toml',case/'validation_complete.json'):
            hashes.append(f'{hashlib.sha256(f.read_bytes()).hexdigest()}  {case.name}/{f.relative_to(case)}\n')
for col,T in enumerate(horizons):
    for i,(title,tol) in enumerate(zip(['Relative cost gap','Every-copy RMS','Original-model maximum violation'],[.001,1e-5,1e-5])):
        ax=axes[i,col];ax.axhline(tol,color=THRESHOLD,ls=':');ax.grid(alpha=.25)
        if col==0:ax.set_ylabel(title)
    axes[0,col].set_title(f'large10k, T={T}, {T} workers',color=SYSTEM_COLORS['large10kC_1ph'])
    axes[0,col].legend(fontsize=8);axes[-1,col].set_xlabel('Core solve time (s)')
solvers=sorted({tomllib.loads((p/'sweep.toml').read_text())['linear_solver'].upper() for p in cases})
fig.suptitle('/'.join(solvers)+': full period parallelism and exact local warm starts',fontsize=14)
for ext in ('png','pdf'):fig.savefig(out/f'implementation-comparison.{ext}',dpi=220)
(out/'summary.json').write_text(json.dumps(dict(single_measured_trajectory=True,arms=summary),indent=2)+'\n')
hashes.append(f'{hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}  {Path(__file__).name}\n')
(out/'analysis-source.sha256').write_text(''.join(hashes))
