"""Validated sweep report: all rho candidates, quality, tuning and memory."""
import csv
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
import sys
import tomllib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import configure_paper_style,TADMM,THRESHOLD,SYSTEM_COLORS

def stamp(s):return datetime.fromisoformat(s.replace('Z','+00:00')).timestamp()
def readcsv(p):return list(csv.DictReader(p.open()))
def central_candidates(case,m):
    original={}
    for line in (case/'central_timed.log').read_text(errors='replace').splitlines():
        match=re.match(r'^\s*(\d+)r?\s+([-+\d.]+e[-+]\d+)\s+([-+\d.]+e[-+]\d+)\s',line)
        if match:original[int(match[1])]=float(match[3])
    trace=readcsv(case/'central_trace.csv');result={}
    for gap in (.005,.001):
        hit=next((r for r in trace if abs(float(r['objective'])-m['central_objective'])/abs(m['central_objective'])<=gap
            and original.get(int(r['iteration']),float('inf'))<=1e-5),None)
        result[str(gap)]=dict(iteration=int(hit['iteration']),seconds=float(hit['elapsed_s'])) if hit else None
    return dict(original_printed_precision_candidates=result,independent_early_dispatch_validation=False)

def memory_phases(case):
    samples=readcsv(case/'process_samples.csv');events=readcsv(case/'phase_events.csv');result={}
    for start in events:
        if not start['event'].endswith('_measured_begin'):continue
        prefix=start['event'][:-len('_measured_begin')]
        end=next(e for e in events if e['event']==prefix+'_measured_end')
        active=[r for r in samples if stamp(start['utc'])<=stamp(r['utc'])<=stamp(end['utc'])]
        if active:result[prefix]=dict(peak_pss_gib=max(int(r['pss_bytes']) for r in active)/2**30,
            peak_summed_rss_gib=max(int(r['working_set_bytes']) for r in active)/2**30,
            max_process_count=max(int(r['process_count']) for r in active))
    return result

out=Path(sys.argv[1]);cases=[Path(p) for p in sys.argv[2:]]
out.mkdir(parents=True,exist_ok=True);configure_paper_style()
fig,axes=plt.subplots(4,len(cases),figsize=(max(7,4.3*len(cases)),10),layout='constrained',squeeze=False,sharex='col')
summary=[];flat=[];hashes=[]
for col,case in enumerate(cases):
    assert (case/'SWEEP_VALIDATED').is_file()
    m=tomllib.loads((case/'sweep.toml').read_text());v=json.loads((case/'validation_complete.json').read_text())
    result=dict(T=m['T'],input_sha256=m['input_sha256'],C_B=m['C_B'],anchor_provenance=m['anchor_provenance'],
        central_full_wall_s=m['central_wall_s'],central_native_solver_s=m['central_solver_s'],
        central_posthoc=central_candidates(case,m),memory=memory_phases(case),arms=[],
        campaign_total_wall_s=m['campaign_total_wall_s'],central_build_s=m['central_build_s'],
        central_warmup_wall_s=m['central_warmup_wall_s'])
    for j,rho in enumerate(m['rhos']):
        name=f'rho{int(rho)}';arm=case/name
        a=tomllib.loads((arm/'assessment.toml').read_text())
        marks=tomllib.loads((arm/'milestones.toml').read_text()) if (arm/'milestones.toml').exists() else {}
        trace=readcsv(arm/'trace.csv');x=[int(r['iteration']) for r in trace]
        a.pop('worker_layout');a['milestones']=marks;a['rho_per_horizon']=rho/m['T']
        a['all_milestone_opendss']=v['arms'][name]['opendss'];result['arms'].append(a)
        label=f'$\\rho_0={rho:g}$';style=['-','--'][j]
        data=[ [max(1e-14,abs(float(r['objective'])-m['central_objective'])/abs(m['central_objective'])) for r in trace],
            [float(r['r_rms']) for r in trace],[float(r['owned_s_rms']) for r in trace],
            [max(1e-14,float(r['original_primal'])) for r in trace]]
        for ax,y in zip(axes[:,col],data):ax.semilogy(x,y,color=TADMM,ls=style,label=label)
        for label,row in marks.items():flat.append(dict(T=m['T'],rho_initial=rho,label=label,**row))
        for f in ('assessment.toml','milestones.toml','trace.csv','opendss_status.toml'):
            if (arm/f).exists():hashes.append(f'{hashlib.sha256((arm/f).read_bytes()).hexdigest()}  {case.name}/{name}/{f}\n')
    for i,(title,threshold) in enumerate([('Relative cost gap',.001),('Every-copy primal RMS',1e-5),
            ('Unique dual RMS',.01),('Maximum original-model violation',1e-5)]):
        ax=axes[i,col];ax.axhline(threshold,color=THRESHOLD,ls=':',label=f'Target {threshold:g}')
        ax.grid(alpha=.25)
        if col==0:ax.set_ylabel(title)
    axes[0,col].set_title(f"large10k, T={m['T']}",color=SYSTEM_COLORS['large10kC_1ph'])
    axes[0,col].legend(fontsize=8);axes[-1,col].set_xlabel('ADMM iteration')
    summary.append(result)
    for f in ('sweep.toml','validation_complete.json','process_samples.csv','phase_events.csv','central_trace.csv','central_timed.log'):
        hashes.append(f'{hashlib.sha256((case/f).read_bytes()).hexdigest()}  {case.name}/{f}\n')
fig.suptitle('Full-warm MA57 rho sweep: convergence and original-model quality\nJuly battery coefficients; free terminal SOC; 24 workers',fontsize=14)
for ext in ('png','pdf'):fig.savefig(out/f'horizon-convergence.{ext}',dpi=220)
(out/'summary.json').write_text(json.dumps(dict(single_measured_trajectory_per_configuration=True,
    benchmark_repeat_medians=False,cases=summary),indent=2)+'\n')
if flat:
    with (out/'milestones.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
for p in (Path(__file__),Path(__file__).with_name('plot_style.py')):
    hashes.append(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n')
(out/'analysis-source.sha256').write_text(''.join(hashes))
