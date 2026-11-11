"""Period distributions and measured orchestration remainder; not MPI timing."""
import csv
import hashlib
import json
from pathlib import Path
import statistics as stats
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import configure_paper_style,TADMM,SYSTEM_COLORS

base=Path(__file__).resolve().parent/'results/july_weight_comparison'
cases=[('ieee2522C_1ph','med2522',base/'med2522/results/ieee2522C_1ph_T24_july_pair_31650588/exact'),
    ('large10kC_1ph','large10k',base/'large10k/results/large10kC_1ph_T24_july_pair_31650589/exact')]
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
configure_paper_style()
fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained',sharex='col')
summary={};hashes=[]
for col,(system,name,case) in enumerate(cases):
    rows=list(csv.DictReader((case/'trace_subproblems.csv').open()))
    trace=list(csv.DictReader((case/'trace.csv').open()))
    post=[r for r in rows if int(r['iteration'])>1]
    per=[]
    for t in range(1,25):
        part=[r for r in post if int(r['period'])==t]
        per.append(dict(period=t,median_ipopt_iterations=stats.median(int(r['ipopt_iterations']) for r in part),
            max_ipopt_iterations=max(int(r['ipopt_iterations']) for r in part),
            median_native_s=stats.median(float(r['solver_s']) for r in part),
            max_native_s=max(float(r['solver_s']) for r in part)))
    groups=[[r for r in rows if int(r['iteration'])==k] for k in range(1,len(trace)+1)]
    native=sum(max(float(r['solver_s']) for r in group) for group in groups)
    stage=sum(max(float(r['stage_elapsed_s']) for r in group) for group in groups)
    core=float(trace[-1]['algorithm_s'])
    summary[name]=dict(periods=per,outer_iterations=len(trace),workers=24,
        native_critical_sum_s=native,local_stage_critical_sum_s=stage,core_s=core,
        nonlocal_core_remainder_s=core-stage,nonlocal_core_remainder_fraction=(core-stage)/core,
        native_balance_efficiency=float(trace[-1]['solver_sum_s'])/(24*native),
        distribution_excludes_first_outer_iteration=True)
    for ax,key,label in [(axes[0,col],'ipopt_iterations','Local Ipopt iterations'),(axes[1,col],'solver_s','Native Ipopt time (s)')]:
        data=[[float(r[key]) for r in post if int(r['period'])==t] for t in range(1,25)]
        ax.boxplot(data,showfliers=True,patch_artist=True,widths=.6,
            boxprops=dict(facecolor=TADMM,alpha=.45,color=TADMM),medianprops=dict(color=TADMM,linewidth=1.5),
            whiskerprops=dict(color=TADMM),capprops=dict(color=TADMM),
            flierprops=dict(marker='.',markersize=3,markerfacecolor=TADMM,markeredgecolor=TADMM,alpha=.3))
        ax.set_ylabel(label);ax.grid(axis='y',alpha=.25)
        ax.set_xticks([1,6,12,18,24],[1,6,12,18,24])
    axes[0,col].set_title(f'{name}, T=24',color=SYSTEM_COLORS[system])
    axes[1,col].set_xlabel('Owned time period')
    for f in ('trace.csv','trace_subproblems.csv'):
        hashes.append(f'{hashlib.sha256((case/f).read_bytes()).hexdigest()}  {case.parent.name}/exact/{f}\n')
fig.suptitle('MA57 local-solve distributions across 199 post-first rounds\nJuly weights; boxes: middle 50%, line: median; points: outliers')
for ext in ('png','pdf'):fig.savefig(out/f'period-distributions.{ext}',dpi=220)
(out/'summary.json').write_text(json.dumps(dict(cases=summary,pure_communication_isolated=False,inter_node_benchmark=False),indent=2)+'\n')
for p in (Path(__file__),Path(__file__).with_name('plot_style.py')):
    hashes.append(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n')
(out/'analysis-source.sha256').write_text(''.join(hashes))
