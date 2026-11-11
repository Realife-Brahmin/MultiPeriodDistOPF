"""Raw local-solve distributions by ADMM iteration, without smoothing."""
import collections
import csv
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).parent
pilot=root/'problems/mpopf/tadmm/review_pilot'
sys.path.insert(0,str(pilot))
from plot_style import configure_july_convergence_style,OBJECTIVE,PRIMAL,DUAL,RHO
out=pilot/'results/historical_anchors_nov09'
cases=[('MA57',pilot/'results/targeted_nov09/results/run_T48_primal_dual_31657407'),
       ('MUMPS',out/'tadmm_review_nov09_mumps_fixed/results/run_T48_primal_dual_31660026')]
configure_july_convergence_style()
fig,axes=plt.subplots(2,2,figsize=(13,8.5))
fig.subplots_adjust(left=.08,right=.98,bottom=.10,top=.90,hspace=.37,wspace=.26)
all_data=[]
for col,(solver,case) in enumerate(cases):
    with (case/'rho120000/trace_subproblems.csv').open() as f:
        raw=list(csv.DictReader(f))
    grouped=collections.defaultdict(list)
    for r in raw:
        if int(r['iteration'])>1:grouped[int(r['iteration'])].append(r)
    ks=sorted(grouped);values=[[float(r['solver_s']) for r in grouped[k]] for k in ks]
    for label,color,y in [('Median',OBJECTIVE,[np.median(v) for v in values]),
                          ('95th percentile',PRIMAL,[np.quantile(v,.95,method='higher') for v in values]),
                          ('Maximum',DUAL,[max(v) for v in values])]:
        axes[0,col].plot(ks,y,label=label,color=color,lw=1.5)
    axes[0,col].set(title=f'{solver}: local solve times across 48 periods',xlabel='ADMM iteration',ylabel='Local Ipopt solve time [s]')
    axes[0,col].legend();axes[0,col].grid(alpha=.25)
    all_data.append((solver,raw,grouped))
for solver,raw,grouped in all_data:
    x=np.sort([float(r['solver_s']) for r in raw if int(r['iteration'])>1])
    axes[1,0].plot(x,np.arange(1,len(x)+1)/len(x),label=solver,color=OBJECTIVE if solver=='MA57' else RHO)
axes[1,0].set_xscale('log');axes[1,0].set(title='Distribution across all local solves',xlabel='Local Ipopt solve time [s, log scale]',ylabel='Fraction of local solves at or below time')
axes[1,0].legend();axes[1,0].grid(alpha=.25)
periods=[14,15,16,44,48]
for i,(solver,raw,grouped) in enumerate(all_data):
    counts=collections.Counter(int(max(v,key=lambda r:float(r['solver_s']))['period']) for v in grouped.values())
    axes[1,1].bar(np.arange(len(periods))+(i-.5)*.36,[100*counts[p]/len(grouped) for p in periods],width=.36,label=solver,color=OBJECTIVE if solver=='MA57' else RHO)
axes[1,1].set_xticks(range(len(periods)),periods)
axes[1,1].set(title='Recurring slowest periods (selected)',xlabel='Time-step subproblem index',ylabel='ADMM iterations where period was slowest [%]')
axes[1,1].legend();axes[1,1].grid(axis='y',alpha=.25)
fig.suptitle('large10k, T=48: local-solve stragglers',fontsize=16)
fig.text(.5,.025,'Measured exact-solve traces; iteration 1 excluded; same rho schedule. Solver time excludes warm-start transfer and coordinator work. No smoothing.',ha='center',fontsize=9)
for ext in ('png','pdf'):fig.savefig(out/f'T48-straggler-distributions.{ext}',dpi=160)
plt.close(fig)
