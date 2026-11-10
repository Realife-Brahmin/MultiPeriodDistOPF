"""Chronological lossless Storage discrepancies, using verified raw tables."""
import csv
import gzip
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).parent;pilot=root/'problems/mpopf/tadmm/review_pilot'
sys.path.insert(0,str(pilot))
from plot_style import configure_july_convergence_style,OBJECTIVE,PRIMAL,DUAL
base=pilot/'results/sequential_dss_nov09'
configure_july_convergence_style()
fig,axes=plt.subplots(1,2,figsize=(12,5.3))
fig.subplots_adjust(left=.085,right=.985,bottom=.20,top=.81,wspace=.27)
series={}
for name,label,color in [('central-model','Centralized',OBJECTIVE),('cap30-model','Cap 30',PRIMAL),('cap10-model','Cap 10',DUAL)]:
    with gzip.open(base/'tables'/name/'battery.csv.gz','rt') as f:rows=list(csv.DictReader(f))
    curve=[max(abs(float(r['engine_energy_balance_error_kwh'])) for r in rows if int(r['t'])==t) for t in range(1,25)]
    axes[0].semilogy(range(1,25),np.maximum(curve,1e-10),color=color,label=label,lw=1.7,marker='o',markersize=3)
    if name=='cap10-model':series['cap10']=rows
for field,label,color in [('requested_bound_error_kwh','Requested powers integrated without clipping',OBJECTIVE),('engine_bound_error_kwh','OpenDSS stored energy after reserve clamping',PRIMAL)]:
    curve=[max(float(r[field]) for r in series['cap10'] if int(r['t'])==t) for t in range(1,25)]
    axes[1].semilogy(range(1,25),np.maximum(curve,1e-10),color=color,label=label,lw=1.7,marker='s',markersize=3)
axes[0].set_title('Energy accounting: maximum across batteries')
axes[0].set_ylabel(r'Maximum $|E_{t}-E_{t-1}+\Delta t\,P_{\rm actual}|$ [kWh]')
axes[1].set_title('Cap 10: operating SOC-bound violations')
axes[1].set_ylabel('Maximum battery energy-bound violation [kWh]')
for ax in axes:
    ax.axhline(.01,color='red',ls='--',lw=1,label='Energy tolerance: 0.01 kWh')
    ax.set_xlabel('End of dispatch interval [h]');ax.set_xlim(.5,24.5);ax.set_xticks([1,4,8,12,16,20,24]);ax.grid(alpha=.25)
    ax.legend(fontsize=8,loc='lower left')
fig.suptitle('Chronological OpenDSS replay: large10k, T=24',fontsize=15,y=.98)
fig.text(.5,.90,'100% charging/discharging efficiency; zero idling loss; live Storage; energy never reset between intervals',ha='center',fontsize=10)
fig.text(.5,.045,'Validation only, outside solve timing. A bounded OpenDSS SOC does not certify energy conservation when the native engine clips at reserve.',ha='center',fontsize=9)
for ext in ('png','pdf'):fig.savefig(base/f'chronological-discrepancies.{ext}',dpi=160)
plt.close(fig)
