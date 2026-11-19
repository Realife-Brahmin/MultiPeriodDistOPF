"""Historical stop versus separately measured original feasibility."""
import csv
import hashlib
import json
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import configure_paper_style,TADMM,THRESHOLD,SYSTEM_COLORS

case=Path(sys.argv[1]); out=Path(sys.argv[2])
assert (case/'ASSESSMENT_VALIDATED').is_file()
v=json.loads((case/'historical_validation.json').read_text())
rows=list(csv.DictReader((case/'trace.csv').open()))
configure_paper_style()
fig,axes=plt.subplots(3,1,figsize=(8,8),sharex=True,layout='constrained')
x=[int(r['iteration']) for r in rows]
for ax,key,tol,label in zip(axes,('owned_r_rms','owned_s_rms','original_primal'),(1e-4,1e-2,1e-5),
        ('Owned-period primal RMS','Unique consensus-change dual RMS','Maximum original-model violation')):
    ax.semilogy(x,[float(r[key]) for r in rows],color=TADMM,label='tADMM')
    ax.axhline(tol,color=THRESHOLD,ls=':',label=f'Threshold {tol:g}')
    if v['historical_iteration']:
        ax.axvline(v['historical_iteration'],color=TADMM,ls='--',alpha=.7)
    ax.set_ylabel(label)
    ax.grid(alpha=.25)
    ax.legend(loc='best')
axes[0].set_title(f"large10k, T=24: historical RMS stop at {v['historical_iteration']}\nJuly battery weight, fixed rho=48, exact MA57 locals",
    color=SYSTEM_COLORS['large10kC_1ph'])
axes[-1].set_xlabel('ADMM iteration; original feasibility is reported separately')
out.parent.mkdir(parents=True,exist_ok=True)
for ext in ('.png','.pdf'): fig.savefig(out.with_suffix(ext),dpi=220)
files=[Path(__file__),Path(__file__).with_name('plot_style.py'),case/'trace.csv',case/'historical_validation.json']
out.with_suffix('.sha256').write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n' for p in files))
