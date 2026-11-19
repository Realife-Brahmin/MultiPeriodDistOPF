"""Actual battery-power-integrated SOC versus the saved consensus trajectory."""
import hashlib
from pathlib import Path
import sys
import tomllib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import configure_paper_style,TADMM,THRESHOLD,SYSTEM_COLORS

src=Path(sys.argv[1]);out=Path(sys.argv[2]);m=tomllib.loads(src.read_text())
configure_paper_style()
fig,axes=plt.subplots(2,1,figsize=(8,6),layout='constrained',sharex=True)
rows=[(key,row) for key,row in m.items() if isinstance(row,dict) and 'integrated_SOC_percent' in row]
endpoints=[(key,row) for key,row in rows if key.endswith('/endpoint.jls')]
if endpoints:rows=endpoints
for i,(key,row) in enumerate(rows):
    rho=key.split('/')[0].replace('rho','');style=['-','--',':','-.'][i%4]
    x=row['hours'];soc=row['integrated_SOC_percent'];consensus=row['consensus_SOC_percent']
    axes[0].plot(x,soc,color=TADMM,ls=style,label=f'Bus {row["max_SOC_error_bus"]}, $\\rho_0={rho}$')
    axes[0].plot(x,consensus,color=TADMM,ls=style,alpha=.25,linewidth=1)
    lo=row['minimum_SOC_percent'];hi=row['maximum_SOC_percent']
    axes[1].plot(x,[max(lo-s,s-hi,0) for s in soc],color=TADMM,ls=style,label=f'$\\rho_0={rho}$')
axes[0].axhline(lo,color=THRESHOLD,ls=':',label=f'Minimum SOC {lo:g}%')
axes[0].axhline(hi,color=THRESHOLD,ls='--',label=f'Maximum SOC {hi:g}%')
axes[0].set_ylabel('Battery SOC (%)')
axes[1].set_ylabel('SOC bound deviation\n(percentage points)')
axes[1].set_xlabel('Hour')
for ax in axes:ax.grid(alpha=.25);ax.legend(fontsize=9)
axes[0].set_title(f"large10k, T={m['T']}: worst SOC-bound trajectories\nLossless model; endpoints; faint curves: consensus",color=SYSTEM_COLORS['large10kC_1ph'])
for ext in ('png','pdf'):fig.savefig(out.with_suffix('.'+ext),dpi=220)
out.with_suffix('.sha256').write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n'
    for p in [src,Path(__file__),Path(__file__).with_name('plot_style.py')]))
