"""Restore the July objective/r/s/rho four-panel view using untouched run data.

Style source: 50b8901c324ba9eefd0d672c6e484e0ededc7a57,
envs/tadmm/Plotter.jl::plot_tadmm_ldf_convergence (July 31, 2026).
July r/s mean owned-period/unique RMS, not original-model feasibility.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import tomllib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
from plot_style import OBJECTIVE,PRIMAL,DUAL,RHO,configure_july_convergence_style

STYLE_COMMIT='50b8901c324ba9eefd0d672c6e484e0ededc7a57'
STYLE_PATH='envs/tadmm/Plotter.jl'

def power_label(value):
    exponent=int(np.floor(np.log10(value)));mantissa=value/10**exponent
    return (rf'10^{{{exponent}}}' if np.isclose(mantissa,1)
            else rf'{mantissa:g}\times10^{{{exponent}}}')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',type=Path)
    parser.add_argument('output',type=Path,help='Output stem; PNG, PDF, CSV and JSON are saved')
    parser.add_argument('--rho',type=float,help='Select an arm when the case contains multiple starts')
    args=parser.parse_args()
    m=tomllib.loads((args.case/'sweep.toml').read_text())
    if args.rho is None:
        if len(m['rhos'])!=1:parser.error('--rho is required for a multi-start case')
        rho=m['rhos'][0]
    else:rho=args.rho
    assert rho in m['rhos']
    arm=args.case/f'rho{int(rho)}';src=arm/'trace.csv'
    with src.open() as f:rows=list(csv.DictReader(f))
    assessment=tomllib.loads((arm/'assessment.toml').read_text())
    assert len(rows)==assessment['iterations']
    x=np.array([int(r['iteration']) for r in rows])
    assert np.array_equal(x,np.arange(1,len(rows)+1))
    columns={'objective':'objective','r_norm':'owned_r_rms','s_norm':'owned_s_rms',
             'all_copy_r_norm':'r_rms','rho':'rho','core_s':'algorithm_s'}
    data={key:np.array([float(r[col]) for r in rows]) for key,col in columns.items()}
    assert all(np.isfinite(v).all() for v in data.values())
    assert all((data[key]>=0).all() for key in ('r_norm','s_norm','all_copy_r_norm'))
    assert (data['rho']>0).all()
    configure_july_convergence_style()
    fig,axes=plt.subplots(4,1,figsize=(9,13),sharex=True)
    fig.subplots_adjust(left=.14,right=.975,bottom=.065,top=.925,hspace=.30)
    n=len(rows);markevery=max(1,n//22)
    common=dict(markevery=markevery,markersize=3.8,markeredgecolor='black',markeredgewidth=.35)
    axes[0].plot(x,data['objective'],color=OBJECTIVE,marker='o',label='tADMM objective',**common)
    central=m['central_objective']
    axes[0].axhline(central,color='#CD6600',lw=1.7,ls='--',label=f'Centralized reference: ${central:,.2f}')
    axes[0].plot([],[],color='none',label=f"tADMM final: ${data['objective'][-1]:,.2f}")
    formatter=ScalarFormatter(useMathText=True,useOffset=False)
    formatter.set_powerlimits((0,0));axes[0].yaxis.set_major_formatter(formatter)
    floor=1e-16
    axes[1].semilogy(x,np.maximum(data['r_norm'],floor),color=PRIMAL,marker='s',label='Owned-period RMS',**common)
    axes[1].semilogy(x,np.maximum(data['all_copy_r_norm'],floor),color=PRIMAL,ls=':',alpha=.65,lw=1.3,label='All-copy RMS')
    primal=m['final_owned_rms'];assert primal==m['final_all_copy_rms']
    axes[1].axhline(primal,color='red',ls='--',lw=1.2,alpha=.7,label=rf'Run target: $\epsilon_r={power_label(primal)}$')
    axes[1].axhline(1e-3,color='#808080',ls=':',lw=1,label=r"$10^{-3}$ reference (not this run's stop)")
    axes[2].semilogy(x,np.maximum(data['s_norm'],floor),color=DUAL,marker='D',label='Unique dual RMS',**common)
    axes[2].axhline(m['final_dual_rms'],color='red',ls='--',lw=1.2,alpha=.7,
                    label=rf"Run target: $\epsilon_s={power_label(m['final_dual_rms'])}$")
    axes[3].step(x,data['rho'],where='post',color=RHO,label=r'Penalty $\rho$')
    axes[3].plot(x[::markevery],data['rho'][::markevery],ls='none',color=RHO,marker='h',
                 markersize=4.3,markeredgecolor='black',markeredgewidth=.35)
    axes[3].set_yscale('log')
    titles=['Objective','Primal Residual','Dual Residual','Penalty Schedule']
    labels=['Objective function [$]',r'Primal residual $\|r\|_{\mathrm{RMS}}$',
            r'Dual residual $\|s\|_{\mathrm{RMS}}$',r'Penalty parameter $\rho$']
    step=max(1,n//5);ticks=[1]+list(range(step,n+1,step))
    if n-ticks[-1]<.3*step and len(ticks)>1:ticks[-1]=n
    else:ticks.append(n)
    for ax,title,label in zip(axes,titles,labels):
        ax.set_title(title,pad=6);ax.set_ylabel(label)
        ax.grid(True,which='major');ax.set_xlim(.5,n+.5)
        ax.legend(loc='upper right',framealpha=.95)
        ax.tick_params(labelbottom=True)
        ax.set_xticks(sorted(set(ticks)))
        ax.set_xlabel(r'Iteration $(k)$',labelpad=3)
    system={'large10kC_1ph':'large10k','ieee2522C_1ph':'med2522','ieee123C_1ph':'ieee123'}.get(m['system'],m['system'])
    policy=m.get('rho_policy','legacy_phase')
    method=(rf"high start $\to$ fixed $\rho={m['rho_floor']:g}$" if policy=='bootstrap_fixed' else 'July adaptive penalty')
    fig.suptitle('tADMM Convergence Summary',fontsize=16,y=.987)
    fig.text(.5,.956,rf"{system}, $T={m['T']}$ | {m['linear_solver'].upper()} | {method} | $\rho_0={rho:g}$",
             ha='center',fontsize=11)
    if m.get('local_solve_policy')=='early_cap30':
        boundary=m['cap_end_iteration']+1
        fig.text(.5,.936,f"Local cap {m['cap_iterations']} in rounds 2 to {boundary-1} ({m.get('cap_screen','feasible')} screen); exact from {boundary}",ha='center',fontsize=9)
        for ax in axes:ax.axvline(boundary,color='#707070',ls=':',lw=1)
    end_status=('RMS/dual target reached' if assessment['final_rms_met'] else 'Iteration budget reached; RMS/dual target not met')
    fig.text(.5,.018,f"{end_status} | Core solve: {assessment['core_s']:.1f} s",ha='center',fontsize=10,color='#505050')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    for ext in ('png','pdf'):fig.savefig(args.output.with_suffix('.'+ext),dpi=160)
    plt.close(fig)
    csvout=args.output.with_suffix('.csv')
    with csvout.open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['iteration',*columns])
        writer.writerows([int(k),*(float(data[key][i]) for key in columns)] for i,k in enumerate(x))
    # Exact numeric round trip: export/plot reads raw columns with no smoothing,
    # residual renormalization, sample dropping or changed convergence flags.
    with csvout.open() as f:roundtrip=list(csv.DictReader(f))
    assert all(float(r[key])==data[key][i] for i,r in enumerate(roundtrip) for key in columns)
    manifest=dict(case=str(args.case),rho_initial=rho,iterations=n,style_source=dict(commit=STYLE_COMMIT,path=STYLE_PATH),
        colors=dict(objective=OBJECTIVE,primal=PRIMAL,dual=DUAL,rho=RHO),font='Computer Modern (cmr10; mathtext cm)',
        source_columns=columns,all_copy_curve='Dotted green; historical owned RMS is solid green',
        local_solve_policy=m.get('local_solve_policy','exact'),
        log_display_floor_for_exact_zeros=floor,
        zero_counts={key:int((data[key]==0).sum()) for key in ('r_norm','s_norm','all_copy_r_norm')},
        quality=dict(final_rms_met=assessment['final_rms_met'],strict_criteria_met=assessment['strict_criteria_met'],
                     original_primal=assessment['original_primal'],objective_gap=assessment['objective_gap']),
        sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                [src,args.case/'sweep.toml',arm/'assessment.toml',Path(__file__),Path(__file__).with_name('plot_style.py'),csvout]})
    args.output.with_suffix('.json').write_text(json.dumps(manifest,indent=2)+'\n')

if __name__=='__main__':main()
