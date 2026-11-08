"""Compare validated historical-stop assessments and diagnostic timing."""
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys
import tomllib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import configure_paper_style,TADMM,THRESHOLD,SYSTEM_COLORS

out=Path(sys.argv[1]); cases=[Path(p) for p in sys.argv[2:]]
out.mkdir(parents=True,exist_ok=True)
configure_paper_style()
fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
variants=[]; provenance=[]; labels=[]; histories=[]
for i,case in enumerate(cases):
    assert (case/'ASSESSMENT_VALIDATED').is_file()
    m=tomllib.loads((case/'historical.toml').read_text())
    v=json.loads((case/'historical_validation.json').read_text())
    trace=list(csv.DictReader((case/'trace.csv').open()))
    histories.append((case,m,trace))
    sub=list(csv.DictReader((case/'trace_subproblems.csv').open()))
    inner=sorted(int(r['ipopt_iterations']) for r in sub if int(r['iteration'])>1)
    if variants:
        assert all(m[k]==variants[0][k] for k in ('input_sha256','C_B','T','linear_solver','stopping_rule','reference_objective'))
    x=[int(r['iteration']) for r in trace]
    mode=m.get('diagnostics_mode','cached')
    tau=m.get('rho_tau_decrease',5)
    policy='fixed' if m['rho_policy']=='fixed' else f'phase /{tau:g}'
    if m.get('rho_watchdog',False):policy+=' + watchdog'
    label=f"rho0={m['rho_initial']:g}, {policy}, {mode} checks"
    labels.append(f"rho0={m['rho_initial']:g}\n{policy}, {mode}")
    style=['-', '--', '-.', ':'][i%4]
    for ax,y in [(axes[0,0],[float(r['owned_r_rms']) for r in trace]),
        (axes[0,1],[float(r['original_primal']) for r in trace]),
        (axes[1,0],[max(1e-14,abs(float(r['objective'])-m['reference_objective'])/abs(m['reference_objective'])) for r in trace])]:
        ax.semilogy(x,y,ls=style,color=TADMM,label=label)
    row={k:m[k] for k in ('rho_initial','iterations','historical_stopped','original_primal','objective','reference_objective',
        'algorithm_s','loop_wall_s','C_B','input_sha256','T','linear_solver','stopping_rule')}
    row.update(rho_policy=m['rho_policy'],rho_tau_decrease=tau,rho_watchdog=m.get('rho_watchdog',False),diagnostics_mode=mode,objective_relative_gap=abs(m['objective']-m['reference_objective'])/abs(m['reference_objective']),
        post_first_inner_median=statistics.median(inner),post_first_inner_p95=inner[int(.95*(len(inner)-1))],post_first_inner_max=max(inner),
        physical_NO_iteration=m['physical_NO_iteration'],dispatch_s=m.get('dispatch_s'),physical_check_s=m.get('physical_check_s'),
        numerical_prefix_matches_prior_exact=v['numerical_prefix_matches_prior_exact'],max_voltage_difference_pu=v['max_voltage_difference_pu'])
    variants.append(row)
    for f in ('historical.toml','historical_validation.json','trace.csv','trace_subproblems.csv'):
        provenance.append(f"{hashlib.sha256((case/f).read_bytes()).hexdigest()}  {case.name}/{f}\n")
for case,m,trace in histories:
    if m.get('rho_watchdog',False):
        ref=next(rows for _,other,rows in histories if other['rho_policy']==m['rho_policy'] and
            other['rho_initial']==m['rho_initial'] and other.get('rho_tau_decrease',5)==m.get('rho_tau_decrease',5) and not other.get('rho_watchdog',False))
        k=int(re.search(r'LEGACY_WATCHDOG k=(\d+)',(case/'solve.log').read_text())[1])
        assert len(ref)>=k and len(trace)>=k
        assert all(math.isclose(float(x[key]),float(y[key]),abs_tol=1e-10,rel_tol=1e-10)
            for x,y in zip(ref[:k],trace[:k]) for key in ('objective','original_primal','r_max','s_max','owned_r_rms','owned_s_rms','rho'))
        (out/'watchdog-prefix-audit.json').write_text(json.dumps(dict(matches_phase5_through_iteration=k,watchdog_fires_at=k),indent=2)+'\n')
        provenance.append(f'{hashlib.sha256((case/"solve.log").read_bytes()).hexdigest()}  {case.name}/solve.log\n')
for ax,threshold,title in [(axes[0,0],1e-4,'Owned-period primal RMS'),(axes[0,1],1e-5,'Maximum original-model violation'),
        (axes[1,0],.001,'Relative objective gap')]:
    ax.axhline(threshold,color=THRESHOLD,ls=':',label=f'Threshold {threshold:g}')
    ax.set(xlabel='ADMM iteration',title=title);ax.grid(alpha=.25)
axes[0,0].legend(fontsize=8)
ax=axes[1,1]
bars=ax.bar(range(len(variants)),[v['algorithm_s'] for v in variants],color=TADMM,alpha=.7)
for bar,hatch,row in zip(bars,['','//','xx','..'],variants):
    bar.set_hatch(hatch)
    ax.annotate(f"{row['algorithm_s']:.1f} s\n{row['iterations']} iters",(bar.get_x()+bar.get_width()/2,bar.get_height()),
        xytext=(0,4),textcoords='offset points',ha='center',fontsize=9)
ax.set_xticks(range(len(labels)),labels,rotation=15,ha='right',fontsize=8)
ax.set(ylabel='Core solve time, excluding diagnostics (s)',title='Single exploratory observations')
ax.set_ylim(0,max(v['algorithm_s'] for v in variants)*1.25)
fig.suptitle('large10k, T=24: rho and diagnostic overhead\nJuly weights, historical RMS stopping, MA57',color=SYSTEM_COLORS['large10kC_1ph'])
for ext in ('png','pdf'):fig.savefig(out/f'rho-diagnostics.{ext}',dpi=220)
(out/'summary.json').write_text(json.dumps(dict(benchmark_timing_claim=False,variants=variants),indent=2)+'\n')
for p in (Path(__file__),Path(__file__).with_name('plot_style.py')):
    provenance.append(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n')
(out/'analysis-source.sha256').write_text(''.join(provenance))
