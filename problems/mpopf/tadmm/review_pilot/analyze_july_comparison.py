"""Reproduce convergence and straggler diagnostics from a validated paired case."""
import csv
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
try:
    import tomllib
except ImportError:
    import tomli as tomllib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from plot_style import TADMM, THRESHOLD, SYSTEM_COLORS, configure_paper_style

configure_paper_style()
case, out = map(Path, sys.argv[1:3])
assert (case / 'ASSESSMENT_VALIDATED').is_file()
meta = tomllib.loads((case / 'comparison.toml').read_text())
validation = json.loads((case / 'validation_complete.json').read_text())
assert validation['identical_initial_primal_starts'] == validation['source_validation'] == 'pass'
out.mkdir(parents=True, exist_ok=True)
ref = meta['reference_objective']
name = {'ieee2522C_1ph': 'med2522', 'large10kC_1ph': 'large10k'}[meta['system']]
datasets = {}
summaries = {}
rounds = []
entries=[('exact',case/'exact',validation['arms']['exact']),('cap30',case/'cap30',validation['arms']['cap30'])]
light_validation=None
if len(sys.argv)>3:
    light=Path(sys.argv[3])
    assert (light/'ASSESSMENT_VALIDATED').is_file()
    lm=tomllib.loads((light/'comparison.toml').read_text())
    assert lm['input_sha256']==meta['input_sha256'] and lm['C_B']==meta['C_B']
    light_validation=json.loads((light/'validation_complete.json').read_text())
    entries.append(('cap30_light',light/'cap30',light_validation['arms']['cap30']))
styles={'exact':'-','cap30':'--','cap30_light':'-.'}
labels={'exact':'Exact local solves','cap30':'30-step cap + full local reports','cap30_light':'30-step cap; global checks'}
for policy,arm,arm_validation in entries:
    info = tomllib.loads((arm / 'assessment.toml').read_text())
    trace = list(csv.DictReader((arm / 'trace.csv').open()))
    stages = list(csv.DictReader((arm / 'trace_subproblems.csv').open()))
    groups = {k:[] for k in range(1,len(trace)+1)}
    for s in stages:
        groups[int(s['iteration'])].append(s)
    stats = []
    bottlenecks = Counter()
    for k, rows in groups.items():
        worst = max(rows, key=lambda s:float(s['solver_s']))
        stage_worst = max(rows, key=lambda s:float(s['stage_elapsed_s']))
        times = [float(s['solver_s']) for s in rows]
        median = statistics.median(times)
        row = dict(policy=policy, iteration=k, median_solver_s=median,
            max_solver_s=max(times), max_to_median=max(times)/max(median,1e-12),
            bottleneck_period=int(worst['period']), bottleneck_ipopt_iterations=int(worst['ipopt_iterations']),
            maximum_ipopt_iterations=max(int(s['ipopt_iterations']) for s in rows),
            max_stage_elapsed_s=float(stage_worst['stage_elapsed_s']),
            framework_at_slowest_stage_s=float(stage_worst['stage_elapsed_s'])-float(stage_worst['solver_s']),
            cutoffs=sum(s['iteration_limit_hit']=='true' for s in rows))
        stats.append(row); rounds.append(row)
        if k>1:
            bottlenecks[row['bottleneck_period']] += 1
    tail = stats[1:] or stats
    limited = [s for s in stages if s['iteration_limit_hit']=='true']
    no = info['no_iteration']
    gap = lambda r:abs(float(r['objective'])-ref)/max(1,abs(ref))
    qualifies=lambda r:gap(r)<=.001 and float(r['original_primal'])<=meta['no_original_tolerance']
    bad=[n for n,r in enumerate(trace) if not qualifies(r)]
    sustained=(bad[-1]+1) if bad else 0
    summaries[policy] = dict(iterations=info['iterations'], no_iteration=no,
        target_under_100=info['target_under_100'], global_strict_criteria_met=info['global_strict_criteria_met'],
        no_qualifying_iterations=sum(qualifies(r) for r in trace),
        sustained_no_from_iteration=(int(trace[sustained]['iteration']) if sustained<len(trace) else None),
        objective_gap_at_end=info['objective_gap'], original_primal_at_end=info['original_primal'],
        no_loop_clock_s=float(trace[no-1]['loop_wall_s']) if no else None,
        no_algorithm_s=float(trace[no-1]['algorithm_s']) if no else None,
        total_loop_s=info['loop_wall_s'], total_algorithm_s=info['algorithm_s'],
        truncated_local_solves=len(limited),
        truncated_local_solves_after_first_iteration=sum(int(s['iteration'])>1 for s in limited),
        largest_truncated_local_primal=max((float(s['local_primal']) for s in limited if math.isfinite(float(s['local_primal']))),default=None),
        max_ipopt_iterations=max(int(s['ipopt_iterations']) for s in stages),
        median_round_max_solver_s=statistics.median(r['max_solver_s'] for r in tail),
        median_round_max_stage_elapsed_s=statistics.median(r['max_stage_elapsed_s'] for r in tail),
        median_framework_at_slowest_stage_s=statistics.median(r['framework_at_slowest_stage_s'] for r in tail),
        largest_round_max_solver_s=max(r['max_solver_s'] for r in tail),
        median_straggler_ratio=statistics.median(r['max_to_median'] for r in tail),
        largest_straggler_ratio=max(r['max_to_median'] for r in tail),
        bottleneck_period_counts=dict(bottlenecks.most_common()),
        cost_only_0p5pct_iteration=next((int(r['iteration']) for r in trace if gap(r)<=.005),None),
        cost_only_0p1pct_iteration=next((int(r['iteration']) for r in trace if gap(r)<=.001),None),
        stage_status_validation=arm_validation['stage_status_validation'],
        no_opendss_validated=('no_opendss' in arm_validation))
    datasets[policy] = (info, trace, stats)
report = dict(system=meta['system'], T=meta['T'], C_B=meta['C_B'],
    input_sha256=meta['input_sha256'], rho_initial=meta['rho_initial'], rho_policy=meta['rho_policy'],
    benchmark_timing_claim=False, timing_note=meta['timing_note'],
    cutoff_diagnostics_note='The instrumented cap arm evaluates full local primal feasibility after each cutoff; its stage elapsed time includes that additional diagnostic work. Raw solver times exclude it.',
    straggler_statistics_exclude_first_iteration=True,
    cost_only_note='Objective-only crossings do not assert feasibility and are not the historical full July criterion.', arms=summaries)
if light_validation is not None:
    report['light_trajectory_matches_full_reports']=light_validation['trajectory_matches_instrumented_cap']
(out / 'summary.json').write_text(json.dumps(report,indent=2)+'\n')
with (out / 'straggler_rounds.csv').open('w',newline='') as stream:
    writer=csv.DictWriter(stream,fieldnames=rounds[0]);writer.writeheader();writer.writerows(rounds)

fig, axes = plt.subplots(5,1,figsize=(10,11),sharex=True,layout='constrained')
for policy,_,_ in entries:
    style=styles[policy]
    info, trace, stats = datasets[policy]
    label = f'{labels[policy]}; NO {info["no_iteration"] or "not reached"}'
    x = [int(r['iteration']) for r in trace]
    values = [[max(abs(float(r['objective'])-ref)/max(1,abs(ref))/.001,1e-12) for r in trace],
        [max(float(r['original_primal'])/meta['no_original_tolerance'],1e-12) for r in trace],
        [max(float(r['r_max'])/2.5e-7,1e-12) for r in trace],
        [max(float(r['s_max'])/1e-4,1e-12) for r in trace]]
    for axis,y in zip(axes,values):
        axis.semilogy(x,y,color=TADMM,linestyle=style,label=label)
    axes[4].step(x,[float(r['rho']) for r in trace],where='post',color=TADMM,linestyle=style)
    if info['no_iteration']:
        for axis in axes:
            axis.axvline(info['no_iteration'],color=THRESHOLD,linestyle=style,alpha=.35)
for axis,label in zip(axes,['Objective gap / NO tolerance','Original violation / NO tolerance',
    'All-copy consensus / strict threshold','Dual residual / strict threshold',r'Penalty $\rho$']):
    axis.set_ylabel(label);axis.grid(which='major');axis.axvline(100,color=THRESHOLD,linestyle=':',linewidth=.8)
for axis in axes[:4]:
    axis.axhline(1,color=THRESHOLD,linestyle=':',linewidth=.8)
axes[0].legend()
axes[-1].set_xlabel('ADMM iteration; NO requires the first two panels at or below 1')
fig.suptitle(f'{name}, $T={meta["T"]}$, July $C_B={meta["C_B"]:.3g}$, MA57\nCorrected losses; current periodic inputs; one short-warmup assessment per arm',
    fontsize=12,color=SYSTEM_COLORS[meta['system']])
for ext in ('png','pdf'):
    fig.savefig(out / f'{name}-july-convergence.{ext}')
plt.close(fig)

fig, axes = plt.subplots(3,1,figsize=(10,8),sharex=True,layout='constrained')
for n,(policy,_,_) in enumerate(entries):
    style=styles[policy];medianstyle=((0,(1,2)),(0,(2,2)),(0,(1,3,3,3)))[n]
    info,trace,stats=datasets[policy]
    x=[r['iteration'] for r in stats]
    label=labels[policy]
    axes[0].semilogy(x,[r['max_solver_s'] for r in stats],color=TADMM,linestyle=style,label=label+': slowest')
    axes[0].semilogy(x,[r['median_solver_s'] for r in stats],color=TADMM,linestyle=medianstyle,linewidth=1,alpha=.65,label=label+': median')
    axes[1].plot(x,[r['max_to_median'] for r in stats],color=TADMM,linestyle=style)
    axes[2].step(x,[r['cutoffs'] for r in stats],where='post',color=TADMM,linestyle=style)
axes[0].set_ylabel('Local Ipopt solver time (s)');axes[0].legend(ncol=2)
axes[1].set_ylabel('Slowest / median local time');axes[1].axhline(1,color=THRESHOLD,linestyle=':',linewidth=.8)
axes[2].set_ylabel('Local iteration-limit hits');axes[2].yaxis.set_major_locator(MaxNLocator(integer=True))
axes[2].set_xlabel('ADMM iteration; local solver times exclude framework/diagnostic overhead')
for axis in axes:
    axis.grid(which='major')
fig.suptitle(f'{name}, $T={meta["T"]}$: local-solve outliers and explicit cutoffs\nMatched July weight, model, initial starts and rho settings',fontsize=12,color=SYSTEM_COLORS[meta['system']])
for ext in ('png','pdf'):
    fig.savefig(out / f'{name}-july-stragglers.{ext}')
plt.close(fig)
source=Path(__file__).resolve()
(out / 'analysis-source.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n'
    for p in (source,source.with_name('plot_style.py'))))
print(json.dumps(report,indent=2))
