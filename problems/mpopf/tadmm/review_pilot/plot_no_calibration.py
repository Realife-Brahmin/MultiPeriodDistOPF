"""Show NO criteria and rho; a capped diagnostic has no benchmark timing claim."""
import csv
import json
from pathlib import Path
import sys
import tomllib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import OBJECTIVE, PRIMAL, DUAL, RHO, THRESHOLD

baseline, candidate, output = map(Path, sys.argv[1:4])
assert (baseline / 'COMPLETE').is_file()
validation = json.loads((candidate / 'calibration_validation.json').read_text())
assert validation['stage_status_validation'] == validation['source_validation'] == 'pass'
meta = tomllib.loads((candidate / 'calibration.toml').read_text())
reference = tomllib.loads((baseline / 'summary.toml').read_text())
assert meta['reference_objective'] == reference['central_objective']
ref = reference['central_objective']
old = list(csv.DictReader((baseline / 'trace.csv').open()))
new = list(csv.DictReader((candidate / 'calibration.csv').open()))
first = lambda rows: next((int(r['iteration']) for r in rows if
    float(r['original_primal']) <= 1e-5 and abs(float(r['objective'])-ref)/max(1,abs(ref)) <= .001), 0)
assert first(new) == meta['no_iteration']
assert not meta['no_iteration'] or validation['no_dispatch_opendss_validation'] == 'pass'
series = [(old, f'Corrected C4 baseline: NO {first(old)}; initial rho=48', '--'),
    (new, f'DDP screening, {meta["rho_policy"]} rho={meta["rho_initial"]:g}: NO {first(new) or "not reached by 100"}', '-')]
for extra_index, extra in enumerate(sys.argv[4:]):
    adaptive = Path(extra)
    av = json.loads((adaptive / 'calibration_validation.json').read_text())
    am = tomllib.loads((adaptive / 'calibration.toml').read_text())
    assert av['stage_status_validation'] == av['source_validation'] == 'pass'
    assert am['reference_objective'] == ref
    ar = list(csv.DictReader((adaptive / 'calibration.csv').open()))
    assert first(ar) == am['no_iteration']
    assert not first(ar) or av['no_dispatch_opendss_validation'] == 'pass'
    series.append((ar, f'DDP screening, {am["rho_policy"]} start={am["rho_initial"]:g}: NO {first(ar) or "not reached by 100"}', ('-.', ':')[extra_index % 2]))
fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True, layout='constrained')
for rows, label, style in series:
    x = [int(r['iteration']) for r in rows]
    gap = [max(abs(float(r['objective'])-ref)/max(1,abs(ref))/.001, 1e-12) for r in rows]
    physical = [max(float(r['original_primal'])/1e-5, 1e-12) for r in rows]
    axes[0].semilogy(x, gap, color=OBJECTIVE, linestyle=style, label=label)
    axes[1].semilogy(x, physical, color=PRIMAL, linestyle=style)
    axes[2].semilogy(x, [float(r['s_max'])/1e-4 for r in rows], color=DUAL, linestyle=style)
    axes[3].step(x, [float(r['rho']) for r in rows], where='post', color=RHO, linestyle=style)
    if first(rows):
        for axis in axes:
            axis.axvline(first(rows), color=THRESHOLD, linestyle=style, alpha=.35)
for axis in axes:
    axis.grid(alpha=.2)
    axis.axvline(100, color=THRESHOLD, linestyle=':', linewidth=1)
    axis.set_xlim(1, 180)
for axis in axes[:3]:
    axis.axhline(1, color=THRESHOLD, linestyle=':', linewidth=1)
axes[0].set_ylabel('Objective gap / NO tolerance')
axes[1].set_ylabel('Original violation / NO tolerance')
axes[2].set_ylabel('Dual residual / strict threshold')
axes[3].set_ylabel('Penalty rho')
axes[3].set_xlabel('ADMM iteration; NO requires the first two panels at or below 1')
axes[0].legend(frameon=False, fontsize=9)
fig.suptitle('large10k T=24, MA57, free terminal SOC\nCapped penalty calibration; exact DDP screening; iteration comparison only', fontsize=11)
fig.savefig(output, dpi=170)
plt.close(fig)
