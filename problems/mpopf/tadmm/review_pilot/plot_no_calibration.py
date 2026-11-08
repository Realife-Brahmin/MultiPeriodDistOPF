"""Show NO criteria and rho; a capped diagnostic has no benchmark timing claim."""
import csv
import json
from pathlib import Path
import sys
import tomllib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import OBJECTIVE, PRIMAL, RHO, THRESHOLD

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
fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True, layout='constrained')
for rows, label, style in ((old, f'Corrected C4 baseline: NO {first(old)}; initial rho=48', '--'),
    (new, f'DDP screening, rho={meta["rho_initial"]:g}: NO {first(new) or "not reached"}', '-')):
    x = [int(r['iteration']) for r in rows]
    gap = [max(abs(float(r['objective'])-ref)/max(1,abs(ref))/.001, 1e-12) for r in rows]
    physical = [max(float(r['original_primal'])/1e-5, 1e-12) for r in rows]
    axes[0].semilogy(x, gap, color=OBJECTIVE, linestyle=style, label=label)
    axes[1].semilogy(x, physical, color=PRIMAL, linestyle=style)
    axes[2].step(x, [float(r['rho']) for r in rows], where='post', color=RHO, linestyle=style)
    if first(rows):
        for axis in axes:
            axis.axvline(first(rows), color=THRESHOLD, linestyle=style, alpha=.35)
for axis in axes:
    axis.grid(alpha=.2)
    axis.axvline(100, color=THRESHOLD, linestyle=':', linewidth=1)
    axis.set_xlim(1, 180)
for axis in axes[:2]:
    axis.axhline(1, color=THRESHOLD, linestyle=':', linewidth=1)
axes[0].set_ylabel('Objective gap / NO tolerance')
axes[1].set_ylabel('Original violation / NO tolerance')
axes[2].set_ylabel('Penalty rho')
axes[2].set_xlabel('ADMM iteration; NO requires both upper panels at or below 1')
axes[0].legend(frameon=False, fontsize=9)
fig.suptitle('large10k T=24, MA57, free terminal SOC\nOne capped calibration; exact screening and initial rho changed; no timing comparison', fontsize=11)
fig.savefig(output, dpi=170)
plt.close(fig)
