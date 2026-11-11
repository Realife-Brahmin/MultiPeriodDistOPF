"""Show NO criteria and rho; a capped diagnostic has no benchmark timing claim."""
import csv
import hashlib
import json
from pathlib import Path
import sys
import tomllib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import OBJECTIVE, PRIMAL, DUAL, RHO, THRESHOLD, SYSTEM_COLORS, configure_paper_style

configure_paper_style()

baseline, candidate, output = map(Path, sys.argv[1:4])
assert (baseline / 'COMPLETE').is_file()
validation = json.loads((candidate / 'calibration_validation.json').read_text())
assert validation['stage_status_validation'] == validation['source_validation'] == 'pass'
meta = tomllib.loads((candidate / 'calibration.toml').read_text())
reference = tomllib.loads((baseline / 'summary.toml').read_text())
assert meta['reference_objective'] == reference['central_objective']
reference_hash = hashlib.sha256((baseline / 'summary.toml').read_bytes()).hexdigest()
assert meta['reference_summary_sha256'] == reference_hash
assert meta['reference_input_sha256'] == tomllib.loads((baseline / 'metadata.toml').read_text())['data_sha256']
ref = reference['central_objective']
no_tolerance = meta['no_original_tolerance']
old = list(csv.DictReader((baseline / 'trace.csv').open()))
new = list(csv.DictReader((candidate / 'calibration.csv').open()))
first = lambda rows: next((int(r['iteration']) for r in rows if
    float(r['original_primal']) <= no_tolerance and abs(float(r['objective'])-ref)/max(1,abs(ref)) <= .001), 0)
assert first(new) == meta['no_iteration']
assert not meta['no_iteration'] or validation['no_dispatch_opendss_validation'] == 'pass'
def setting(m):
    if m['rho_policy'] == 'fixed':
        return f'fixed rho={m["rho_initial"]:g}'
    return f'rho0={m["rho_initial"]:g}, balance={m.get("rho_balance_ratio", 10):g}'

series = [(old, f'C4 baseline: rho0={float(old[0]["rho"]):g}, balance=10; NO {first(old)}', '--'),
    (new, f'Redundant-bound removal: {setting(meta)}; NO {first(new) or "not reached"}', '-')]
for extra_index, extra in enumerate(sys.argv[4:]):
    adaptive = Path(extra)
    av = json.loads((adaptive / 'calibration_validation.json').read_text())
    am = tomllib.loads((adaptive / 'calibration.toml').read_text())
    assert av['stage_status_validation'] == av['source_validation'] == 'pass'
    assert am['reference_objective'] == ref
    assert am['reference_summary_sha256'] == reference_hash
    assert am['no_original_tolerance'] == no_tolerance
    ar = list(csv.DictReader((adaptive / 'calibration.csv').open()))
    assert first(ar) == am['no_iteration']
    assert not first(ar) or av['no_dispatch_opendss_validation'] == 'pass'
    series.append((ar, f'Redundant-bound removal: {setting(am)}; NO {first(ar) or "not reached"}', ('-.', ':', (0, (5, 1, 1, 1, 1, 1)))[extra_index % 3]))
fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True, layout='constrained')
for rows, label, style in series:
    x = [int(r['iteration']) for r in rows]
    gap = [max(abs(float(r['objective'])-ref)/max(1,abs(ref))/.001, 1e-12) for r in rows]
    physical = [max(float(r['original_primal'])/no_tolerance, 1e-12) for r in rows]
    axes[0].semilogy(x, gap, color=OBJECTIVE, linestyle=style, label=label)
    axes[1].semilogy(x, physical, color=PRIMAL, linestyle=style)
    axes[2].semilogy(x, [float(r['s_max'])/1e-4 for r in rows], color=DUAL, linestyle=style)
    axes[3].step(x, [float(r['rho']) for r in rows], where='post', color=RHO, linestyle=style)
    if first(rows):
        for axis in axes:
            axis.axvline(first(rows), color=THRESHOLD, linestyle=style, alpha=.35)
for axis in axes:
    axis.grid(which='major')
    axis.axvline(100, color=THRESHOLD, linestyle=':', linewidth=1)
    axis.set_xlim(1, max(110, first(old)+12))
for axis in axes[:3]:
    axis.axhline(1, color=THRESHOLD, linestyle=':', linewidth=1)
axes[0].set_ylabel('Objective gap / NO tolerance')
axes[1].set_ylabel('Original violation / NO tolerance')
axes[2].set_ylabel('Dual residual / strict threshold')
axes[3].set_ylabel(r'Penalty $\rho$')
axes[3].set_xlabel('ADMM iteration; NO requires the first two panels at or below 1')
axes[0].legend(fontsize=9)
system_label = 'med2522' if meta['system'] == 'ieee2522C_1ph' else 'large10k'
fig.suptitle(f'{system_label}, $T={meta["T"]}$, MA57, free terminal SOC\n100-step assessment; redundant-bound removal; iteration comparison only',
             fontsize=12, color=SYSTEM_COLORS[meta['system']])
fig.savefig(output)
fig.savefig(output.with_suffix('.pdf'))
plt.close(fig)
