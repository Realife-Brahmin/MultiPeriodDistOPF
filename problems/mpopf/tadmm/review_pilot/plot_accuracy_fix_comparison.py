"""Compare the preserved failed warm-up with a fully validated corrected run."""
import csv
import json
from pathlib import Path
import sys
import tomllib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import PRIMAL, DUAL, RHO, THRESHOLD

old, new, out = map(lambda x: Path(x).resolve(), sys.argv[1:4])
assert (new / 'COMPLETE').exists(), 'Do not promote an incomplete corrected run'
validation = json.loads((new / 'validation.json').read_text())
assert all(validation.get(key) == 'pass' for key in (
    'optimization_validation', 'warmup_trace_validation', 'opendss_validation', 'source_snapshot_validation'))
old_meta = tomllib.loads((old / 'metadata.toml').read_text())
new_meta = tomllib.loads((new / 'metadata.toml').read_text())
for key in ('data_sha256', 'C_B', 'terminal_soc_penalty', 'rho_policy', 'rho_initial', 'linear_solver'):
    assert old_meta[key] == new_meta[key], (key, old_meta[key], new_meta[key])
summary = tomllib.loads((new / 'summary.toml').read_text())
old_rows = list(csv.DictReader((old / 'warmup.csv').open()))
new_rows = list(csv.DictReader((new / 'trace.csv').open()))
fig, axes = plt.subplots(4, 1, figsize=(10, 10), layout='constrained', sharex=True)
metrics = [('original_primal', 5e-7, 'Original feasibility / threshold', PRIMAL),
    ('r_max', 2.5e-7, 'All-copy consensus / threshold', PRIMAL),
    ('s_max', 1e-4, 'Dual residual / threshold', DUAL),
    ('rho', 1, 'Penalty rho (cost / (p.u. h)^2)', RHO)]
for rows, label, style in ((old_rows, 'Original harness: warm-up failed at 5000', '--'),
        (new_rows, f"Corrected harness: {summary['iterations']} iterations, validated", '-')):
    for axis, (key, scale, ylabel, color) in zip(axes, metrics):
        x = [int(r['iteration']) for r in rows]
        y = [float(r[key])/scale for r in rows]
        if key == 'rho':
            axis.step(x, y, where='post', color=color, linestyle=style, label=label)
            axis.set_xscale('log'); axis.set_yscale('log')
        else:
            axis.loglog(x, y, color=color, linestyle=style, label=label)
        axis.set_ylabel(ylabel)
        axis.grid(True, alpha=.2)
for axis in axes[:3]:
    axis.axhline(1, color=THRESHOLD, linestyle=':', linewidth=.8)
axes[0].legend(frameon=False, fontsize=9)
axes[-1].set_xlabel('ADMM iteration (log scale)')
fig.suptitle('med2522 T=24: source-voltage redundancy and local-accuracy fix\nSame frozen input, MA57, initial state and rho controller; corrected dispatch passes OpenDSS', fontsize=11)
fig.savefig(out / 'accuracy-fix-comparison.png', dpi=170)
plt.close(fig)
