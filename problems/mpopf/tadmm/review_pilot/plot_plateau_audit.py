"""Controlled saved-state replay: metric colors follow the repository palette."""
import csv
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import PRIMAL, DUAL, THRESHOLD

root = Path(sys.argv[1]).resolve()
series = [('plateau_audit_31645686/default', 'Original local accuracy', '--', None),
    ('plateau_audit_31645686/strict_comp', 'Explicit complementarity limit 1e-8', '-', 'o'),
    ('plateau_audit_31645728/root_bounds_removed', 'Redundant root voltage bounds removed', ':', 's')]
plt.rcParams.update({'font.size': 10})
fig, axes = plt.subplots(3, 1, figsize=(10, 8), layout='constrained', sharex=True)
metrics = [('original_primal', 5e-7, 'Original feasibility / threshold', PRIMAL),
    ('r_max', 2.5e-7, 'All-copy consensus / threshold', PRIMAL),
    ('s_max', 1e-4, 'Dual residual / threshold', DUAL)]
for rel, label, style, marker in series:
    rows = list(csv.DictReader((root / 'diagnostics' / rel / 'replay.csv').open()))
    for axis, (key, threshold, ylabel, color) in zip(axes, metrics):
        axis.semilogy([int(r['step']) for r in rows], [float(r[key])/threshold for r in rows],
            label=label, color=color, linestyle=style, marker=marker, markersize=4, linewidth=1.7)
        axis.set_ylabel(ylabel + '\n(dimensionless)')
for axis in axes:
    axis.axhline(1, color=THRESHOLD, linestyle=':', linewidth=.8)
    axis.grid(True, alpha=.2)
axes[0].legend(frameon=False, fontsize=9)
axes[-1].set_xlabel('Additional ADMM steps from the same iteration-5000 checkpoint')
fig.suptitle('med2522 T=24: damping restored while rho stays fixed at 7,600\nMA57 only; original consensus and dual updates; short diagnostic replay', fontsize=12)
fig.savefig(root / 'controlled-replay.png', dpi=170)
plt.close(fig)
