"""Scientific plots of the frozen T=6 verification trajectories, not timings."""
import csv
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullLocator
from plot_style import PRIMAL, DUAL, RHO, THRESHOLD

campaign = Path(sys.argv[1]).resolve()
out = Path(sys.argv[2]).resolve()
fig, axes = plt.subplots(4, 2, figsize=(10, 11), layout='constrained')
metrics = [('original_primal', 'Original-model maximum violation (p.u.)', 5e-7),
    ('r_max', 'All-copy consensus maximum (p.u. h)', 2.5e-7),
    ('s_max', 'Dual residual maximum (cost / p.u. h)', 1e-4),
    ('rho', 'Penalty rho (cost / (p.u. h)^2)', None)]
for col, (system, name) in enumerate((('ieee2522C_1ph', 'med2522'), ('large10kC_1ph', 'large10k'))):
    for backend, style in (('ma57', '-'), ('ma97', '--')):
        cases = sorted((campaign / 'results').glob(f'{system}_T6_{backend}_verify_r0_*'))
        cases = [p for p in cases if (p / 'COMPLETE').exists()]
        assert len(cases) == 1, (system, backend, cases)
        trace = list(csv.DictReader((cases[0] / 'trace.csv').open()))
        for row, (key, label, threshold) in enumerate(metrics):
            axis = axes[row, col]
            color = DUAL if key == 's_max' else RHO if key == 'rho' else PRIMAL
            x = [int(r['iteration']) for r in trace]
            y = [max(1e-16, float(r[key])) for r in trace]
            if key == 'rho':
                axis.step(x, y, where='post', color=color, linestyle=style, linewidth=1.5, label=backend.upper())
                axis.set_xscale('log')
                axis.set_yscale('log', base=2)
            else:
                axis.loglog(x, y, color=color, linestyle=style, linewidth=1.5, label=backend.upper())
            axis.set(xlabel='ADMM iteration (log scale)', ylabel=label)
            axis.grid(True, which='major', alpha=.2)
    for row, (_, _, threshold) in enumerate(metrics):
        if threshold is not None:
            axes[row, col].axhline(threshold, color=THRESHOLD, linewidth=.8, linestyle=':')
    axes[0, col].set_title(f'{name}: corrected free-terminal MPOPF, T = 6')
    axes[0, col].legend(frameon=False)
    rho_axis = axes[3, col]
    levels = sorted(set(y for line in rho_axis.lines for y in line.get_ydata()))
    rho_axis.set_yticks(levels)
    rho_axis.yaxis.set_major_formatter(FuncFormatter(lambda value, pos: f'{value:g}'))
    rho_axis.yaxis.set_minor_locator(NullLocator())
fig.suptitle('Penalty and convergence diagnosis: same relative-residual controller\nMeasured restarts after full warm-up; dashed lines are full-convergence thresholds', fontsize=11)
fig.text(.5, -.012, 'Source: Kamiak frozen verification cases, 7 November 2026. Shared-node runs: diagnose trajectories, not benchmark timings.', ha='center', fontsize=8)
fig.savefig(out / 'penalty-convergence.png', dpi=180, bbox_inches='tight')
fig.savefig(out / 'penalty-convergence.svg', bbox_inches='tight')
plt.close(fig)
