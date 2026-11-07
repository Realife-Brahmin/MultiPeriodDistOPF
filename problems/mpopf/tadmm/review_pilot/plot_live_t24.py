"""Plot immutable partial-run snapshots, with no claim of full validation."""
import csv
import json
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullLocator
import numpy as np
from plot_style import OBJECTIVE, PRIMAL, DUAL, RHO, CENTRAL_REFERENCE, THRESHOLD

root = Path(sys.argv[1]).resolve()
meta = json.loads((root / 'snapshot.json').read_text())
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9})
fig, axes = plt.subplots(4, 2, figsize=(12, 10), layout='constrained', sharex='col')
diagnoses = []
for col, case in enumerate(meta['cases']):
    rows = list(csv.DictReader((root / case['name'] / 'warmup.csv').open()))
    x = np.array([int(r['iteration']) for r in rows])
    values = lambda key: np.array([float(r[key]) for r in rows])
    ref = case['central_objective']
    original, rmax, smax = (values(k) for k in ('original_primal', 'r_max', 's_max'))
    objective, rho = values('objective'), values('rho')
    for data, threshold, label, color, style in ((original, 5e-7, 'Original feasibility / 5e-7', PRIMAL, '-'),
            (rmax, 2.5e-7, 'All-copy consensus / 2.5e-7', PRIMAL, '--'),
            (smax, 1e-4, 'Dual residual / 1e-4', DUAL, '-')):
        axes[0, col].loglog(x, np.maximum(data / threshold, 1e-16), label=label, color=color, linestyle=style)
    axes[0, col].axhline(1, color=THRESHOLD, linestyle='--', linewidth=.8)
    axes[0, col].set(title=f"{case['name']}: T=24, MA57; k={x[-1]}, {case.get('warmup_status', 'snapshot')}")
    axes[0, col].legend(frameon=False, fontsize=8)
    gap = np.abs(objective - ref) / max(1, abs(ref))
    axes[1, col].loglog(x, np.maximum(gap, 1e-16), color=OBJECTIVE)
    axes[1, col].axhline(.001, color=THRESHOLD, linestyle='--', linewidth=.8, label='0.1% objective criterion')
    axes[1, col].legend(frameon=False, fontsize=8)
    axes[2, col].step(x, rho, where='post', color=RHO)
    axes[2, col].set_xscale('log')
    axes[2, col].set_yscale('log', base=2)
    axes[2, col].set_yticks(sorted(set(rho)))
    axes[2, col].yaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
    axes[2, col].yaxis.set_minor_locator(NullLocator())
    relative_ratio = values('r_relative') / np.maximum(values('s_relative'), 1e-16)
    axes[3, col].loglog(x, np.maximum(relative_ratio, 1e-16), color=RHO)
    axes[3, col].axhspan(.1, 10, color='.9', zorder=-1)
    for boundary in (.1, 10):
        axes[3, col].axhline(boundary, color=THRESHOLD, linestyle='--', linewidth=.8)
    axes[3, col].text(.03, .93, 'Shaded: no rho update (0.1 to 10); checked every 5 iterations', transform=axes[3, col].transAxes, fontsize=8)
    for row in range(4):
        axes[row, col].grid(True, which='major', alpha=.2)
    axes[3, col].set_xlabel('ADMM iteration (log scale)')
    changes = []
    for i in np.where(rho[1:] != rho[:-1])[0] + 1:
        changes.append(dict(decision_after_iteration=int(x[i-1]), rho_old=float(rho[i-1]), rho_new=float(rho[i])))
    n = min(100, len(x))
    tail = slice(-n, None)
    near = (gap <= .001) & (original <= 1e-5)
    bad = np.where(~near)[0]
    sustained = int(bad[-1] + 1) if len(bad) else 0
    diagnoses.append(dict(system=case['name'], snapshot_iteration=int(x[-1]),
        rho_changes=changes, rho_final=float(rho[-1]),
        first_near_optimal_iteration=int(x[np.where(near)[0][0]]) if np.any(near) else None,
        near_optimal_sustained_through_snapshot_from=int(x[sustained]) if sustained < len(x) else None,
        tail_iterations=n, last_original=float(original[-1]), last_consensus=float(rmax[-1]), last_dual=float(smax[-1]),
        tail_original_min=float(min(original[tail])), tail_original_max=float(max(original[tail])),
        tail_dual_min=float(min(smax[tail])), tail_dual_max=float(max(smax[tail])),
        objective_gap_last=float(gap[-1]), relative_ratio_last=float(relative_ratio[-1]),
        loop_wall_s=float(rows[-1]['loop_wall_s']),
        warmup_status=case.get('warmup_status', 'snapshot'),
        note='Warm-up diagnosis only; sustained qualification refers to this saved trajectory, not a separate measured run or OpenDSS validation'))
    if case['name'] == 'med2522':
        zoom, zz = plt.subplots(3, 1, figsize=(10, 7), layout='constrained', sharex=True)
        for ax, y, label, color in zip(zz, (objective-ref, original, smax/1e-4),
                ('Objective minus reference\n(objective units)', 'Original-model violation\n(per-unit units)', 'Dual residual / threshold\n(dimensionless)'),
                (OBJECTIVE, PRIMAL, DUAL)):
            ax.plot(x[tail], y[tail], color=color, linewidth=.8, alpha=.6)
            for parity, marker, filled, label_parity in ((1, 'o', True, 'Odd iterations'), (0, 's', False, 'Even iterations')):
                mask = x[tail] % 2 == parity
                ax.scatter(x[tail][mask], y[tail][mask], s=14, marker=marker,
                    facecolors=color if filled else 'white', edgecolors=color, linewidths=.8, label=label_parity)
            ax.set_ylabel(label)
            ax.grid(True, alpha=.2)
        zz[0].axhline(0, color=CENTRAL_REFERENCE, linewidth=.8, linestyle='--')
        zz[0].legend(frameon=False, ncol=2)
        zz[2].axhline(1, color=THRESHOLD, linewidth=.8, linestyle='--')
        zz[2].set_ylim(0, max(2, max(smax[tail]/1e-4)*1.2))
        zz[-1].set_xlabel('ADMM iteration')
        zoom.suptitle(f'med2522 T=24: last {n} warm-up iterations; rho = {rho[-1]:g}\nObserved alternation is a diagnostic, not proof of its cause', fontsize=11)
        zoom.savefig(root / 'med2522-tail.png', dpi=170)
        plt.close(zoom)
for row, label in enumerate(('Residual / threshold\n(dimensionless)', 'Absolute relative\nobjective gap', 'Penalty rho\n(cost / (p.u. h)^2)', 'Relative primal /\nrelative dual')):
    axes[row, 0].set_ylabel(label)
fig.suptitle('T=24 tADMM: initial convergence passes and penalty-controller behavior\nSaved warm-up trajectories; measured timing and OpenDSS completion are separate', fontsize=11)
fig.savefig(root / 'T24-convergence.png', dpi=170)
plt.close(fig)
(root / 'diagnosis.json').write_text(json.dumps(diagnoses, indent=2) + '\n')
print(json.dumps(diagnoses, indent=2))
