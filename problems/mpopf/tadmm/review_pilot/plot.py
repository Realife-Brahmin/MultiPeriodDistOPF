"""Static scientific convergence figure, only for completed validated pilots."""
import csv
from pathlib import Path
import tomllib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent / 'results'
fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout='constrained')
names = {'ieee123C_1ph': 'IEEE123', 'ieee2522C_1ph': 'med2522', 'large10kC_1ph': 'large10k'}
panels = [('gap', 'Objective agreement', 'Absolute relative objective gap', 1e-5),
          ('original_primal', 'Original-model feasibility', 'Maximum violation (per unit)', 5e-7),
          ('r_max', 'All SOC copies', 'Maximum consensus error (p.u. h)', 2.5e-7),
          ('s_max', 'ADMM stationarity', 'Maximum dual residual (objective / p.u. h)', 1e-4)]
found = False
for out in sorted(root.iterdir()):
    if not out.name.endswith(('_phase', '_relative_ma57')):
        continue
    if not (out / 'validation.json').is_file():
        continue
    found = True
    summary = tomllib.loads((out / 'summary.toml').read_text())
    trace = list(csv.DictReader((out / 'trace.csv').open()))
    iterations = [int(r['iteration']) for r in trace]
    ref = summary['central_objective']
    for axis, (key, title, ylabel, tol) in zip(axes.flat, panels):
        values = [abs(float(r['objective']) - ref) / max(1, abs(ref)) if key == 'gap' else float(r[key]) for r in trace]
        axis.semilogy(iterations, [max(v, 1e-16) for v in values], label=names[summary['system']], linewidth=1.6)
        axis.set(title=title, xlabel='ADMM iteration', ylabel=ylabel)
        axis.axhline(tol, color='0.6', linestyle='--', linewidth=0.7)
        axis.grid(True, which='major', alpha=0.2)
assert found, 'No validated measured traces available'
axes[0, 0].legend(frameon=False)
fig.suptitle('Corrected tADMM: free terminal SOC, T = 6, sequential local solves\nIEEE123/med2522: MUMPS, phase-aware rho; large10k: MA57, normalized residual balancing\nMeasured restarts after full warm-up; dashed lines show pilot validation tolerances', fontsize=10)
fig.savefig(root / 'convergence.png', dpi=180)
fig.savefig(root / 'convergence.svg')
plt.close(fig)
