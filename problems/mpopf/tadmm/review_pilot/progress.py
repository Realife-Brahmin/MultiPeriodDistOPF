import csv
from pathlib import Path
root = Path(__file__).resolve().parent / 'results'
for out in sorted(root.iterdir()):
    if not out.is_dir() or not out.name.endswith(('_phase', '_fixed_ma57', '_relative_ma57')) or (out / 'diagnostic_status.json').exists():
        continue
    if (out / 'validation.json').exists():
        continue
    trace = out / 'trace.csv'
    phase = 'measured'
    if not trace.exists():
        trace = out / 'warmup.csv'
        phase = 'warm-up'
    if not trace.exists():
        print(out.name, 'preparing centralized reference')
        continue
    rows = list(csv.DictReader(trace.open()))
    if rows:
        r = rows[-1]
        print(f"{out.name}: {phase} iteration {r['iteration']}, original-model residual {float(r['original_primal']):.2e}, dual residual {float(r['s_max']):.2e}")
