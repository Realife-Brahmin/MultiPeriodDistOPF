import csv
import json
from pathlib import Path
import sys
try:
    import tomllib
except ImportError:
    import tomli as tomllib

root = Path(sys.argv[1])
serial = tomllib.loads((root / 'results/ieee123C_1ph_T6_ma57_serial_gate/summary.toml').read_text())
results = {}
for solver in ('ma57', 'ma97'):
    out = root / f'results/ieee123C_1ph_T6_{solver}_parallel_gate'
    assert (out / 'COMPLETE').is_file()
    par = tomllib.loads((out / 'summary.toml').read_text())
    gap = abs(par['tadmm_objective'] - serial['tadmm_objective']) / max(1, abs(serial['tadmm_objective']))
    assert gap <= 1e-7
    val = json.loads((out / 'validation.json').read_text())
    assert val['worker_layout_validation'] == val['opendss_validation'] == 'pass'
    results[solver] = dict(relative_objective_difference_from_serial=gap,
        serial_iterations=serial['iterations'], parallel_iterations=par['iterations'],
        workers=par['period_workers'], validation='pass')
print(json.dumps(results, indent=2))
