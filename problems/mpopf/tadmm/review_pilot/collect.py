"""Package frozen pilot evidence; never turn an incomplete case into a result."""
import csv
import gzip
import json
from datetime import datetime
from pathlib import Path
import tarfile
import tomllib

pilot = Path(__file__).resolve().parent
root = pilot / 'results'
rows = []
ledger = []
for out in sorted(root.iterdir()):
    if not out.is_dir():
        continue
    if (out / 'metadata.toml').is_file() and (out / 'source.sha256').is_file():
        meta = tomllib.loads((out / 'metadata.toml').read_text())
        diagnostic = (out / 'diagnostic_status.json').exists()
        item = dict(configuration=out.name, linear_solver=meta['linear_solver'],
            status='interrupted diagnostic' if diagnostic else ('validated' if (out / 'validation.json').exists() else 'incomplete'))
        for stem in ('warmup', 'trace'):
            file = out / f'{stem}.csv'
            trace = list(csv.DictReader(file.open())) if file.exists() else []
            item[stem + '_iterations'] = int(trace[-1]['iteration']) if trace else 0
            item[stem + '_algorithm_s'] = float(trace[-1]['algorithm_s']) if trace else 0.0
        samples = list(csv.DictReader((out / 'process_samples.csv').open()))
        item['sampled_process_span_s'] = (datetime.fromisoformat(samples[-1]['utc'].replace('Z', '+00:00')) - datetime.fromisoformat(samples[0]['utc'].replace('Z', '+00:00'))).total_seconds() if len(samples) > 1 else 0
        ledger.append(item)
    if (out / 'source').is_dir():
        with tarfile.open(out / 'source_snapshot.tgz', 'w:gz') as archive:
            archive.add(out / 'source', arcname='source')
    for file in list(out.glob('*.log')) + list(out.glob('*.err')):
        (out / (file.name + '.gz')).write_bytes(gzip.compress(file.read_bytes(), mtime=0))
    if not out.name.endswith(('_phase', '_fixed_ma57', '_relative_ma57')) or not (out / 'validation.json').is_file():
        continue
    val = json.loads((out / 'validation.json').read_text())
    s = tomllib.loads((out / 'summary.toml').read_text())
    meta = tomllib.loads((out / 'metadata.toml').read_text())
    assert val['optimization_validation'] == val['opendss_validation'] == 'pass'
    row = {key: s[key] for key in ('system', 'T', 'iterations', 'central_objective',
        'tadmm_objective', 'relative_objective_gap', 'central_wall_s', 'tadmm_algorithm_s',
        'tadmm_wall_with_diagnostics_s', 'ideal_parallel_estimate_s', 'tadmm_original_primal',
        'tadmm_consensus_max', 'tadmm_dual_max', 'tadmm_soc_equality_gap')}
    row['opendss_voltage_difference_pu'] = val['tadmm']['max_voltage_difference_pu']
    row['sampled_process_peak_gib'] = val['sampled_process_peak_gib']
    row['evidence_directory'] = out.name
    row['linear_solver'] = meta['linear_solver']
    row['rho_policy'] = meta['rho_policy']
    row['rho_initial'] = meta.get('rho_initial', 3000 * (s['T'] / 24) ** 0.5)
    rows.append(row)
if rows:
    with (root / 'completed_pilots.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
if ledger:
    with (root / 'experiment_ledger.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=ledger[0].keys())
        writer.writeheader()
        writer.writerows(ledger)
print(json.dumps(rows, indent=2))
