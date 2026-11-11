"""Summarize a completed compact Kamiak tuning archive."""
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile
import tomllib

for argument in sys.argv[1:]:
    archive = Path(argument)
    with archive.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    with tarfile.open(archive, 'r:gz') as evidence:
        root = next(name.removesuffix('sweep.toml') for name in evidence.getnames() if name.endswith('/sweep.toml'))

        def read(name):
            with evidence.extractfile(root + name) as stream:
                return stream.read().decode()

        meta = tomllib.loads(read('sweep.toml'))
        center = meta['central_wall_s']
        result = {
            'archive_sha256': digest,
            'central_s': center,
            'settings': {key: meta[key] for key in ('local_solve_policy', 'cap_iterations', 'rho_floor', 'bootstrap_decrease')},
            'arms': {},
            'memory': json.loads(read('posthoc-memory-peaks.json')),
        }
        for rho in meta['rhos']:
            name = f'rho{int(rho)}'
            assessment = tomllib.loads(read(name + '/assessment.toml'))
            milestones = tomllib.loads(read(name + '/milestones.toml'))
            audit = json.loads(read(name + '_endpoint_sequential/socp-storage-audit.json'))
            trace = list(csv.DictReader(io.StringIO(read(name + '/round-timing.csv'))))
            assert len(trace) == assessment['iterations']
            within_budget = [row for row in trace if float(row['loop_wall_s']) <= center / 2]
            physical = audit['physical_error_statistics']
            result['arms'][str(int(rho))] = {
                **{key: assessment[key] for key in ('iterations', 'core_s', 'loop_wall_s', 'original_primal', 'objective_gap', 'sequential_schedule_pass')},
                'physical_NO_milestone': milestones.get('physical_NO'),
                'at_half_central_wall': within_budget[-1] if within_budget else None,
                'cap_hits': sum(int(row['capped_local_solves']) for row in trace),
                'cone': audit['SOCP'],
                'battery_maxima': {
                    key: {'maximum': row['maximum'], 'worst_location': row['worst_location']}
                    for key, row in physical['battery'].items()
                },
                'voltage_maximum_difference_pu': physical['bus_voltage']['voltage_agreement_pu']['maximum'],
            }
    archive.with_suffix('.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(archive.name, 'central', center, [(rho, arm['core_s'], arm['sequential_schedule_pass']) for rho, arm in result['arms'].items()])
