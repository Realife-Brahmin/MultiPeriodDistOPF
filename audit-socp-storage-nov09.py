"""Independent table audit; never participates in optimization timing.

Usage: python audit-socp-storage-nov09.py REPLAY_FOLDER Storage.dss
Writes explicit per-branch tightness and rating-normalized storage errors.
"""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys
try:
    import tomllib
except ImportError:
    import tomli as tomllib

folder, storage = map(Path, sys.argv[1:3])
summary = tomllib.loads((folder / 'summary.toml').read_text())
base = summary.get('base_kva', 1000.0)
ratings = {}
for line in storage.read_text().splitlines():
    if not re.match(r'\s*new\s+storage\.', line, re.I):
        continue
    fields = dict((k.lower(), v.strip('"')) for k, v in re.findall(r'([\w%]+)\s*=\s*(\S+)', line))
    bus = int(fields['bus1'].split('.')[0])
    ratings[bus] = (float(fields['kwhrated']), float(fields['kwrated']))

def rows(name):
    path = folder / (name + '.csv.gz')
    raw = gzip.open(path, 'rb').read() if path.exists() else (folder / (name + '.csv')).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == summary['table_sha256'][name + '.csv']
    return csv.DictReader(raw.decode().splitlines())

voltage = {(int(r['t']), int(r['bus'])): float(r['opf_voltage_pu']) for r in rows('bus')}
abs_gaps, rel_gaps, negative = [], [], 0.0
worst_cone = None
with gzip.open(folder / 'socp-tightness.csv.gz', 'wt', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['t', 'from_bus', 'to_bus', 'v_ell_pu2', 'P2_Q2_pu2', 'signed_gap_pu2', 'relative_gap', 'relative_denominator_floor_pu2'])
    for r in rows('branch'):
        t, bus = int(r['t']), int(r['from_bus'])
        lhs = voltage[t, bus] ** 2 * float(r['opf_ell_pu'])
        rhs = (float(r['opf_P_kw']) / base) ** 2 + (float(r['opf_Q_kvar']) / base) ** 2
        gap = lhs - rhs
        assert all(math.isfinite(x) for x in (lhs, rhs, gap))
        relative = abs(gap) / max(abs(lhs), rhs, 1e-12)
        abs_gaps.append(abs(gap)); rel_gaps.append(relative)
        negative = max(negative, -gap)
        record = [t, bus, int(r['to_bus']), lhs, rhs, gap, relative, 1e-12]
        writer.writerow(record)
        if worst_cone is None or abs(gap) > abs(worst_cone[5]):
            worst_cone = record

maxima = {k: 0.0 for k in ('SOC_state_error_pp', 'energy_accounting_error_pct_capacity', 'requested_bound_error_pp', 'engine_bound_error_pp', 'dispatch_error_pct_rated_power')}
worst_storage = {}
for r in rows('battery'):
    capacity, power = ratings[int(r['bus'])]
    assert capacity > 0 and power > 0
    values = {
        'SOC_state_error_pp': abs(float(r['engine_minus_opf_kwh'])) / capacity * 100,
        'energy_accounting_error_pct_capacity': abs(float(r['engine_energy_balance_error_kwh'])) / capacity * 100,
        'requested_bound_error_pp': float(r['requested_bound_error_kwh']) / capacity * 100,
        'engine_bound_error_pp': float(r['engine_bound_error_kwh']) / capacity * 100,
        'dispatch_error_pct_rated_power': abs(float(r['power_difference_kw'])) / power * 100,
    }
    for key, value in values.items():
        assert math.isfinite(value)
        if value >= maxima[key]:
            maxima[key] = value
            worst_storage[key] = {'t': int(r['t']), 'bus': int(r['bus']), 'capacity_kwh': capacity, 'power_kw': power}

def percentile(values, p):
    data = sorted(values)
    x = (len(data) - 1) * p
    a, b = math.floor(x), math.ceil(x)
    return data[a] + (data[b] - data[a]) * (x - a)

result = {
    'source': str(folder), 'summary_sha256': hashlib.sha256((folder / 'summary.toml').read_bytes()).hexdigest(),
    'ratings_sha256': hashlib.sha256(storage.read_bytes()).hexdigest(), 'base_kva': base,
    'SOCP': {'definition': 'v_i_squared*ell_ij - P_ij_squared - Q_ij_squared, per-unit squared',
             'max_absolute_gap_pu2': max(abs_gaps), 'p95_absolute_gap_pu2': percentile(abs_gaps, .95),
             'max_cone_inequality_violation_pu2': negative, 'max_relative_gap': max(rel_gaps),
             'p95_relative_gap': percentile(rel_gaps, .95), 'relative_denominator_floor_pu2': 1e-12,
             'absolute_tightness_target_pu2': 1e-5, 'absolute_tightness_pass': max(abs_gaps) <= 1e-5,
             'worst_record': worst_cone, 'branch_period_records': len(abs_gaps),
             'note': 'Inequality feasibility alone does not establish tightness. Relative gaps near zero flow are floor-sensitive; inspect absolute and relative results together.'},
    'storage_maxima': maxima, 'storage_worst_locations': worst_storage,
    'engineering_target': {'SOC_pp_and_dispatch_pct_rated_target': .05, 'concerning_level': .5,
                           'all_storage_metrics_within_target': max(maxima.values()) <= .05,
                           'note': 'Engineering agreement target is separate from unchanged strict physical feasibility and sequential schedule flags.'},
    'strict_sequential_schedule_pass': summary['sequential_schedule_pass'],
}
(folder / 'socp-storage-audit.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
