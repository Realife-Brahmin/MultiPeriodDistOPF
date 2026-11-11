"""Independent table audit; never participates in optimization timing.

Usage: python audit-socp-storage-nov09.py REPLAY_FOLDER Storage.dss [REPORT_NAME]
Writes explicit per-branch tightness and physical error distributions.
"""
import csv
from array import array
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
report_name = sys.argv[3] if len(sys.argv) > 3 else 'socp-storage-audit.json'
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

class ErrorStats:
    def __init__(self):
        self.values = {}
        self.worst = {}

    def add(self, location, **errors):
        for name, value in errors.items():
            assert math.isfinite(value) and value >= 0, (name, location, value)
            if name not in self.values:
                self.values[name] = array('d')
                self.worst[name] = (value, location)
            self.values[name].append(value)
            if value > self.worst[name][0]:
                self.worst[name] = (value, location)

    def summarize(self):
        result = {}
        for name, values in self.values.items():
            ordered = sorted(values)
            result[name] = {
                'samples': len(ordered),
                'nonzero_samples': sum(value > 0 for value in ordered),
                'mean': math.fsum(ordered) / len(ordered),
                'median': percentile_sorted(ordered, .5),
                'p95': percentile_sorted(ordered, .95),
                'maximum': ordered[-1],
                'worst_location': self.worst[name][1],
            }
        return result


voltage = {}
voltage_stats = ErrorStats()
for r in rows('bus'):
    t, bus = int(r['t']), int(r['bus'])
    opf, dss = float(r['opf_voltage_pu']), float(r['dss_voltage_pu'])
    lower, upper = float(r['lower_pu']), float(r['upper_pu'])
    voltage[t, bus] = opf
    voltage_stats.add(
        {'t': t, 'bus': bus}, voltage_agreement_pu=abs(opf - dss),
        opf_undervoltage_pu=max(0.0, lower - opf),
        opf_overvoltage_pu=max(0.0, opf - upper),
        dss_undervoltage_pu=max(0.0, lower - dss),
        dss_overvoltage_pu=max(0.0, dss - upper),
    )
abs_gaps, rel_gaps, negative = [], [], 0.0
worst_cone = None
branch_stats = ErrorStats()
with gzip.open(folder / 'socp-tightness.csv.gz', 'wt', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['t', 'from_bus', 'to_bus', 'v_ell_pu2', 'P2_Q2_pu2', 'signed_gap_pu2', 'relative_gap', 'relative_denominator_floor_pu2'])
    for r in rows('branch'):
        t, bus = int(r['t']), int(r['from_bus'])
        to_bus = int(r['to_bus'])
        opf_ell, dss_ell = float(r['opf_ell_pu']), float(r['dss_ell_pu'])
        opf_current = math.sqrt(max(0.0, opf_ell))
        dss_current = math.sqrt(max(0.0, dss_ell))
        lhs = voltage[t, bus] ** 2 * opf_ell
        rhs = (float(r['opf_P_kw']) / base) ** 2 + (float(r['opf_Q_kvar']) / base) ** 2
        gap = lhs - rhs
        assert all(math.isfinite(x) for x in (lhs, rhs, gap))
        relative = abs(gap) / max(abs(lhs), rhs, 1e-12)
        branch_stats.add(
            {'t': t, 'from_bus': bus, 'to_bus': to_bus},
            cone_absolute_equality_gap_pu2=abs(gap),
            cone_positive_relaxation_slack_pu2=max(0.0, gap),
            cone_inequality_violation_pu2=max(0.0, -gap),
            cone_gap_normalized_by_at_least_one=abs(gap) / max(1.0, abs(lhs), rhs),
            cone_gap_relative_to_flow=relative,
            active_power_agreement_kw=abs(float(r['P_difference_kw'])),
            reactive_power_agreement_kvar=abs(float(r['Q_difference_kvar'])),
            squared_current_agreement_pu=abs(opf_ell - dss_ell),
            opf_current_magnitude_pu=opf_current,
            dss_current_magnitude_pu=dss_current,
            current_magnitude_agreement_pu=abs(opf_current - dss_current),
            opf_squared_current_lower_bound_violation_pu=max(0.0, -opf_ell),
        )
        abs_gaps.append(abs(gap)); rel_gaps.append(relative)
        negative = max(negative, -gap)
        record = [t, bus, to_bus, lhs, rhs, gap, relative, 1e-12]
        writer.writerow(record)
        if worst_cone is None or abs(gap) > abs(worst_cone[5]):
            worst_cone = record

maxima = {k: 0.0 for k in ('SOC_state_error_pp', 'energy_accounting_error_pct_capacity', 'requested_bound_error_pp', 'engine_bound_error_pp', 'dispatch_error_pct_rated_power')}
worst_storage = {}
storage_stats = ErrorStats()
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
    location = {'t': int(r['t']), 'bus': int(r['bus']), 'capacity_kwh': capacity, 'power_kw': power}
    opf_energy = float(r['opf_energy_kwh'])
    opf_bound = max(0.0, float(r['SOC_low_kwh']) - opf_energy, opf_energy - float(r['SOC_high_kwh']))
    integrated_error = abs(float(r['requested_model_energy_kwh']) - opf_energy)
    opf_power_bound = max(0.0, abs(float(r['command_kw'])) - power)
    dss_power_bound = max(0.0, abs(float(r['actual_kw'])) - power)
    storage_stats.add(
        location, **values,
        SOC_state_error_kwh=abs(float(r['engine_minus_opf_kwh'])),
        energy_accounting_error_kwh=abs(float(r['engine_energy_balance_error_kwh'])),
        requested_bound_error_kwh=float(r['requested_bound_error_kwh']),
        engine_bound_error_kwh=float(r['engine_bound_error_kwh']),
        opf_bound_error_kwh=opf_bound, opf_bound_error_pp=100 * opf_bound / capacity,
        integrated_dispatch_vs_opf_energy_kwh=integrated_error,
        integrated_dispatch_vs_opf_SOC_pp=100 * integrated_error / capacity,
        dispatch_error_kw=abs(float(r['power_difference_kw'])),
        opf_dispatch_bound_error_kw=opf_power_bound,
        opf_dispatch_bound_error_pct_rated=100 * opf_power_bound / power,
        dss_dispatch_bound_error_kw=dss_power_bound,
        dss_dispatch_bound_error_pct_rated=100 * dss_power_bound / power,
    )
    for key, value in values.items():
        assert math.isfinite(value)
        if value >= maxima[key]:
            maxima[key] = value
            worst_storage[key] = {'t': int(r['t']), 'bus': int(r['bus']), 'capacity_kwh': capacity, 'power_kw': power}

def percentile_sorted(data, p):
    x = (len(data) - 1) * p
    a, b = math.floor(x), math.ceil(x)
    return data[a] + (data[b] - data[a]) * (x - a)


system_stats = ErrorStats()
for r in rows('period'):
    system_stats.add(
        {'t': int(r['t'])},
        substation_active_power_agreement_kw=abs(float(r['substation_difference_kw'])),
        substation_reactive_power_agreement_kvar=abs(float(r['dss_substation_kvar']) - float(r['opf_substation_kvar'])),
        line_loss_agreement_kw=abs(float(r['dss_line_loss_kw']) - float(r['opf_line_loss_kw'])),
        total_battery_dispatch_agreement_kw=abs(float(r['actual_battery_kw']) - float(r['command_battery_kw'])),
    )

result = {
    'source': str(folder), 'summary_sha256': hashlib.sha256((folder / 'summary.toml').read_bytes()).hexdigest(),
    'ratings_sha256': hashlib.sha256(storage.read_bytes()).hexdigest(), 'base_kva': base,
    'SOCP': {'definition': 'v_i_squared*ell_ij - P_ij_squared - Q_ij_squared, per-unit squared',
             'max_absolute_gap_pu2': max(abs_gaps), 'p95_absolute_gap_pu2': percentile_sorted(sorted(abs_gaps), .95),
             'max_cone_inequality_violation_pu2': negative, 'max_relative_gap': max(rel_gaps),
             'p95_relative_gap': percentile_sorted(sorted(rel_gaps), .95), 'relative_denominator_floor_pu2': 1e-12,
             'absolute_tightness_target_pu2': 1e-5, 'absolute_tightness_pass': max(abs_gaps) <= 1e-5,
             'worst_record': worst_cone, 'branch_period_records': len(abs_gaps),
             'note': 'Inequality feasibility alone does not establish tightness. Relative gaps near zero flow are floor-sensitive; inspect absolute and relative results together.'},
    'storage_maxima': maxima, 'storage_worst_locations': worst_storage,
    'physical_error_statistics': {
        'aggregation': 'Every recorded device and period, including zero errors. Nonzero counts include floating-point roundoff.',
        'battery': storage_stats.summarize(),
        'bus_voltage': voltage_stats.summarize(),
        'branch': branch_stats.summarize(),
        'system': system_stats.summarize(),
        'thermal_ampacity_violation': None,
        'thermal_ampacity_scope': 'The pilot has no modeled thermal upper bounds. Current agreement and the squared-current lower bound are measured.',
    },
    'engineering_target': {'SOC_pp_and_dispatch_pct_rated_target': .05, 'concerning_level': .5,
                           'all_storage_metrics_within_target': max(maxima.values()) <= .05,
                           'note': 'Engineering agreement target is separate from unchanged strict physical feasibility and sequential schedule flags.'},
    'strict_sequential_schedule_pass': summary['sequential_schedule_pass'],
}
(folder / report_name).write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
