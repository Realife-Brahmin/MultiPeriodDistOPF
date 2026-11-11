"""Derive per-round timing distributions from the measured trace files."""
import csv
import math
from pathlib import Path
import statistics
import sys
from collections import defaultdict

arm = Path(sys.argv[1])
periods_by_round = defaultdict(list)
with (arm / 'trace_subproblems.csv').open(newline='') as source:
    for row in csv.DictReader(source):
        periods_by_round[int(row['iteration'])].append(row)

with (arm / 'trace.csv').open(newline='') as source:
    rounds = list(csv.DictReader(source))
assert rounds and len(rounds) == len(periods_by_round)
expected_periods = set(range(1, len(periods_by_round[1]) + 1))

derived = []
previous = {'algorithm_s': 0.0, 'loop_wall_s': 0.0, 'solver_sum_s': 0.0}
for round_row in rounds:
    iteration = int(round_row['iteration'])
    workers = periods_by_round[iteration]
    assert len(workers) == len(expected_periods)
    assert {int(row['period']) for row in workers} == expected_periods
    values = dict(round_row)
    for key, before in previous.items():
        cumulative = float(round_row[key])
        assert math.isfinite(cumulative) and cumulative >= before
        values[f'round_{key}'] = cumulative - before
        previous[key] = cumulative
    for key in ('solver_s', 'stage_elapsed_s', 'ipopt_iterations'):
        samples = [float(row[key]) for row in workers]
        assert all(math.isfinite(value) and value >= 0 for value in samples)
        values[f'local_min_{key}'] = min(samples)
        values[f'local_median_{key}'] = statistics.median(samples)
        values[f'local_mean_{key}'] = math.fsum(samples) / len(samples)
        values[f'local_max_{key}'] = max(samples)
    slowest = max(workers, key=lambda row: float(row['stage_elapsed_s']))
    values['slowest_period'] = int(slowest['period'])
    values['slowest_period_ipopt_iterations'] = int(slowest['ipopt_iterations'])
    values['capped_local_solves'] = sum(row['iteration_limit_hit'] == 'true' for row in workers)
    derived.append(values)

with (arm / 'round-timing.csv').open('w', newline='') as output:
    writer = csv.DictWriter(output, fieldnames=list(derived[0]))
    writer.writeheader()
    writer.writerows(derived)
print(f'Summarized {len(derived)} rounds with {len(expected_periods)} periods each.')
