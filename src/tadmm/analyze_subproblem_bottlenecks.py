#!/usr/bin/env python3
"""
Analyze which time periods are bottlenecks in tADMM subproblem solving.
Pure Python - no dependencies required.
"""

import csv
import statistics
from pathlib import Path
from collections import Counter

# Load data
csv_path = Path(r"c:\Users\arjha\OneDrive - Tesla\Documents\documents_general_addendum\MultiPeriodDistOPF\envs\tadmm\processedData\ieee2552_1ph_T12\subproblem_timing_details.csv")

data = []
with open(csv_path, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        data.append({
            'iteration': int(row['iteration']),
            'subproblem': int(row['subproblem']),
            'solve_time_sec': float(row['solve_time_sec'])
        })

# Group by iteration and find worst subproblem per iteration
iterations = {}
for row in data:
    k = row['iteration']
    if k not in iterations:
        iterations[k] = []
    iterations[k].append((row['subproblem'], row['solve_time_sec']))

# Find the worst (slowest) subproblem for each iteration
worst_per_iter = []
worst_times = []
for k in sorted(iterations.keys()):
    subproblems = iterations[k]
    worst_t0, worst_time = max(subproblems, key=lambda x: x[1])
    worst_per_iter.append(worst_t0)
    worst_times.append(worst_time)

# Count how often each time period is the bottleneck
bottleneck_counts = Counter(worst_per_iter)

# Print summary statistics
print("\n" + "="*70)
print("BOTTLENECK ANALYSIS: ieee2552_1ph, T=12")
print("="*70)
print("\nWhich time periods are slowest most often?")
print("-" * 70)
print(f"{'Time Period':>15} | {'Count':>8} | {'Percentage':>12}")
print("-" * 70)
for t0 in sorted(bottleneck_counts.keys()):
    pct = 100 * bottleneck_counts[t0] / len(worst_per_iter)
    print(f"  t0={t0:2d}           | {bottleneck_counts[t0]:4d}     | {pct:6.1f}%")
print("-" * 70)

print(f"\nTotal iterations analyzed: {len(worst_per_iter)}")

mean_time = statistics.mean(worst_times)
median_time = statistics.median(worst_times)
std_time = statistics.stdev(worst_times)

print(f"\nWorst time distribution (makespan per iteration):")
print(f"  Mean:   {mean_time:.3f}s")
print(f"  Median: {median_time:.3f}s")
print(f"  Min:    {min(worst_times):.3f}s")
print(f"  Max:    {max(worst_times):.3f}s")
print(f"  Std:    {std_time:.3f}s")

# Check if any time period is a persistent culprit (>2x expected frequency)
expected_freq = len(worst_per_iter) / 12
threshold = 2 * expected_freq
culprits = {t0: count for t0, count in bottleneck_counts.items() if count > threshold}

print(f"\nExpected frequency (uniform distribution): {expected_freq:.2f} times per time period")
print(f"Outlier threshold (>2x expected):          {threshold:.2f} times")

if len(culprits) > 0:
    print(f"\nWARNING: PERSISTENT BOTTLENECKS DETECTED:")
    for t0 in sorted(culprits.keys()):
        count = culprits[t0]
        pct = 100 * count / len(worst_per_iter)
        print(f"  t0={t0:2d}: {count} times ({pct:.1f}%) - {count/expected_freq:.2f}x expected")
else:
    print(f"\nNO PERSISTENT BOTTLENECKS")
    print(f"  All time periods are within 2x expected frequency")
    print(f"  Range: {min(bottleneck_counts.values())} to {max(bottleneck_counts.values())} times")
    spread = max(bottleneck_counts.values()) - min(bottleneck_counts.values())
    print(f"  Spread: {spread} occurrences ({100*spread/expected_freq:.1f}% of expected)")

print("="*70)
print("\nAnalysis complete")
