"""Summarize raw Ipopt termination tables and controlled replay diagnostics."""
import csv
import json
from pathlib import Path
import re
import sys

root = Path(sys.argv[1]).resolve()
rows, summary = [], {}
labels = {'dual_unscaled': 'Dual infeasibility', 'constraint_unscaled': 'Constraint violation',
          'complementarity_unscaled': 'Complementarity', 'nlp_error_unscaled': 'Overall NLP error'}
for variant in sorted(p for p in root.iterdir() if p.is_dir()):
    local = []
    for file in sorted(variant.glob('step_*_period_*.log')):
        match = re.fullmatch(r'step_(\d+)_period_(\d+)\.log', file.name)
        text = file.read_text()
        row = dict(variant=variant.name, step=int(match[1]), period=int(match[2]))
        for key, label in labels.items():
            m = re.search(r'^' + re.escape(label) + r'\.+:\s+(\S+)\s+(\S+)', text, re.MULTILINE)
            row[key] = float(m[2]) if m else None
        exit_match = re.search(r'^EXIT: (.*)$', text, re.MULTILINE)
        row['exit'] = exit_match[1] if exit_match else 'missing termination table'
        local.append(row)
    rows.extend(local)
    replay = list(csv.DictReader((variant / 'replay.csv').open()))
    item = dict(replayed_steps=len(replay), first=replay[0] if replay else None,
        last=replay[-1] if replay else None, locally_rejected=(variant / 'failure.toml').exists())
    for key in labels:
        candidates = [r for r in local if r[key] is not None]
        if candidates:
            worst = max(candidates, key=lambda r: r[key])
            item['max_' + key] = dict(value=worst[key], step=worst['step'], period=worst['period'])
    summary[variant.name] = item
if rows:
    with (root / 'local_kkt.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
(root / 'audit_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
