"""Read-only, immutable diagnostic snapshots of the two current T=24 runs."""
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import sys

campaign, out = map(lambda p: Path(p).resolve(), sys.argv[1:3])
out.mkdir(parents=True, exist_ok=False)
cases = [('med2522', 'ieee2522C_1ph_T24_ma57_benchmark_r1_31643458_6_a0'),
         ('large10k', 'large10kC_1ph_T24_ma57_firstpass_31643647_a0')]
report = dict(snapshot_utc=datetime.now(timezone.utc).isoformat(), cases=[])
for name, dirname in cases:
    case = campaign / 'results' / dirname
    dest = out / name
    dest.mkdir()
    for filename in ('metadata.toml', 'source.sha256', 'input.sha256', 'solve.log',
            'solve.err', 'central_trace.csv', 'phase_events.csv', 'central_ipopt_timed.log'):
        file = case / filename
        if file.exists():
            (dest / filename).write_bytes(file.read_bytes())
    # Keep only whole flushed lines; never infer completion from a live CSV.
    raw = (case / 'warmup.csv').read_bytes()
    raw = raw[:raw.rfind(b'\n') + 1]
    rows = list(csv.DictReader(io.StringIO(raw.decode())))
    assert rows and all(all(v is not None for v in r.values()) for r in rows)
    (dest / 'warmup.csv').write_bytes(raw)
    text = (dest / 'solve.log').read_text()
    match = re.search(r'CENTRAL_DONE objective=([\deE+.-]+)', text)
    assert match, name
    item = dict(name=name, source_case=dirname, phase='warmup',
        iterations=len(rows), last_iteration=int(rows[-1]['iteration']),
        central_objective=float(match[1]),
        completed_case_at_snapshot=(case / 'COMPLETE').exists(),
        note='Partial initial warm-up trajectory; not a completed measured run or independently validated dispatch')
    done = re.search(r'TADMM_DONE phase=warmup converged=(true|false) iterations=(\d+)', text)
    item['warmup_status'] = ('converged' if done[1] == 'true' else 'not converged') if done else 'in progress'
    item['note'] = 'Initial warm-up trajectory; measured-run and independent validation completion must be checked separately'
    report['cases'].append(item)
hashes = []
for file in sorted(out.rglob('*')):
    if file.is_file():
        hashes.append(hashlib.sha256(file.read_bytes()).hexdigest() + '  ' + file.relative_to(out).as_posix())
(out / 'snapshot.sha256').write_text('\n'.join(hashes) + '\n')
(out / 'snapshot.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
