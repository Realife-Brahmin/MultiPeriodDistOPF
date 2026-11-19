"""Administrative watcher for the two predeclared split T96 batches.

Does not modify numerical code or interrupt either requested first trajectory.
Run on the login node: two filesystem checks per minute, no numerical work.
"""
import datetime
import json
from pathlib import Path
import subprocess
import sys
import time

campaign = Path(sys.argv[1]).resolve()
pending = {31652342: 60000, 31653264: 240000}
deadline = time.monotonic() + 9 * 3600
while pending and time.monotonic() < deadline:
    for job, rho in list(pending.items()):
        case = campaign / 'results' / f'large10k_T96_{job}'
        validation = case / f'validation_rho{rho}.json'
        if validation.exists():
            try:
                record = json.loads(validation.read_text())
                arm = record['arms'][f'rho{rho}']
                assert arm['opendss']['endpoint.jls']['converged']
                assert (case / f'rho{rho}' / 'assessment.toml').exists()
            except (ValueError, KeyError, AssertionError):
                continue
            note = dict(job=job, selected_rho=rho,
                        utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                        reason='Selected first arm fully validated; prevent duplicate second arm. '
                               'This cancellation is not a solver convergence flag.',
                        validation=str(validation), arm=arm)
            dest = campaign / 'campaign' / f'stop-after-first-{job}.json'
            dest.write_text(json.dumps(note, indent=2) + '\n')
            result = subprocess.run(['scancel', str(job)], capture_output=True, text=True)
            note.update(scancel_returncode=result.returncode, scancel_stdout=result.stdout,
                        scancel_stderr=result.stderr)
            dest.write_text(json.dumps(note, indent=2) + '\n')
            print(json.dumps(note), flush=True)
            if result.returncode == 0:
                del pending[job]
    if pending:
        time.sleep(60)
if pending:
    raise SystemExit(f'Watcher expired without cancelling incomplete selected arms: {pending}')
