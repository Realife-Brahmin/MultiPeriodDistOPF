"""Sample the optimizer process and its descendants; Slurm peak is authoritative."""
import csv
from datetime import datetime, timezone
import os
from pathlib import Path
import sys
import time

pid = int(sys.argv[1])
ticks = os.sysconf('SC_CLK_TCK')
seen = {}
def descend(p):
    result = [p]
    try:
        children = Path(f'/proc/{p}/task/{p}/children').read_text().split()
        for child in children:
            result.extend(descend(int(child)))
    except (FileNotFoundError, ProcessLookupError, PermissionError):
        pass
    return result

with open(sys.argv[2], 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['utc', 'working_set_bytes', 'private_bytes', 'cpu_seconds', 'pss_bytes', 'process_count'])
    while Path(f'/proc/{pid}').exists():
        rss = private = pss = count = 0
        for p in descend(pid):
            try:
                status = dict(line.split(':', 1) for line in Path(f'/proc/{p}/status').read_text().splitlines() if ':' in line)
                if status.get('State', '').strip().startswith('Z'):
                    continue
                stat = Path(f'/proc/{p}/stat').read_text().rsplit(')', 1)[1].split()
                seen[(p, stat[19])] = (int(stat[11]) + int(stat[12])) / ticks
                rss += int(status.get('VmRSS', '0 kB').split()[0]) * 1024
                roll = dict(line.split(':', 1) for line in Path(f'/proc/{p}/smaps_rollup').read_text().splitlines() if ':' in line)
                private += sum(int(roll.get(k, '0 kB').split()[0]) * 1024 for k in ('Private_Clean', 'Private_Dirty'))
                pss += int(roll.get('Pss', '0 kB').split()[0]) * 1024
                count += 1
            except (FileNotFoundError, ProcessLookupError, PermissionError):
                continue
        writer.writerow([datetime.now(timezone.utc).isoformat(), rss, private, sum(seen.values()), pss, count])
        f.flush()
        time.sleep(0.5)
