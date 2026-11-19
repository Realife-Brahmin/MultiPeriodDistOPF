"""Node-local PSS/RSS/CPU sampling of explicitly registered compute processes."""
import csv
from datetime import datetime, timezone
import os
from pathlib import Path
import socket
import sys
import time

folder=Path(sys.argv[1]);host=socket.gethostname();seen={};ticks=os.sysconf('SC_CLK_TCK')
with (folder.parent/f'process_samples.{host}.csv').open('w',newline='') as f:
    writer=csv.writer(f)
    writer.writerow(['utc','working_set_bytes','private_bytes','cpu_seconds','pss_bytes','process_count'])
    while not (folder.parent/'SAMPLING_DONE').exists():
        rss=private=pss=count=0
        for entry in folder.glob(f'{host}_*.pid'):
            try:
                pid=int(entry.read_text());base=Path(f'/proc/{pid}')
                status=dict(line.split(':',1) for line in (base/'status').read_text().splitlines() if ':' in line)
                if status.get('State','').strip().startswith('Z'):continue
                stat=(base/'stat').read_text().rsplit(')',1)[1].split()
                seen[(pid,stat[19])]=(int(stat[11])+int(stat[12]))/ticks
                roll=dict(line.split(':',1) for line in (base/'smaps_rollup').read_text().splitlines() if ':' in line)
                rss+=int(status.get('VmRSS','0 kB').split()[0])*1024
                pss+=int(roll.get('Pss','0 kB').split()[0])*1024
                private+=sum(int(roll.get(k,'0 kB').split()[0])*1024 for k in ('Private_Clean','Private_Dirty'))
                count+=1
            except (FileNotFoundError,ProcessLookupError,PermissionError):pass
        writer.writerow([datetime.now(timezone.utc).isoformat(),rss,private,sum(seen.values()),pss,count])
        f.flush();time.sleep(1)
