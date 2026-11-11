"""Phase-bounded resident memory from preserved one-second node samples."""
import csv
from datetime import datetime
import json
from pathlib import Path
import tomllib

base=Path(__file__).parent/'problems/mpopf/tadmm/review_pilot/results/historical_anchors_nov09'
def rows(p):
    with p.open() as f:return list(csv.DictReader(f))
def stamp(s):return datetime.fromisoformat(s.replace('Z','+00:00')).timestamp()
result=[]
for c in sorted(base.glob('tadmm_review_*/results/run_*')):
    m=tomllib.loads((c/'sweep.toml').read_text())
    phases=[]
    events=rows(c/'phase_events.csv')
    for start in events:
        if not start['event'].endswith('_measured_begin'):continue
        phase=start['event'].removesuffix('_measured_begin')
        end=next(e for e in events if e['event']==phase+'_measured_end')
        nodes=[]
        for p in sorted(c.glob('process_samples.*.csv')):
            active=[r for r in rows(p) if stamp(start['utc'])<=stamp(r['utc'])<=stamp(end['utc'])]
            if not active:continue
            nodes.append(dict(node=p.name,samples=len(active),
                peak_pss_gib=max(int(r['pss_bytes']) for r in active)/2**30,
                peak_rss_gib=max(int(r['working_set_bytes']) for r in active)/2**30,
                peak_private_gib=max(int(r['private_bytes']) for r in active)/2**30,
                max_compute_processes=max(int(r['process_count']) for r in active)))
        phases.append(dict(phase=phase,nodes=nodes,
            sum_node_peak_pss_gib_upper_bound=sum(n['peak_pss_gib'] for n in nodes),
            sum_node_peak_rss_gib_upper_bound=sum(n['peak_rss_gib'] for n in nodes)))
    result.append(dict(case=c.name,system=m['system'],T=m['T'],solver=m['linear_solver'],
        allocated_nodes=m['nodes'],configured_memory_gb=m['configured_memory_gb'],phases=phases))
payload=dict(note='PSS apportions shared pages; RSS can double count them. Peaks sampled once per second during measured solves include registered coordinator/worker runtime. Multi-node sums are upper bounds on synchronized aggregate peaks. No swap or old-PC paging measurements are available in these sample files.',cases=result)
(base/'memory-peaks.json').write_text(json.dumps(payload,indent=2)+'\n')
for c in result:
    print(c['system'],c['T'],c['solver'])
    for p in c['phases']:
        print(p['phase'], 'PSS_GiB',round(p['sum_node_peak_pss_gib_upper_bound'],3),
              'RSS_GiB',round(p['sum_node_peak_rss_gib_upper_bound'],3),
              'nodes',[(n['node'],round(n['peak_pss_gib'],3),n['max_compute_processes']) for n in p['nodes']])
