"""Strict exact-cache pilot collection. Timing wins never override science gates."""
import argparse,csv,hashlib,json,math,re,statistics,tomllib
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
r=a.root;d=r/'results';o=a.out;assert not o.exists();o.mkdir(parents=True)
def rows(name):
 with (d/(name+'.csv')).open(newline='') as f:return list(csv.DictReader(f))
def write(name,data):
 if data:
  with (o/(name+'.csv')).open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for x in data for k in x)));w.writeheader();w.writerows(data)
assert (r/'pilot.finished').exists()
assert 'CACHE_UNIT_PASS' in (r/'unit.log').read_text()
log=(r/'pilot.log').read_text();assert 'EXACT_BARRIER_CACHE_PILOT_COMPLETE' in log
assert not (r/'stderr.log').read_text().strip()
assert re.search(r'(?m)^\d+\|COMPLETED\|[^|]+\|0:0\|',(r/'accounting.txt').read_text())
hardware=(r/'hardware.txt').read_text();assert all(x in hardware for x in ('Gold 6438Y+','NumCPUs=64','mem=300G','Requeue=0'))
for name in ('source-check.log','source-check-after.log','runtime-input-check.log','runtime-input-check-after.log'):
 x=(r/name).read_text().splitlines();assert x and all(t.endswith(': OK') for t in x),name
for line in (r/'source.sha256').read_text().splitlines():
 h,name=line.split(None,1);assert hashlib.sha256((r/name.strip()).read_bytes()).hexdigest()==h,name
meta=tomllib.loads((d/'validation.toml').read_text());assert meta['complete'] and meta['global_leaf_matrix_invariance'] and meta['global_68_merge_invariance'] and meta['signature_invariance']
assert meta['measured_pairs']==5 and meta['scaled_gate']==1e-7 and meta['relative_gate']==1e-5 and meta['relative_floor']==1e-10
assert meta['process_workers']==24 and meta['blas_threads']==1 and meta['warmup_compilation_cleared']
inv=rows('invariance');assert len(inv)==521,len(inv)
assert sum(x['kind']=='leaf' for x in inv)==72 and sum(x['kind']=='terminal' for x in inv)==1
assert len({(x['depth'],x['left'],x['right']) for x in inv if x['kind']=='merge'})==68
assert all(x['exact_bits']=='true' and float(x['absolute_fro'])==0 for x in inv)
e=rows('errors');negative=[x for x in e if x['scope']=='stale_gradient_negative'];positive=[x for x in e if x['scope']!='stale_gradient_negative']
assert negative and any(x['pass']=='false' for x in negative)
assert all(x['pass']=='true' and float(x['unit_scaled'])<1e-7 and float(x['relative_floor1e10'])<1e-5 for x in positive)
for arm in ('baseline','candidate'):
 for rep in range(1,6):
  selected=[x for x in e if x['scope']=='measured_'+arm and int(x['repeat'])==rep]
  assert len(selected)==24*49
  assert {int(x['stage']) for x in selected}==set(range(1,25))
  for stage in range(1,25):
   q=[x for x in selected if int(x['stage'])==stage]
   assert len(q)==49 and len({(x['probe'],x['field']) for x in q})==49
layout=meta.get('layout_version',1)
controls=rows('negative_controls');assert len(controls)==(7 if layout==2 else 6)
if layout==2:
 assert meta['population_full_output_gate']
 population=[x for x in e if x['scope']=='cache_population_preflight'];assert len(population)==24*49
assert all(x['rejected_or_invalidated']=='true' for x in controls)
assert all(int(x['hits'])==0 and x['correct_fresh_outputs']=='true' for x in controls if x['control']!='stale_gradient_consumer')
s=tomllib.loads((d/'storage.toml').read_text());assert s['merge_entries']==68 and s['leaf_entries']==(0 if layout==2 else 24)
if layout==2:assert s['closed_values']==24
assert s['numeric_bytes']==(204 if layout==2 else 408)*1020*1020*8
assert (s['skipped_cholesky_per_hit'],s['skipped_matrix_triangular_solves_per_hit'],s['skipped_matrix_products_per_hit'])==(68,112,156)
aff=rows('affinity');assert len(aff)==24 and len({x['allowed'] for x in aff})==24
assert all(x['blas_threads']=='1' and ',' not in x['allowed'] and '-' not in x['allowed'] for x in aff)
groups=[list(map(int,x.split(','))) for x in (r/'cpu_groups.csv').read_text().splitlines()]
assert len(groups)==25 and all(len(x)==1 for x in groups)
assert {int(x['allowed']) for x in aff}=={x[0] for x in groups[1:]}
top={}
for line in (r/'cpu_topology.csv').read_text().splitlines():
 if line.startswith('#'):continue
 cpu,core,socket=map(int,line.split(','));top[cpu]=(socket,core)
assert len({top[g[0]] for g in groups})==25
workers=tomllib.loads((d/'workers.toml').read_text())['workers'];assert len(workers)==24
assert all(x['blas']==x['threads']==1 and 'libopenblas' in x['backend'] for x in workers)
t=rows('totals');ph=rows('timings');assert len(t)==10 and len(ph)==20
expected=[(str(rep),arm) for rep in range(1,6) for arm in (('baseline','candidate') if rep%2 else ('candidate','baseline'))]
assert [(x['repeat'],x['arm']) for x in t]==expected
assert all(float(x['worker_compile_ns'])==float(x['worker_recompile_ns'])==0 for x in t)
assert all(float(x['coordinator_compile_s'])==float(x['coordinator_recompile_s'])==0 for x in ph)
for x in ph:
 assert int(x['fallbacks'])==0
 hit=x['arm']=='candidate' and x['phase']=='hit_or_fresh_second'
 pop=x['arm']=='candidate' and x['phase']=='population_or_fresh_first'
 assert int(x['cache_hits'])==(68 if hit else 0)
 assert int(x['populated'])==(68 if pop else 0)
 assert all(math.isfinite(float(x[k])) and float(x[k])>=0 for k in ('wall_s','detection_s','leaf_s','scan_s','consumer_s','forward_s'))
summary=[]
for arm in ('baseline','candidate'):
 z=[x for x in t if x['arm']==arm];res=dict(arm=arm,repeats=5)
 for key in ('charged_s','first_pass_s','second_pass_s','retirement_gc_s','worker_alloc_bytes','worker_gc_s'):
  res[key]=statistics.median(float(x[key]) for x in z)
 res['charged_min_s']=min(float(x['charged_s']) for x in z);res['charged_max_s']=max(float(x['charged_s']) for x in z)
 summary.append(res)
write('summary',summary);write('paired',[dict(repeat=rep,baseline_s=float(next(x['charged_s'] for x in t if x['repeat']==str(rep) and x['arm']=='baseline')),candidate_s=float(next(x['charged_s'] for x in t if x['repeat']==str(rep) and x['arm']=='candidate'))) for rep in range(1,6)])
phase_summary=[]
for arm in ('baseline','candidate'):
 for phase in ('population_or_fresh_first','hit_or_fresh_second'):
  z=[x for x in ph if x['arm']==arm and x['phase']==phase]
  phase_summary.append(dict(arm=arm,phase=phase,**{k:statistics.median(float(x[k]) for x in z) for k in ('wall_s','detection_s','leaf_s','scan_s','consumer_s','forward_s')}))
write('phases',phase_summary)
write('negative_controls',controls)
worst=[]
for scope in sorted({x['scope'] for x in e}):
 z=[x for x in e if x['scope']==scope]
 for key in ('absolute_fro','unit_scaled','relative_floor1e10'):
  x=max(z,key=lambda x:float(x[key]));worst.append(dict(scope=scope,metric=key,stage=x['stage'],probe=x['probe'],field=x['field'],value=x[key]))
write('worst_errors',worst)
validation=dict(validated=True,global_exact_matrix_rows=len(inv),clean_two_pass_workloads=10,measured_pairs=5,current_output_gates_unchanged=True,negative_control_count=len(controls),cache_layout_version=layout,population_full_output_gate=(layout==2),cache_numeric_bytes=s['numeric_bytes'],scope=meta['scope'],timing_includes_cache_population=True,worker_allocation_delta_includes_counter_rpc_overhead=True,worker_rss_is_shared_mapping_highwater_not_private_memory=True)
(o/'validation.json').write_text(json.dumps(validation,indent=2));print(json.dumps(summary,indent=2));print(json.dumps(validation,indent=2))
