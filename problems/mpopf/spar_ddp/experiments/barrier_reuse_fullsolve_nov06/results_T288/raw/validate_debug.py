"""Posthoc validator for actual complete diagnostic output; not a benchmark."""
import argparse,collections,csv,hashlib,json,math,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('collection',type=Path);a=p.parse_args()
r=a.collection;raw=r/'raw';prov=r/'provenance';job=prov/'jobs/31612890'
def read_csv(path):
 with path.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
def checks_ok(path):
 s=path.read_text(encoding='utf-8');assert s.strip() and all(x.endswith(': OK') for x in s.splitlines()),path
for f in ('inherited-source-check.log','inherited-input-check.log','source-check-after-debug.log','harness-check-after-debug.log','runtime-check-after-debug.log'):checks_ok(prov/f)
for f in ('source-check.log','harness-check.log'):checks_ok(job/f)
for f,marker in (('unit.log','BR_UNIT_PASS'),('worker-unit.log','BR_WORKER_UNIT_PASS'),('shadow-unit.log','BR_SHADOW_UNIT_PASS'),('parse.log','BR_PARSE_PASS')):assert marker in (job/f).read_text(encoding='utf-8')
s=(raw/'solve.log').read_text(encoding='utf-8');warm=(raw/'warmup-debug.log').read_text(encoding='utf-8')
assert (raw/'COMPLETE').is_file()
assert not any(x in s+warm for x in ('BR_DEBUG_FAIL','AssertionError','ERROR: LoadError','FILTERDDP_PROCS_FAIL','WARMUP_STREAM_MISMATCH'))
assert re.search(r'solve complete: [\d.]+ s, iterations=14, status=0\b',s)
final=re.search(r'final residuals: primal=(\S+) dual=(\S+) complementarity=(\S+)',s);assert final
assert all(math.isfinite(float(x)) and 0<=float(x)<=1e-5 for x in final.groups())
physical=float(re.search(r'original_max=(\S+)',(raw/'physical.log').read_text())[1]);assert 0<=physical<=1e-5
assert 'BOUND_CHECK ' in s and 'worst_violation=0.000e+00' in s
dss=(raw/'opendss.log').read_text();assert 'below limit: 0' in dss and 'above limit: 0' in dss
checks=[dict(label=m[0],absolute=float(m[1]),scaled=float(m[2]),relative=float(m[3]),maxabs=float(m[4]),exact=m[5]=='true') for m in re.findall(r'BR_DEBUG_CHECK label=(\S+) absolute=(\S+) scaled=(\S+) relative=(\S+) maxabs=(\S+) exact=(\S+)',s)]
assert checks and all(all(math.isfinite(x[k]) for k in ('absolute','scaled','relative','maxabs')) and x['scaled']<1e-7 and x['relative']<1e-5 for x in checks)
cache=[dict(re.findall(r'(\w+)=(\S+)',line)) for line in re.findall(r'^BR_PASS solve=1 .*',s,re.M)]
passes=read_csv(raw/'passes.csv');timing=read_csv(raw/'iteration_timing.csv');quality=read_csv(raw/'quality_trace.csv')
assert len(cache)==len(passes)==len(timing)==len(quality)==23
assert [(x['iteration'],x['barrier_iteration']) for x in cache]==[(x['iteration'],x['barrier_iteration']) for x in quality]
assert all(x['debug']=='true' and x['valid']=='true' and int(x['entries'])==52 and int(x['bytes'])==1298419200 and int(x['fallback'])==0 for x in cache)
assert sum(x['hit']=='true' for x in cache)==8
for i,x in enumerate(cache):
 eligible=i>0 and x['iteration']==cache[i-1]['iteration'] and int(x['barrier_iteration'])==int(cache[i-1]['barrier_iteration'])+1 and float(passes[i]['mu'])<float(passes[i-1]['mu'])
 assert (x['requested']=='true')==eligible and (x['hit']=='true')==eligible
 assert int(x['reused'])==(52 if eligible else 0) and int(x['populated'])==(0 if eligible else 52)
counts=collections.Counter(x['label'] for x in checks)
ops=[];h=1
while h<=24:
 ops.append([(t,i+h) for i in range(1,26-h,2*h) for t in range(i,(i if h<2 else i+h-1)+1,2)]);h*=2
ops.append([(i,i+1) for i in range(2,25,2)]);assert sum(map(len,ops))==52
val={25:True}
for depth,pairs in enumerate(ops,1):
 for t,mid in pairs:
  closed=val.get(mid,False)
  for f in (('V','v') if closed else ('Maa','Mab','Mbb','ma','mb')):assert counts[f'shadow_{depth}_{t}_{f}']==23
  for f in (('ma',) if closed else ('ma','mb')):assert counts[f'merge_{depth}_{t}_{f}']==23
  val[t]=closed
for t in range(1,25):
 for f in ('ma','mb'):assert counts[f'leaf_{t}_{f}']==23
 for f in ('xE','Q'+chr(251),'B_active','XB'):assert counts[f'recovery_{t}_{f}']==23
 for f in (chr(945),chr(968),chr(967)+'l',chr(967)+'u',chr(931)+'_L',chr(931)+'_U'):assert counts[f'rule_{t}_{f}']==23
 assert counts[f'feedback_{t}']==23
 actions=collections.Counter()
 for label,count in counts.items():
  m=re.fullmatch(r'(control|equality|lower|upper)_action_(\d+)_(\d+)_(\S+)',label)
  if m and int(m[2])==t:actions[(m[1],int(m[3]),float(m[4]))]+=count
 assert len(actions)==36 and all(x==23 for x in actions.values())
 assert {x[2] for x in actions}=={1e-4,1.,10.} and {x[1] for x in actions}=={1,2,3}
 for f in ('control','equality','lower','upper'):assert counts[f'actual_{t}_{f}']==14
 if t<24:assert counts[f'actual_forward_{t}']==14
neg=re.findall(r'BR_NEGATIVE stage=(\d+) kind=(\S+) rejected=(\S+)',s)
assert len(neg)==168 and all(x[2]=='true' for x in neg)
assert all(len({x[1] for x in neg if int(x[0])==t})==7 for t in range(1,25))
stale=re.findall(r'BR_STALE_VECTOR depth=(\d+) stage=(\d+) rejected=(\S+) scaled=(\S+) relative=(\S+)',s)
assert len(stale)==416 and all(x[2]=='true' and (float(x[3])>=1e-7 or float(x[4])>=1e-5) for x in stale)
groups=[set(map(int,x.split(','))) for x in (job/'24_debug_1/cpu_groups.csv').read_text().splitlines()]
assert len(groups)==25 and len(groups[0])==1 and all(len(x)==2 for x in groups[1:])
assert len({x for g in groups for x in g})==49
topology={}
for line in (job/'24_debug_1/cpu_topology.csv').read_text().splitlines():
 if not line.startswith('#'):
  cpu,core,socket=map(int,line.split(','));topology[cpu]=(socket,core)
assert len({topology[c] for g in groups for c in g})==49
assert all(len({topology[c][0] for c in g})==1 for g in groups)
hw=(job/'hardware.txt').read_text();assert '6760P' in hw and 'NumCPUs=128' in hw
for tag in ('after-setup','after-solve'):
 audit=(raw/f'worker-audit-{tag}.txt').read_text().splitlines();assert len(audit)==24
 for i,line in enumerate(audit):
  assert 'libopenblas' in line and int(re.search(r'blas\s*=\s*(\d+)',line)[1])==2 and int(re.search(r'blas_default\s*=\s*(\d+)',line)[1])==2
  masks=re.findall(r'allowed\s*=\s*"([^"]+)"',line);assert masks
  for mask in masks:
   cpus=set()
   for span in mask.split(','):
    ends=list(map(int,span.split('-')));cpus.update(range(ends[0],ends[-1]+1))
   assert cpus==groups[i+1]
assert 'coordinator_threads=1' in s and 'BR_RETIRE bytes=' in s and 'collected=true' in s
for name,data in (('debug_checks.csv',checks),('cache_passes.csv',cache)):
 with (r/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
result=dict(validated=True,scope='diagnostic full T24 solve only; NOT clean timing comparison',stages=24,passes=23,accepted_steps=14,hits=8,fresh_passes=15,merges_per_pass=52,full_shadow_nodes=1196,debug_checks=len(checks),actual_negatives=len(neg),stale_gradient_rejections=len(stale),max_scaled=max(x['scaled'] for x in checks),max_relative=max(x['relative'] for x in checks),max_absolute=max(x['absolute'] for x in checks),all_checks_bitexact=all(x['exact'] for x in checks),physical_original_residual=physical,source_harness_runtime_checks=True,cache_bytes=1298419200,archive_sha256=hashlib.sha256((r/'evidence.tgz').read_bytes()).hexdigest())
(r/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
