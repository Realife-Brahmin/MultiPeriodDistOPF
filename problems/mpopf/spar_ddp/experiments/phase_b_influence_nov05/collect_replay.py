"""Strict posthoc collector for frozen influence diagnostics, never a benchmark."""
import argparse,csv,hashlib,json,math,re,tomllib
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
root,out=a.root,a.out
assert not out.exists();out.mkdir(parents=True)
r=root/'results'
def read(name):
 with (r/(name+'.csv')).open(newline='') as f:return list(csv.DictReader(f))
def write(name,rows):
 if not rows:return
 with (out/(name+'.csv')).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for x in rows for k in x)));w.writeheader();w.writerows(rows)
def checks(path):
 lines=path.read_text().splitlines();assert lines and all(x.endswith(': OK') for x in lines),path
for name in ('source-check.log','source-check-after.log','capture-runtime-input-check.log','capture-runtime-input-check-after.log'):checks(root/name)
for line in (root/'source.sha256').read_text().splitlines():
 digest,name=line.split(None,1);path=Path(name.strip())
 if path.is_absolute():
  # Runtime pin files use root-relative paths in this harness.
  raise AssertionError('Unexpected absolute replay source path')
 assert hashlib.sha256((root/path).read_bytes()).hexdigest()==digest,path
assert 'PRODUCTION_RUNTIME_UNIT_PASS' in (root/'runtime-unit.log').read_text()
assert 'INFLUENCE_SYNTHETIC_PASS' in (root/'unit.log').read_text()
assert (root/'replay.finished').exists()
account=(root/'accounting.txt').read_text()
assert re.search(r'(?m)^\d+\|COMPLETED\|[^|]+\|0:0\|',account)
assert 'Gold 6438Y+' in (root/'hardware.txt').read_text()
assert not (root/'stderr.log').read_text().strip()
log=(root/'replay.log').read_text()
assert 'INFLUENCE_AUDIT_COMPLETE' in log
meta=tomllib.loads((r/'audit_metadata.toml').read_text());assert meta['completed'] and meta['baseline_gate'] and meta['blas_threads']==1
assert not meta['timing_is_benchmark'] and meta['interventions_planned']==24
base=tomllib.loads((r/'baseline_gate.toml').read_text());assert base['pass'] and base['states']==8 and base['stages_per_state']==24
assert base['unit_scaled_gate']==1e-7 and base['relative_gate']==1e-5 and base['relative_floor']==1e-10
states={f'{name}_{when}' for name in ('early','middle','late','barrier') for when in ('before','after')}
status=read('baseline_status');assert len(status)==8 and {x['state'] for x in status}==states
assert all(x['production_pass']=='true' and x['negative_rejected']=='true' for x in status)
baseline=read('baseline');assert len(baseline)==2536,len(baseline)
for x in baseline:
 assert math.isfinite(float(x['unit_scaled'])) and float(x['unit_scaled'])<1e-7
 assert math.isfinite(float(x['relative_floor1e10'])) and float(x['relative_floor1e10'])<1e-5
 assert x['unchanged_error_gate']=='true'
negative=read('negative');assert len(negative)==16 and all(x['unchanged_error_gate']=='false' for x in negative)
failures=read('failures')
# Failures remain reportable adverse outcomes; never erase or promote them.
# Initial complete scientific influence collection requires all planned arms,
# or explicit classification in a later failure assessment, not this happy path.
assert not failures,failures
modes=('H','BL','BR','h','Hdiag','BLoffdiag')
contexts={f'{pair}_{mode}' for pair in ('early','middle','late','barrier') for mode in modes}
dep=read('dependencies');assert len(dep)==24*68
for context in contexts:
 rows=[x for x in dep if x['context']==context];assert len(rows)==68
 hit=[x for x in rows if x['injected']=='true'];assert len(hit)==1
 assert (hit[0]['depth'],hit[0]['left'],hit[0]['right'])==('2','5','7')
 assert any(x['right']=='5' and x['right_dependency']=='true' for x in rows)
 assert {int(x['left']) for x in rows if x['affected']=='true'}==set(range(1,6))
tables={n:read(n) for n in ('merge_errors','suffix_errors','directions','trajectory','actions','changes','condition','ranks','spectra','costs','independent','input_errors','quantiles','blocks','outliers')}
assert len(tables['directions'])==24*24*7*4
assert len(tables['trajectory'])==24*24*6
assert len(tables['actions'])==24*24*7*2
assert len(tables['suffix_errors'])==24*25*2
assert len(tables['changes'])==4*3*4
assert len(tables['condition'])==4*3*2
assert len(tables['ranks'])==4*3*3*2*4
assert {x['kind'] for x in tables['spectra']}=={'delta','current'}
for name in ('directions','trajectory','actions','suffix_errors','merge_errors'):
 assert all(math.isfinite(float(x['unit_scaled'])) and math.isfinite(float(x['relative_floor1e10'])) for x in tables[name]),name
summary=[];first=[]
for context in sorted(contexts):
 result={'context':context,'interpretation':'mixed_input_sensitivity_not_deployable_reuse'}
 for name in ('merge_errors','suffix_errors','actions','directions','trajectory'):
  rows=[x for x in tables[name] if x['context']==context or x['context'].startswith(context+'_')]
  assert rows,(context,name)
  result[name+'_max_unit_scaled']=max(float(x['unit_scaled']) for x in rows)
  result[name+'_max_relative_floor1e10']=max(float(x['relative_floor1e10']) for x in rows)
  bad=[x for x in rows if x['unchanged_error_gate']=='false']
  result[name+'_threshold_exceedances']=len(bad)
  if bad:
   # CSV merge order is exact global chronological scan order. Stage outputs
   # are stage-ordered; no artificial scalar amplification factor is inferred.
   x=bad[0];first.append(dict(intervention=context,table=name,first_context=x['context'],stage=x['stage'],field=x['field'],unit_scaled=x['unit_scaled'],relative_floor1e10=x['relative_floor1e10']))
 result['all_recovered_direction_probes_within_original_gates']=result['directions_threshold_exceedances']==result['trajectory_threshold_exceedances']==0
 summary.append(result)
write('interventions',summary);write('first_amplification',first)
# Independent consumer/reference errors are explicit separate outcomes.
ind=[]
for x in tables['independent']:
 if '_before_' in x['context'] or '_after_' in x['context']:ind.append(x)
write('independent_reference',ind)
write('largest_changes',sorted(tables['changes'],key=lambda x:float(x['unit_scaled']),reverse=True))
write('delta_ranks',[x for x in tables['ranks'] if x['kind']=='delta'])
operations={}
for x in tables['costs']:
 key=x['operation'];v=operations.setdefault(key,dict(operation=key,calls=0,diagnostic_s=0.0,allocation_bytes=0,gc_s=0.0))
 v['calls']+=1;v['diagnostic_s']+=float(x['seconds']);v['allocation_bytes']+=int(x['allocation_bytes']);v['gc_s']+=float(x['gc_s'])
write('diagnostic_costs',list(operations.values()))
validation=dict(validated=True,baseline_probe_rows=len(baseline),negative_control_rows=len(negative),interventions=24,full_stage_count=24,original_scaled_gate=1e-7,original_relative_gate=1e-5,relative_floor=1e-10,scientific_failures=failures,timing_is_benchmark=False,cost_limitation='retained-input copies and summary/outlier processing not individually timed; whole job includes them',memory_limitation='process RSS includes captures/references/JIT/validation; no private-memory savings claim')
(out/'validation.json').write_text(json.dumps(validation,indent=2))
print(json.dumps(validation,indent=2))
