"""Chronology plus source eligibility; never infer exact numerical keys from logs."""
import argparse,csv,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
rows=[];totals=[]
for rep in (1,2):
 d=a.root/f'r{rep}'
 with (d/'iteration_timing.csv').open() as f:timing=list(csv.DictReader(f))
 with (d/'passes.csv').open() as f:passes=list(csv.DictReader(f))
 assert len(timing)==len(passes)==23
 settings=(d/'settings.txt').read_text();assert 'FILTERDDP_PROCS_WORKERS=24' in settings and 'OPENBLAS_NUM_THREADS=1' in settings
 assert not any(x.startswith('FILTERDDP_HYBRID') or x.startswith('FILTERDDP_STALE') for x in settings.splitlines())
 eligible=[]
 for index,(prev,cur) in enumerate(zip(timing,timing[1:]),start=1):
  source_eligible=(prev['outcome']=='barrier_update' and cur['iteration']==prev['iteration'] and int(cur['barrier_iteration'])==int(prev['barrier_iteration'])+1 and float(cur['mu'])<float(prev['mu']))
  if source_eligible:
   eligible.append(index)
   rows.append(dict(repeat=rep,previous_pass=index,current_pass=index+1,k=cur['iteration'],j_before=prev['barrier_iteration'],j_after=cur['barrier_iteration'],mu_before=prev['mu'],mu_after=cur['mu'],source_fixed_unshifted_path=True,accepted_step_between=False,hybrid_switch_configured=False,numerically_fingerprinted_in_full_capture=(cur['iteration']=='4' and cur['barrier_iteration']=='1'),classification='source_eligible_under_pinned_path_exact_numeric_key_only_known_for_captured_pair',fresh_scan_s=passes[index]['rounds_s']))
 summary=json.loads((d/'timing_summary.json').read_text())
 totals.append(dict(repeat=rep,total_passes=len(passes),source_eligible=len(eligible),other_fresh_passes=len(passes)-len(eligible),eligible_scan_s=sum(float(passes[i]['rounds_s']) for i in eligible),all_scan_s=sum(float(x['rounds_s']) for x in passes),optimizer_s=summary['full_optimizer_s'],exact_captured_transition_count=1))
a.out.mkdir(exist_ok=False)
for name,data in (('transitions',rows),('frequency',totals)):
 with (a.out/(name+'.csv')).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
print(json.dumps(totals,indent=2))
