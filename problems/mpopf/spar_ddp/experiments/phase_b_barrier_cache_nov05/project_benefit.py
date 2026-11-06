"""Explicit cost-balance scenarios, not new optimizer measurements."""
import argparse,csv,statistics,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('analysis',type=Path);p.add_argument('frequency',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
def read(path):
 with path.open() as f:return list(csv.DictReader(f))
ph=read(a.analysis/'phases.csv');freq=read(a.frequency/'frequency.csv')
def row(arm,phase):return next(x for x in ph if x['arm']==arm and x['phase']==phase)
bf=row('baseline','population_or_fresh_first');bh=row('baseline','hit_or_fresh_second');cf=row('candidate','population_or_fresh_first');ch=row('candidate','hit_or_fresh_second')
result=[]
for f in freq:
 n=int(f['source_eligible']);fresh=int(f['other_fresh_passes']);available=float(f['eligible_scan_s']);whole=float(f['optimizer_s'])
 detect=fresh*float(cf['detection_s'])+n*float(ch['detection_s'])
 scan_ratio=float(ch['scan_s'])/float(bh['scan_s'])
 pop=fresh*(float(cf['scan_s'])-float(bf['scan_s']))
 gross=available*(1-scan_ratio)
 result.append(dict(repeat=f['repeat'],source_conditional_hits=n,fresh_population_passes=fresh,numerically_verified_representative_transitions=1,original_optimizer_s=whole,eligible_fresh_scan_s=available,zero_cost_eligible_scan_ceiling_fraction=available/whole,representative_hit_scan_ratio=scan_ratio,gross_kernel_saving_scenario_s=gross,population_scan_overhead_scenario_s=pop,exact_detection_cost_scenario_s=detect,net_scan_and_detection_saving_scenario_s=gross-pop-detect,assumptions='all8keysunchanged_and_Choleskycacheable_representativeperpasscosts_samehardwareclass;not_a_solve_measurement;excludes_additional_GC_and_consumer_effects'))
a.out.mkdir(exist_ok=False)
with (a.out/'projections.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=result[0]);w.writeheader();w.writerows(result)
print(json.dumps(result,indent=2))
