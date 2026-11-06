"""Audit legacy target-clock allocation without altering raw evidence or gates."""
import argparse,csv,json,statistics,math
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('case',type=Path);a=p.parse_args();r=a.case
with (r/'analysis/per_run.csv').open(encoding='utf-8',newline='') as f:per=list(csv.DictReader(f))
rows=[]
for x in per:
 T=int(x['T']);folder=r/f'raw/results/10k_1020_T{T}_ddp_r{x["arm"]}_{x["repeat"]}'
 summary=json.loads((folder/'timing_summary.json').read_text());legacy=float(x['quality_s']);rem=float(summary['timer_remainder_s']);prefix=legacy-max(0.,rem)
 assert all(math.isfinite(z) for z in (legacy,rem,prefix)) and prefix>0 and rem>=-.005
 rows.append(dict(arm=x['arm'],repeat=int(x['repeat']),T=T,legacy_frontloaded_target_s=legacy,entire_solve_unallocated_wall_s=rem,logged_iteration_prefix_s=prefix,target_elapsed_lower_bound_s=prefix,target_elapsed_upper_bound_s=legacy+.001))
med={}
for arm in ('base','cache'):
 armrows=[x for x in rows if x['arm']==arm]
 med[arm]={k:statistics.median(x[k] for x in armrows) for k in rows[0] if k not in ('arm','repeat','T')}
b,c=med['base'],med['cache']
result=dict(target_direct_wall_clock_recorded=False,legacy_estimator='entire full optimizer minus summed iteration durations is allocated before first pass; includes post-loop exit/retirement',logged_prefix_reduction_percent=100*(1-c['logged_iteration_prefix_s']/b['logged_iteration_prefix_s']),wall_target_reduction_lower_bound_percent=100*(1-c['target_elapsed_upper_bound_s']/b['target_elapsed_lower_bound_s']),wall_target_reduction_upper_bound_percent=100*(1-c['target_elapsed_lower_bound_s']/b['target_elapsed_upper_bound_s']),medians=med,full_optimizer_results_unchanged=True,raw_evidence_unchanged=True,scientific_thresholds_unchanged=True)
with (r/'analysis/quality_clock_audit.csv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(r/'analysis/quality_clock_audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
