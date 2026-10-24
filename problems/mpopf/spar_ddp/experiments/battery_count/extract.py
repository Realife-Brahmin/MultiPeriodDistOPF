import csv,json,math,pathlib,re,sys
p=pathlib.Path(sys.argv[1]); method=sys.argv[2];text=(p/'solve.log').read_text()
kv=lambda s:dict(re.findall(r'(\w+)=([^\s]+)',s))
def lines(tag):return [kv(x) for x in text.splitlines() if x.startswith(tag+' ')]
def save(name,rows):
 assert rows,name
 with (p/name).open('w',newline='') as f:
  keys=list(dict.fromkeys(k for r in rows for k in r));w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
if method=='ddp':
 passes=lines('FILTERDDP_PROCS_PASS');original=lines('FILTERDDP_ORIGINAL_PRIMAL');timing=lines('FILTERDDP_ITER_TIMING')
 assert len(passes)==len(original)==len(timing)>0
 total=float(re.search(r'solve complete: ([\d.]+) s',text)[1]);setup=total-sum(float(r['total_s']) for r in timing)
 assert setup>=-.005,setup
 elapsed=0;trace=[]
 for a,b,t in zip(passes,original,timing):
  assert (a['iteration'],a['barrier_iteration'])==(b['iteration'],b['barrier_iteration'])==(t['iteration'],t['barrier_iteration'])
  assert all(math.isfinite(float(b[k])) for k in ('original_max','objective'))
  trace.append(dict(iteration=a['iteration'],barrier_iteration=a['barrier_iteration'],seconds=max(0,setup)+elapsed+float(t['backward_s']),objective=b['objective'],original=b['original_max'],dual=a['dual_inf'],complementarity=a['cs_inf']))
  elapsed+=float(t['total_s'])
 save('passes.csv',passes);save('iteration_timing.csv',timing);save('quality_trace.csv',trace)
 save('directions.csv',lines('FILTERDDP_PROCS_DIRECTIONS'));save('trials.csv',lines('FILTERDDP_PROCS_TRIAL'))
 summary={'full_optimizer_s':total,'timer_remainder_s':setup,'phase_A_s':sum(float(a['step1_s']) for a in passes),'scan_wall_s':sum(float(a['rounds_s']) for a in passes),'phase_C_s':sum(float(a['step3_s']) for a in passes),'memory':lines('FILTERDDP_PROCS_MEMORY')}
else:
 rows=lines('IPOPT_ITER');assert rows
 orig={}
 for line in (p/'ipopt_internal.log').read_text().splitlines():
  m=re.match(r'^\s*(\d+)r?\s+([-0-9.]+e[-+]\d+)\s+([0-9.]+e[-+]\d+)\s',line)
  if m:orig[int(m[1])]=m[3]
 trace=[dict(iteration=r['iteration'],seconds=r['elapsed_s'],objective=r['objective'],original=orig[int(r['iteration'])],dual=r['inf_du'],mu=r['mu']) for r in rows]
 save('quality_trace.csv',trace)
 summary={'central':lines('CENTRAL_IPOPT'),'memory':lines('CENTRAL_IPOPT_MEMORY'),'preprocessing':lines('CENTRAL_IPOPT_PREPROCESS_TIMING')}
(p/'timing_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
