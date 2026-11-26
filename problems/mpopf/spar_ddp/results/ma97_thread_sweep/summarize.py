"""Generate comparisons from preserved text evidence; never edit raw logs."""
from pathlib import Path
import re,json,csv
root=Path(__file__).resolve().parent
kv=lambda line:dict(re.findall(r'(\w+)=([^\s]+)',line))
rows=[];baseline=None;reference=None
for job,threads in [(31251874,1),(31251875,4),(31251876,16),(31251877,64)]:
 p=next((root/'evidence_v2/logs').glob(f'{job}_*'))
 s=(p/'solve.log').read_text(); raw=(p/'ipopt_internal.log').read_text()
 assert (p/'COMPLETE').exists()
 assert json.loads((p/'logging_check.json').read_text())['logging_validated']
 assert not re.search(r'ERROR:|UNHANDLED TASK|ProcessExitedException',s+(p/'stderr.log').read_text())
 runtime=kv(next(x for x in s.splitlines() if x.startswith('MA97_THREAD_RUNTIME ')))
 assert int(runtime['omp_max_threads'])==threads and runtime['blas_threads']=='1'
 assert '6438Y+' in (p/'metadata.txt').read_text()
 assert 'below limit: 0' in (p/'opendss.log').read_text() and 'above limit: 0' in (p/'opendss.log').read_text()
 c=kv(next(x for x in s.splitlines() if x.startswith('CENTRAL_IPOPT system=')))
 assert c['status']=='LOCALLY_SOLVED' and int(c['iterations'])==18
 records=[kv(x) for x in s.splitlines() if x.startswith('IPOPT_ITER ')]
 trajectory=[tuple(x[k] for k in ('iteration','objective','inf_pr','inf_du')) for x in records]
 if baseline is None: baseline=trajectory;reference=float(c['objective'])
 assert trajectory==baseline, 'Numerical trajectory changed'
 original={}
 for line in raw.splitlines():
  m=re.match(r'^\s*(\d+)r?\s+([-0-9.]+e[-+]\d+)\s+([0-9.]+e[-+]\d+)\s',line)
  if m:original[int(m[1])]=float(m[3])
 def quality(gap):
  return next(float(x['elapsed_s']) for x in records if original[int(x['iteration'])]<=1e-5 and abs(float(x['objective'])-reference)/abs(reference)<=gap)
 def wall(category):
  m=re.findall(re.escape(category)+r'\.+:\s+[0-9.]+\s+\(sys:\s+[0-9.]+\s+wall:\s+([0-9.]+)\)',raw)
  return float(m[-1]) if m else None
 last=records[-1]
 maxrss=int(kv(next(x for x in s.splitlines() if x.startswith('MA97_COMPLETE_WALL ')))['maxrss_bytes'])
 row=dict(job=job,omp_threads=threads,blas_threads=1,solve_s=float(c['solve_time_s']),build_s=float(c['build_s']),runner_wall_s=float(c['wall_s']),iterations=18,objective=float(c['objective']),original_primal=original[18],callback_primal=float(last['inf_pr']),callback_dual=float(last['inf_du']),quality_05percent_s=quality(.005),quality_01percent_s=quality(.001),process_maxrss_gib=maxrss/2**30,trajectory_identical=True)
 for category in ['OverallAlgorithm','ComputeSearchDirection','PDSystemSolverTotal','StdAugSystemSolverMultiSolve','LinearSystemSymbolicFactorization','LinearSystemFactorization','LinearSystemBackSolve','Function Evaluations']:
  row[category+'_wall_s']=wall(category)
 comp=re.findall(r'Complementarity\.+:\s+([-+0-9.e]+)\s+([-+0-9.e]+)',raw)
 row['complementarity']=float(comp[-1][1])
 rows.append(row)
(root/'comparison.json').write_text(json.dumps(rows,indent=2)+'\n')
with (root/'comparison.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps(rows,indent=2))
