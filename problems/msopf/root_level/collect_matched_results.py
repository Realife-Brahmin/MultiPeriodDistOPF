"""Validate checkpointed matched-study evidence and prepare publication data.
No simulation, no interpolation of missing rays, no modification of raw CSVs.
"""
import csv,json,math,hashlib,sys
from pathlib import Path
import numpy as np
from prediction_polytope import prediction_mesh,attach_prediction_meshes
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/matched_secondary_oct25'
CLUSTER=OUT/'cluster/results'
LABELS=['primary2','secondary2','primary3','secondary3','primary4','secondary4','sf4']
GOOD={'LOCALLY_SOLVED','OPTIMAL','ALMOST_LOCALLY_SOLVED'}
RESOLVED={'direct_checked','bracketed','coarse_bracket','oracle_bracket','oracle_coarse_bracket'}
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def complete(p):return (p/'metadata.txt').exists() and 'exit_code=0' in (p/'metadata.txt').read_text(encoding='utf-8')
def validate(p):
 good=[];failed=[]
 for filename in ('solves.csv','direct_replays.csv'):
  if not (p/filename).exists():continue
  for row in read(p/filename):
   if row['status'] not in GOOD:failed.append(row['status']);continue
   assert row['replay_ok']=='true',(p,filename,row)
   for key,limit in [('dP_kW',1e-3),('dQ_kvar',1e-3),('dV_pu',1e-7),('dB_kW',1e-3)]:assert math.isfinite(float(row[key])) and float(row[key])<limit,(p,key,row)
   for key in ('terminal_kWh','soc_violation_kWh','energy_balance_kWh','battery_power_violation_kW','pv_q_violation_kvar','original_violation'):
    if key in row:assert float(row[key])<.01,(p,key,row[key])
   good.append(row)
 return {'verified_solves':len(good),'unsuccessful_statuses':failed,'max_dP_kW':max([float(x['dP_kW']) for x in good] or [0]),'max_dV_pu':max([float(x['dV_pu']) for x in good] or [0])}

def validate_oracle(p):
 rows=read(p/'oracle_trials.csv');good=[]
 for row in rows:
  if row['replay_ok']!='true':continue
  for key,limit in [('dP_kW',1e-3),('dQ_kvar',1e-3),('dV_pu',1e-7),('dB_kW',1e-3),('terminal_kWh',.01),('soc_violation_kWh',.01),('energy_balance_kWh',.01),('battery_power_violation_kW',.01),('pv_q_violation_kvar',.01)]:
   assert math.isfinite(float(row[key])) and float(row[key])<limit,(p,key,row)
  if row['physical_feasible']=='true':
   assert float(row['vmin'])>=.95-1e-6 and float(row['vmax'])<=1.05+1e-6 and float(row['min_import_kW'])>=-.01
  if row['local_outside']=='true':assert row['status'] in GOOD and row['physical_feasible']=='false'
  good.append(row)
 return dict(verified_solves=len(good),unsuccessful_statuses=[r['status'] for r in rows if r['status'] not in GOOD],max_dP_kW=max([float(x['dP_kW']) for x in good] or [0]),max_dV_pu=max([float(x['dV_pu']) for x in good] or [0]))
def arrays(folder,n):
 d=read(next(folder.glob('sensitivity_*.csv')))
 P=np.array([[float(x['import_kW']) for x in d if x['source']==f'subs{i}'] for i in range(1,n+1)])
 K=np.array([[[float(x[f'K_subs{j}']) for j in range(2,n+1)] for x in d if x['source']==f'subs{i}'] for i in range(1,n+1)])
 assert P.shape==(n,24) and K.shape==(n,24,n-1)
 return P,K
def prediction(P,K,case,u):
 a=K@u
 candidates=np.where(a.sum(axis=1)<0,P.sum(axis=1)/np.maximum(-a.sum(axis=1),1e-300),np.inf) if case=='cyclic' else np.min(np.where(a<0,P/np.maximum(-a,1e-300),np.inf),axis=1)
 return float(np.min(candidates)),int(np.argmin(candidates))
def directions(n):
 if n==2:return np.array([[-1.],[1.]])
 if n==3:return np.array([[math.cos(2*math.pi*i/32),math.sin(2*math.pi*i/32)] for i in range(32)])
 out=[]
 for i in range(78):
  z=1-2*(i+.5)/78;r=math.sqrt(1-z*z);phi=i*math.pi*(3-math.sqrt(5));out.append([r*math.cos(phi),r*math.sin(phi),z])
 return np.array(out)
def main(strict=False):
 results={};progress=[];proof={};profile_hashes={};publication_edges=[];contradictions=[]
 for label in LABELS:
  folder=CLUSTER/label;n=int(label[-1]);expected={2:2,3:32,4:78}[n];base=read(folder/'prep/baseline.csv') if (folder/'prep/baseline.csv').exists() else []
  ready=complete(folder/'prep');refready=complete(folder/'reference');edges={};written=0
  if ready:
   proof[str((folder/'prep').relative_to(OUT))]=validate(folder/'prep')
   profile_hashes[label]=hashlib.sha256((folder/'prep/profiles.csv').read_bytes()).hexdigest()
  for p in sorted(folder.glob('*'))+sorted((OUT/'cluster/oracle_results'/label).glob('*')):
   if not p.is_dir() or not (p/'edges.csv').exists():continue
   raw=read(p/'edges.csv');written+=len(raw)
   if not complete(p):continue
   isoracle=(p/'oracle_trials.csv').exists()
   proof[str(p.relative_to(OUT))]=validate_oracle(p) if isoracle else validate(p)
   for e in raw:
    assert int(e['T'])==24 and e['configuration']==label
    i=int(e['direction']);u=np.array([float(e[f'u{j}']) for j in range(1,n)])
    assert np.allclose(u,directions(n)[i-1],atol=1e-12)
    assert float(e['vmin'])>=.95-1e-6 and float(e['vmax'])<=1.05+1e-6 and float(e['min_import_kW'])>=-.01
    key=(e['case'],i)
    if isoracle and float(e['measured_deg'])>0:
     trials=read(p/'oracle_trials.csv')
     assert any(t['physical_feasible']=='true' and abs(float(t['radius_deg'])-float(e['measured_deg']))<1e-8 for t in trials),(p,'unverified lower')
     if e['method'] in RESOLVED:assert any(t['local_outside']=='true' and abs(float(t['radius_deg'])-float(e['upper_deg']))<1e-8 for t in trials),(p,'unverified upper')
    if key in edges:
     old=edges[key];oldwidth=float(old['upper_deg'])-float(old['measured_deg']);newwidth=float(e['upper_deg'])-float(e['measured_deg'])
     if old['method'] in RESOLVED and float(e['measured_deg'])>float(old['upper_deg'])+.01:contradictions.append(dict(configuration=label,case=e['case'],direction=i,old_upper=old['upper_deg'],new_feasible=e['measured_deg']))
     # Retain the furthest independently feasible endpoint. A narrower local
     # bracket must not suppress a better feasible dispatch from another start.
     if float(old['measured_deg'])>float(e['measured_deg'])+1e-8:continue
     if abs(float(old['measured_deg'])-float(e['measured_deg']))<=1e-8 and old['method'] in RESOLVED and (e['method'] not in RESOLVED or oldwidth<=newwidth):continue
    edges[key]={**e,'evidence_folder':str(p.relative_to(OUT))}
  P=K=None
  if refready:
   proof[str((folder/'reference').relative_to(OUT))]=validate(folder/'reference')
   P,K=arrays(folder/'reference',n)
   assert proof[str((folder/'reference').relative_to(OUT))]['verified_solves']==2*n-1
  for case in ('nobatt','cyclic'):
   b=next((x for x in base if x['case']==case),None)
   selected=[e for (cs,i),e in sorted(edges.items()) if cs==case]
   resolved=[e for e in selected if e['method'] in RESOLVED and math.isfinite(float(e['upper_deg']))]
   for e in selected:
    row={**e,'search_reference_predicted_deg':e['predicted_deg']}
    if refready:
     row['published_predicted_deg'],row['predicted_face']=prediction(P,K,case,np.array([float(e[f'u{i}']) for i in range(1,n)]))
    publication_edges.append(row)
   progress.append(dict(configuration=label,case=case,expected_directions=expected,written_rows=written,completed_directions=len(selected),resolved_directions=len(resolved),positive_feasible_points=sum(float(e['measured_deg'])>0 for e in selected),unresolved_directions=sum(e['method'] not in RESOLVED for e in selected),at_limit=sum(e['method']=='at_limit' for e in selected),origin_feasible=b and b['origin_feasible']=='true',reference_ready=refready))
   if strict:
    assert ready and refready,label
    assert b and b['origin_feasible']=='true',(label,case,'no common admissible origin')
    assert len(selected)==expected,(label,case,len(selected),expected)
  results[label]=dict(base=base,edges=list(edges.values()),reference_ready=refready)
 assert len(set(profile_hashes.values()))<=1,'Profiles differ within matched campaign'
 validation=dict(complete_attempt_coverage=strict,progress=progress,phase_checks=proof,profile_hashes=profile_hashes,local_bracket_disagreements=contradictions,total_verified_solves=sum(x['verified_solves'] for x in proof.values()),max_dP_kW=max([x['max_dP_kW'] for x in proof.values()] or [0]),max_dV_pu=max([x['max_dV_pu'] for x in proof.values()] or [0]),published_reference='Voltage-banded idle cost OPFs, h=0.1 degrees; original unbanded guesses retained unchanged in raw CSVs')
 (OUT/'validation_progress.json').write_text(json.dumps(validation,indent=2)+'\n',encoding='utf-8')
 (OUT/'publication_edges.json').write_text(json.dumps(publication_edges,indent=2)+'\n',encoding='utf-8')
 summary={'coverage':progress,'paired_comparisons':[],'sampling':'2 signed rays, 32 evenly spaced planar rays, 78 Fibonacci spatial rays; T=24','reference':validation['published_reference'],'local_bracket_disagreements':contradictions}
 for n in (2,3,4):
  for case in ('nobatt','cyclic'):
   aa={int(e['direction']):e for e in publication_edges if e['configuration']==f'primary{n}' and e['case']==case and e['method'] in RESOLVED and float(e['measured_deg'])>0}
   bb={int(e['direction']):e for e in publication_edges if e['configuration']==f'secondary{n}' and e['case']==case and e['method'] in RESOLVED and float(e['measured_deg'])>0}
   paired=[]
   for i in sorted(aa.keys()&bb.keys()):
    a,b=aa[i],bb[i];al,au=float(a['measured_deg']),float(a['upper_deg']);bl,bu=float(b['measured_deg']),float(b['upper_deg'])
    paired.append(dict(direction=i,primary_lower_deg=al,primary_upper_deg=au,secondary_lower_deg=bl,secondary_upper_deg=bu,relative_lower_percent=100*(bl/au-1),relative_upper_percent=100*(bu/al-1),relative_endpoint_percent=100*(bl/al-1)))
   summary['paired_comparisons'].append(dict(sources=n,case=case,paired_brackets=paired))
 (OUT/'comparison_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({k:validation[k] for k in ['total_verified_solves','max_dP_kW','max_dV_pu']},indent=2))
 for x in progress:print(x['configuration'],x['case'],str(x['resolved_directions'])+'/'+str(x['expected_directions']),'resolved; reference',x['reference_ready'])
 if strict:
  data=json.loads((ROOT/'results/multiperiod/four_viewer_data_fib78.json').read_text(encoding='utf-8'));attach_prediction_meshes(data,ROOT)
  for label,title in [('primary4','IEEE9500 primary (matched, T24)'),('secondary4','IEEE9500 + secondaries (matched, T24)'),('sf4','San Francisco P1U (T24, 78 rays)')]:
   P,K=arrays(CLUSTER/label/'reference',4);names=[f'subs{i}' for i in range(1,5)]
   ss={'id':'matched_'+label,'title':title,'sources':names,'normals':(-K.sum(axis=1)).tolist(),'offsets':P.sum(axis=1).tolist(),'linearity_max':None,'sparse':False,'study_note':'T = 24; Fibonacci-78 directions attempted for each battery case. Exact affine import-only prediction, no display cap. Predictions use voltage-banded idle references at +/-0.1 deg. Hard cyclic SOC, PV reactive dispatch, source EMFs fixed at 1.05 pu. Endpoints and local inside/outside checks replayed in OpenDSS. Local results, not global certificates. No random linearity test.','cases':{}}
   for case in ('nobatt','cyclic'):
    points=[]
    for e in publication_edges:
     if e['configuration']!=label or e['case']!=case or float(e['measured_deg'])<=0:continue
     u=np.array([float(e[f'u{i}']) for i in range(1,4)]);r=float(e['measured_deg']);bind=e['binding'];bound=-1 if bind=='voltage' else names.index(bind) if bind in names else -2
     bracketed=e['method'] in RESOLVED;up=float(e['upper_deg'])
     labeltext=('Voltage proximity' if bound==-1 else 'Minimum-import proximity: '+bind)+('; local bracket '+f"{r:.4f}-{up:.4f} deg" if bracketed else '; verified feasible lower bound, boundary unresolved')
     points.append([*(u*r).tolist(),e['predicted_face'],r,e['published_predicted_deg'],False,bound,labeltext,up if math.isfinite(up) else None,not bracketed])
    coverage=next(x for x in progress if x['configuration']==label and x['case']==case)
    ss['cases'][case]={'mesh':prediction_mesh(P,K,case),'points':points,'coverage':coverage}
   ss['study_note']+=' Verified feasible points: '+', '.join(f"{cs} {len(ss['cases'][cs]['points'])}/78 ({ss['cases'][cs]['coverage']['resolved_directions']} locally bracketed)" for cs in ('nobatt','cyclic'))+'. Crosses are feasible lower bounds with unresolved boundaries. Proximity labels are not exclusive active-set certificates.'
   data.append(ss)
  encoded=json.dumps(data,separators=(',',':'),allow_nan=False)
  (OUT/'four_viewer_data.json').write_text(encoded,encoding='utf-8')
  viewer=ROOT/'results/multiperiod/viewer';html=(viewer/'viewer_head.html').read_text(encoding='utf-8')+encoded+(viewer/'viewer_tail.html').read_text(encoding='utf-8')
  (viewer/'four_regions.html').write_text(html,encoding='utf-8');(ROOT.parents[1]/'docs/four_regions.html').write_text(html,encoding='utf-8')
if __name__=='__main__':main('--strict' in sys.argv)
