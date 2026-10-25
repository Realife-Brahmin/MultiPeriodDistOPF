"""Exact static controls for the IEEE9500 primary/service-transformer pair."""
from pathlib import Path
import json,collections,hashlib
from check_secondary_8500_inputs import definitions,fields,ROOT
A='ieee9500_1ph';B='ieee9500xf_1ph'
a=definitions(A,'BranchData.dss');b=definitions(B,'BranchData.dss')
assert a==[x for x in b if not x.lower().startswith('new line.svc_')]
for f in ('PVSystem.dss','Storage.dss','Ties.dss'):assert definitions(A,f)==definitions(B,f),f
services={}
for l in b:
 if l.lower().startswith('new line.svc_'):
  d=fields(l);services[d['bus2'].split('.')[0].lower()]=d['bus1'].split('.')[0].lower()
def load(deck):
 out=collections.defaultdict(complex)
 for l in definitions(deck,'Loads.dss'):
  d=fields(l);out[d['bus'].split('.')[0].lower()]+=complex(float(d['kw']),float(d['kvar']))
 return out
pa=load(A);raw=load(B);agg=collections.defaultdict(complex)
for bus,s in raw.items():agg[services.get(bus,bus)]+=s
err=max(abs(pa[x]-agg[x]) for x in pa.keys()|agg.keys());assert err<1e-3,err
result={'validated':True,'primary_lines':len(a),'service_transformers':len(services),'primary_delivery_proxies':len(pa),'explicit_delivery_points':len(raw),'load_aggregation_max_error_kVA':err,'identical_primary_lines_DERs_and_ties':True,'fourth_source_site':'l2804256','fourth_source_series_ohm':[.104186,.617339],'voltage_contract':'0.95-1.05 pu at all actual load delivery buses; primary attachments act as proxies when service transformers are omitted; includes direct-primary deliveries','legacy_reference':'Historical 2705 primary data additionally banded all primary transit buses; retain separately, not as the matched control','T':24,'four_source_directions':'same Fibonacci-78 for both models and San Francisco','three_source_directions':32}
out=ROOT/'problems/msopf/results/matched_secondary_oct25';out.mkdir(exist_ok=True)
(out/'input_check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
