import pathlib,re,sys,json,math
C=pathlib.Path(sys.argv[1]);V=sys.argv[2];T=int(sys.argv[3]);repeats=(0,) if T==3 else (1,2)
vals={}
for r in repeats:
 for m in ('ddp','ma57','ma97'):
  p=C/'results'/f'{V}_T{T}_{m}_r{r}';assert (p/'COMPLETE').is_file(),p
  text=(p/'physical.log').read_text();hit=re.search(r'original_max=(\S+) objective=(\S+)',text);assert hit
  assert math.isfinite(float(hit[2])) and float(hit[1])<=1e-5
  vals[f'{m}_r{r}']=float(hit[2])
central=[v for k,v in vals.items() if k.startswith('ma')]
spread=(max(central)-min(central))/max(abs(max(central)),abs(min(central)),1)
assert spread<=1e-5,('central disagreement',vals)
assert (max(vals.values())-min(vals.values()))/max(abs(max(vals.values())),abs(min(vals.values())),1)<=.001,('DDP/central disagreement',vals)
print(json.dumps({'variant':V,'T':T,'objectives':vals,'central_relative_range':spread,'matched_quality_gate':True},indent=2))
