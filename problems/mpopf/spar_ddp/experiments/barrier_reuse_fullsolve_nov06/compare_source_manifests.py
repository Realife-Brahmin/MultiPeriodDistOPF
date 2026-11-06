from pathlib import Path
import sys,json
def manifest(path):
 result={}
 for line in path.read_text().splitlines():
  digest,name=line.split('  ',1);result[name]=digest
 return result
old=manifest(Path(sys.argv[1]));new=manifest(Path(sys.argv[2]));T=int(sys.argv[3])
assert set(old)<=set(new)
assert all(old[k]==new[k] for k in old), 'Numerical source/input differs from validated v4'
extra=set(new)-set(old)
allowed={f'spar_ddp/results/network_filterddp/network_data_large10kC_1ph_T{T}_periodic.jls'}
assert extra<=allowed,extra
print(json.dumps(dict(exact_common_source_and_inputs=True,added_inputs=sorted(extra),T=T)))
