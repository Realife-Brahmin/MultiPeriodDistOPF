# Required scientific gate between stages of the immutable launcher.
import json,math,re,sys
from pathlib import Path
p=Path(sys.argv[1]);debug=len(sys.argv)>2 and sys.argv[2]=='debug'
assert (p/'COMPLETE').is_file()
s=(p/'solve.log').read_text();summary=json.loads((p/'timing_summary.json').read_text())
assert 'WARMUP full' in s and 'BR_DEBUG_FAIL' not in s and 'FILTERDDP_PROCS_FAIL' not in s
assert not re.search(r'ERROR:|AssertionError|LoadError',s)
m=re.search(r'final residuals: primal=(\S+) dual=(\S+) complementarity=(\S+)',s)
assert m and all(math.isfinite(float(x)) and 0<=float(x)<=1e-5 for x in m.groups())
physical=(p/'physical.log').read_text();x=float(re.search(r'original_max=(\S+)',physical)[1]);assert 0<=x<=1e-5
s2=(p/'opendss.log').read_text();assert 'below limit: 0' in s2 and 'above limit: 0' in s2
passes=re.findall(r'^BR_PASS solve=1 .*$',s,re.M);assert passes
candidate='requested=true' in '\n'.join(passes)
if 'cache' in p.name or debug:
 assert candidate and any('hit=true' in x for x in passes)
 assert all('debug='+str(debug).lower() in x for x in passes)
 assert 'BR_RETIRE bytes=' in s and 'collected=true' in s
else:
 assert all('hit=false' in x and 'reused=0' in x for x in passes)
if debug:
 assert 'BR_DEBUG_CHECK label=actual_forward_' in s
 assert 'BR_DEBUG_CHECK label=actual_' in s
 assert 'BR_NEGATIVE stage=24' in s
 assert re.search(r'BR_STALE_VECTOR .*rejected=true',s)
print('BR_SCIENTIFIC_GATE_PASS',p.name,'debug='+str(debug))
