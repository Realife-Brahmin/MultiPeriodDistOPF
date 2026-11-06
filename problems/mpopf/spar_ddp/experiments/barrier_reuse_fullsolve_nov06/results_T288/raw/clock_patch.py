"""Add passive elapsed timestamps to a new private benchmark snapshot."""
from pathlib import Path
import difflib,hashlib,json,sys
root=Path(sys.argv[1]);source=root/'variants/10k_1020/source';NL=chr(10)
paths=list(source.glob('spar_ddp/*/src/solve.jl'));assert len(paths)==1
module='DDP4OPF' if 'DDP4OPF' in str(paths[0]) else 'SPARDDP'
changed=[]
def rewrite(path,old,new):
 assert old!=new
 relative=path.relative_to(source).as_posix()
 saved=root/'clock-originals'/relative;saved.parent.mkdir(parents=True,exist_ok=True)
 assert not saved.exists();saved.write_bytes(path.read_bytes())
 path.write_text(new,encoding='utf-8')
 changed.append(dict(file=relative,before_sha256=hashlib.sha256(old.encode()).hexdigest(),after_sha256=hashlib.sha256(new.encode()).hexdigest(),diff=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=relative,tofile=relative))))
p=paths[0];old=p.read_text(encoding='utf-8');assert 'STUDY_OPTIMIZER_START_NS' not in old
header=NL.join(('const STUDY_OPTIMIZER_START_NS=Ref{UInt64}(0)','study_clock_start!()=(STUDY_OPTIMIZER_START_NS[]=time_ns();time())',''))
needle='        backward_s = (time_ns() - backward_start_ns) / 1e9'+NL
assert old.count(needle)==1
marker=NL.join(('        if timing_diagnostic && STUDY_OPTIMIZER_START_NS[] != 0','            @printf("FILTERDDP_ELAPSED_CLOCK iteration=%d barrier_iteration=%d elapsed_s=%.9f'+chr(92)+'n",','                    data.k, data.j, (time_ns()-STUDY_OPTIMIZER_START_NS[])/1e9)','        end',''))
rewrite(p,old,header+old.replace(needle,needle+marker))
p=source/'spar_ddp/examples/power_system/ieee123c_filterddp.jl';old=p.read_text(encoding='utf-8');new=old
for variable,indent in (('tw','    '),('t1','')):
 needle=f'{indent}{variable} = time()'+NL;assert new.count(needle)==1
 new=new.replace(needle,f'{indent}{variable} = {module}.study_clock_start!()'+NL)
rewrite(p,old,new)
(root/'clock-source-changes.json').write_text(json.dumps(dict(scope='passive elapsed clock only; identical baseline/candidate source',changed=changed),indent=2),encoding='utf-8')
(root/'clock-source.patch').write_text(''.join(x['diff'] for x in changed),encoding='utf-8')
old_entries={}
for line in (root/'source-before-clock.sha256').read_text().splitlines():
 digest,name=line.split('  ',1);old_entries[name]=digest
allowed={x['file'] for x in changed};different=[]
for name,digest in old_entries.items():
 current=hashlib.sha256((source/name).read_bytes()).hexdigest()
 if current!=digest:different.append(name)
assert set(different)==allowed,(different,allowed)
print(json.dumps(dict(exactly_two_timing_only_source_sites=True,changed_files=sorted(different),all_other_input_model_kernel_files_identical=True)))
