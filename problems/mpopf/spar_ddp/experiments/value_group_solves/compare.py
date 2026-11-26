import json,pathlib,re,sys,math
C=pathlib.Path(sys.argv[1]);systems=('ieee123C_1ph','ieee2522C_1ph','large10kC_1ph') if sys.argv[2]=='all' else (sys.argv[2],);T=int(sys.argv[3]);out=[]
for system in systems:
    base=C/'results'/f'{system}_T{T}_dense';assert (base/'COMPLETE').is_file()
    ref=float(re.search(r'objective=(\S+)',(base/'physical.log').read_text())[1])
    for mode in ('dense','root','chunks16'):
        p=C/'results'/f'{system}_T{T}_{mode}';row=dict(system=system,T=T,mode=mode,validated=(p/'COMPLETE').is_file(),reference_objective=ref)
        if row['validated']:
            objective=float(re.search(r'objective=(\S+)',(p/'physical.log').read_text())[1]);gap=abs(objective-ref)/max(1,abs(ref));assert math.isfinite(gap)
            row.update(objective=objective,relative_gap=gap,within_0_1_percent=gap<=.001)
            if T==3:assert gap<=.001,'smoke objective gate'
        out.append(row)
print(json.dumps(out,indent=2))
