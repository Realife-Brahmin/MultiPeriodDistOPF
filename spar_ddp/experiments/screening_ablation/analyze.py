"""Read-only revalidation; write new analysis directory outside frozen source."""
import csv,itertools,json,math,pathlib,re,statistics,subprocess,sys
C=pathlib.Path(sys.argv[1]); O=pathlib.Path(sys.argv[2]); O.mkdir(exist_ok=False)
H=C/'source/spar_ddp/experiments'; rows=[]
def match(p,t):
    m=re.search(p,t); assert m,p; return m
for p in sorted(list((C/'smoke').glob('m*'))+list((C/'results').glob('*_T96_r*_m*'))):
    assert (p/'arm.complete').is_file(), str(p)
    text=(p/'solve.log').read_text()
    h=match(r'SCREEN_ABLATION_RUN system=(\S+) T=(\d+) repeat=(\d+) mask=(\d+) screen=(\S+)',text)
    subprocess.run([sys.executable,str(H/'screening_ablation/validate.py'),str(p/'solve.log'),h[1]],check=True,capture_output=True)
    subprocess.run([sys.executable,str(H/'tree_ablation/validate_opendss.py'),str(p/'opendss.log')],check=True,capture_output=True)
    orig=re.findall(r'^FILTERDDP_ORIGINAL_PRIMAL iteration=(\d+) barrier_iteration=(\d+) original_max=(\S+) primal_inf=(\S+)',text,re.M)
    passes=re.findall(r'^FILTERDDP_PROCS_PASS ',text,re.M)
    assert len(orig)==len(passes)>0 and all(math.isfinite(float(v)) for r in orig for v in r[2:]),str(p)
    assert float(orig[-1][2]) <= (1e-6 if h[1]=='ieee123C_1ph' else 1e-5)
    if int(h[2])!=96: continue
    s=match(r'solve complete: ([\d.]+) s, iterations=(\d+), status=0\b',text)
    r=match(r'final residuals: primal=(\S+) dual=(\S+) complementarity=(\S+)',text)
    b=match(r'BOUND_SCREEN .*voltage limits dropped (\d+) lower, (\d+) upper of (\d+) buses; ell >= 0 dropped on (\d+) lines; P_Subs >= 0 dropped: (\w+)',text)
    ph=[{k:float(v) for k,v in re.findall(r'\b(step1_s|rounds_s|step3_s)=([\deE.+-]+)',line)} for line in text.splitlines() if line.startswith('FILTERDDP_PROCS_PASS ')]
    rows.append(dict(system=h[1],T=int(h[2]),repeat=int(h[3]),mask=int(h[4]),screen=h[5],optimizer_s=float(s[1]),iterations=int(s[2]),objective=float(match(r'FilterDDP objective=(\S+)',text)[1]),primal=float(r[1]),dual=float(r[2]),complementarity=float(r[3]),original_max=float(orig[-1][2]),lower_dropped=int(b[1]),upper_dropped=int(b[2]),ell_dropped=int(b[4]),import_dropped=b[5],phase_A_s=sum(x['step1_s'] for x in ph),phase_B_s=sum(x['rounds_s'] for x in ph),phase_C_s=sum(x['step3_s'] for x in ph)))
systems=sorted({r['system'] for r in rows})
assert len(rows)==48 and len(systems)==3
assert {(r['system'],r['repeat'],r['mask']) for r in rows}==set(itertools.product(systems,(1,2),range(8)))
def csvout(name,data):
    with (O/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
csvout('per_run.csv',rows)
med=[]; margins=[]; summary={}
for sysname in systems:
    subset=[r for r in rows if r['system']==sysname]; objs=[r['objective'] for r in subset]
    gap=(max(objs)-min(objs))/max(1,abs(min(objs)),abs(max(objs)))
    summary[sysname]={'objective_relative_range':gap,'objective_range':[min(objs),max(objs)]}
    times={m:statistics.median(r['optimizer_s'] for r in subset if r['mask']==m) for m in range(8)}
    for m in range(8):
        a=[r for r in subset if r['mask']==m]
        med.append(dict(system=sysname,mask=m,screen=a[0]['screen'],optimizer_s=times[m],iterations=statistics.median(r['iterations'] for r in a),repeat_min_s=min(r['optimizer_s'] for r in a),repeat_max_s=max(r['optimizer_s'] for r in a),baseline_over_mask=times[0]/times[m]))
    for bit,label in ((1,'vupper'),(2,'ell'),(4,'psubs')):
        for m in range(8):
            if m&bit: continue
            margins.append(dict(system=sysname,rule=label,other_mask=m,off_s=times[m],on_s=times[m|bit],speedup=times[m]/times[m|bit],saved_percent=100*(1-times[m|bit]/times[m])))
csvout('medians.csv',med);csvout('conditional_effects.csv',margins)
(O/'validation.json').write_text(json.dumps({'validated_main_arms':48,'validated_smoke_arms':8,'per_system':summary},indent=2)+'\n')
print(json.dumps(summary,indent=2))
for r in med: print(r)
