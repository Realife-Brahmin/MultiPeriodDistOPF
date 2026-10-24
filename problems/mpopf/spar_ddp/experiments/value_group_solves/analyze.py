import csv,json,math,pathlib,re,subprocess,sys
C=pathlib.Path(sys.argv[1]);ROOT=pathlib.Path(sys.argv[2]);OUT=pathlib.Path(sys.argv[3]);OUT.mkdir(parents=True,exist_ok=True)
kv=lambda s:dict(re.findall(r'(\w+)=([^\s]+)',s))
assert (C/'campaign.finished').is_file() and (C/'smoke.complete').is_file()
assert 'Gold 6438Y+' in (C/'hardware.txt').read_text() and 'MinMemoryNode=300G' in (C/'hardware.txt').read_text()
assert not (C/'logs/slurm-31356155.err').read_text().strip()
assert all(x.endswith(': OK') for x in (C/'logs/source-check.log').read_text().splitlines())
comparison=json.loads((C/'comparison.json').read_text());assert len(comparison)==9 and all(r['validated'] and r['within_0_1_percent'] for r in comparison)
per=[];audits=[];passes=[];validated=0
for P in sorted((C/'results').iterdir()):
    if not P.is_dir():continue
    m=re.fullmatch(r'(.+)_T(\d+)_(dense|root|chunks16)',P.name);assert m
    system,T,mode=m[1],int(m[2]),m[3];assert (P/'COMPLETE').is_file()
    for script,args in [('screening_ablation/validate.py',[str(P/'solve.log'),system]),('tree_ablation/validate_opendss.py',[str(P/'opendss.log')])]:
        r=subprocess.run([sys.executable,str(ROOT/'spar_ddp/experiments'/script),*args],capture_output=True,text=True)
        assert r.returncode==0,(P,r.stdout,r.stderr);validated+=1
    text=(P/'solve.log').read_text();sol=re.search(r'solve complete: ([\d.]+) s, iterations=(\d+), status=(\d+)',text);assert sol and sol[3]=='0'
    physical=kv((P/'physical.log').read_text());tol=1e-6 if system=='ieee123C_1ph' else 1e-5
    assert math.isfinite(float(physical['objective'])) and 0<=float(physical['original_max'])<=tol
    tags=[kv(x) for x in re.findall(r'VALUE_GROUP_DIAGNOSTIC ([^\n]+)',text)];assert tags
    g=tags[0];assert all(x==g for x in tags)
    checks=[kv(x) for x in text.splitlines() if x.startswith('VALUE_REFERENCE_CHECK ')];assert checks
    for x in checks:
        for f in ('relative_Vxx','relative_Vx'):assert math.isfinite(float(x[f]))
        if mode=='dense':assert max(float(x['relative_Vxx']),float(x['relative_Vx']))<1e-7
        audits.append(dict(system=system,T=T,**x))
    ps=[kv(x) for x in text.splitlines() if x.startswith('FILTERDDP_PROCS_PASS ')];assert ps
    for p in ps:passes.append(dict(system=system,T=T,mode=mode,**p))
    per.append(dict(system=system,T=T,mode=mode,states=int(g['states']),groups=int(g['groups']),max_group=int(g['max_group']),iterations=int(sol[2]),backward_passes=len(ps),logged_audit_contaminated_s=float(sol[1]),objective=float(physical['objective']),original_primal=float(physical['original_max']),reference_checks=len(checks),audited_iterations=';'.join(str(x['iteration']) for x in checks),max_observed_relative_Vxx=max(float(x['relative_Vxx']) for x in checks),max_observed_relative_Vx=max(float(x['relative_Vx']) for x in checks),scan_wall_sum_including_reference_build_s=sum(float(p['rounds_s']) for p in ps)))
assert len(per)==12 and validated==24
for name,data in [('per_run.csv',per),('reference_checks.csv',audits),('passes.csv',passes)]:
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in data for k in r)));w.writeheader();w.writerows(data)
summary=dict(full_cases=9,smokes=3,raw_validation_checks=validated,all_pass=True,diagnostic_times_not_reportable_benchmarks=True,late_audit_coverage_incomplete_due_j0_filter=True)
(OUT/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
for r in per:
    if r['T']==24:print(json.dumps(r))
print('VALUE_GROUP_COLLECTION_VALIDATED')
