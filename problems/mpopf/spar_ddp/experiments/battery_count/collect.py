"""Read-only scoring of an extracted immutable battery-count evidence archive."""
import csv, json, math, pathlib, re, statistics, sys

C=pathlib.Path(sys.argv[1]); OUT=pathlib.Path(sys.argv[2]); OUT.mkdir(parents=True,exist_ok=True)
def rows(p):
    with p.open(newline='') as f: return list(csv.DictReader(f))
def save(name,data):
    assert data
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in data for k in r)));w.writeheader();w.writerows(data)
def finite(x):
    v=float(x);assert math.isfinite(v);return v
def fields(line):return dict(re.findall(r'(\w+)=([^\s]+)',line))
submitted=rows(C/'submitted.csv'); variants={(r['variant'],int(r['T'])):r for r in rows(C/'variants.csv')}
assert (C/'smoke.complete').is_file() and (C/'construction.complete').is_file()
assert len(rows(C/'battery_selection.csv'))==1020
assert sum(r['in_X']=='true' for r in rows(C/'battery_selection.csv'))==250
assert sum(r['in_Y']=='true' for r in rows(C/'battery_selection.csv'))==100
jobids={(r['variant'],int(r['T'])):r['job'] for r in submitted if r['variant']!='prepare'}
assert len(jobids)==8
for v in ('10k_X','10k_Y'):
    s=json.loads((C/'jobs'/submitted[0]['job']/f'{v}-smoke-validation.json').read_text())
    assert s['matched_quality_gate']
    for m in ('ma97','ma57','ddp'):assert (C/'results'/f'{v}_T3_{m}_r0'/'COMPLETE').is_file()

per=[]; mem=[]; validation=[]; medians=[]; comparisons=[]; traces={}
for (v,T),jid in jobids.items():
    J=C/'jobs'/jid; assert (J/'COMPLETE').is_file()
    group=json.loads((J/'group-validation.json').read_text());assert group['matched_quality_gate']
    hardware=(J/'hardware.txt').read_text();assert 'Gold 6438Y+' in hardware and 'NumCPUs=64' in hardware
    assert 'MinMemoryNode=300G' in hardware
    for name in ('source-check.txt','input-check.txt'):
        text=(J/name).read_text();assert text.strip() and all(x.endswith(': OK') for x in text.splitlines()),name
    assert not (J/'stderr.log').read_text().strip()
    cm=dict(re.findall(r'(memory\.\w+)\s+(\d+)',(J/'node_memory.txt').read_text()))
    peak=int(cm['memory.max_usage_in_bytes']);limit=int(cm['memory.limit_in_bytes']);assert peak<limit==300*1024**3
    mem.append(dict(variant=v,T=T,job=jid,job_cgroup_peak_bytes=peak,job_cgroup_peak_GiB=peak/1024**3,allocation_GiB=300,scope='entire six-arm allocation; not backend-specific'))
    central=[finite(val) for k,val in group['objectives'].items() if k.startswith('ma')]
    ref=statistics.mean(central);assert (max(central)-min(central))/max(abs(ref),1)<=1e-5
    initial=[]
    for method in ('ddp','ma97','ma57'):
        for repeat in (1,2):
            P=C/'results'/f'{v}_T{T}_{method}_r{repeat}';assert (P/'COMPLETE').is_file()
            text=(P/'solve.log').read_text();summary=json.loads((P/'timing_summary.json').read_text())
            phys=fields((P/'physical.log').read_text());assert finite(phys['original_max'])<=1e-5
            assert abs(finite(phys['objective'])-ref)/max(abs(ref),1)<=.001
            op=(P/'opendss.log').read_text();assert 'OPENDSS_SUMMARY' in op
            assert not re.search(r'UNHANDLED TASK ERROR|LoadError|OutOfMemoryError|Segmentation fault',text)
            trace=rows(P/'quality_trace.csv'); assert trace
            times=[finite(r['seconds']) for r in trace];assert all(b>=a for a,b in zip(times,times[1:]))
            for r in trace:
                for k in ('objective','original','dual'):finite(r[k])
                assert finite(r['original'])>=0
            initial.append(finite(trace[0]['objective']))
            if method=='ddp':
                match=re.search(r'solve complete: ([\d.]+) s, iterations=(\d+), status=(\d+)',text);assert match and match[3]=='0'
                full=finite(summary['full_optimizer_s']);iterations=int(match[2]);passes=rows(P/'passes.csv');clock=rows(P/'iteration_timing.csv')
                assert len(passes)==len(clock)==len(trace)
                assert abs(full-sum(finite(r['total_s']) for r in clock)-finite(summary['timer_remainder_s']))<1e-7
                assert finite(summary['timer_remainder_s'])>=-.005
                for a,b,t in zip(passes,clock,trace):
                    assert (a['iteration'],a['barrier_iteration'])==(b['iteration'],b['barrier_iteration'])==(t['iteration'],t['barrier_iteration'])
                scan=[finite(a['rounds_s']) for a in passes]
                extra=dict(backward_passes=len(passes),phase_A_s=finite(summary['phase_A_s']),scan_wall_s=sum(scan),phase_C_s=finite(summary['phase_C_s']),scan_per_pass_mean_s=statistics.mean(scan),scan_per_pass_median_s=statistics.median(scan),scan_per_pass_min_s=min(scan),scan_per_pass_max_s=max(scan),scan_fraction_of_full=sum(scan)/full,coordinator_maxrss_MiB=summary['memory'][0]['coordinator_maxrss_mib'],worker_maxrss_MiB_max=summary['memory'][0]['worker_maxrss_mib_max'],worker_maxrss_MiB_sum=summary['memory'][0]['worker_maxrss_mib_sum'])
            else:
                centralrow=summary['central'][0];assert centralrow['status']=='LOCALLY_SOLVED'
                full=finite(centralrow['solve_time_s']);iterations=int(centralrow['iterations'])
                assert 'IPOPT_WARMUP_CHECK' in text and summary['preprocessing'][0]['matched_physical_start']=='1'
                internal=(P/'ipopt_internal.log').read_text();assert 'ComputeSearchDirection' in internal
                extra=dict(build_s=finite(centralrow['build_s']),runner_wall_s=finite(centralrow['wall_s']),process_maxrss_MiB=summary['memory'][0]['maxrss_mib'])
            assert times[-1]<=full+.01
            hits=[r for r in trace if abs(finite(r['objective'])-ref)/max(abs(ref),1)<=.001 and finite(r['original'])<=1e-5]
            assert hits,(v,T,method,repeat)
            hit=hits[0];quality=finite(hit['seconds'])
            tracekey=(v,T,method);trajectory=[(r['iteration'],r.get('barrier_iteration'),finite(r['objective']),finite(r['original'])) for r in trace]
            if tracekey in traces:
                a=traces[tracekey];assert len(a)==len(trajectory)
                for x,y in zip(a,trajectory):
                    assert x[:2]==y[:2] and abs(x[2]-y[2])<=1e-9*max(abs(x[2]),1) and abs(x[3]-y[3])<=1e-9*max(abs(x[3]),1)
            traces[tracekey]=trajectory
            totalwall=finite((P/'wall.txt').read_text().strip().split('=')[1])
            per.append(dict(variant=v,T=T,batteries=int(variants[v,T]['batteries']),method=method,repeat=repeat,job=jid,reference_objective=ref,final_objective=finite(phys['objective']),final_original_primal=finite(phys['original_max']),iterations=iterations,full_optimizer_s=full,quality_s=quality,quality_iteration=hit['iteration'],quality_barrier_iteration=hit.get('barrier_iteration',''),quality_original_primal=finite(hit['original']),quality_relative_objective_gap=abs(finite(hit['objective'])-ref)/max(abs(ref),1),total_process_wall_including_validation_s=totalwall,**extra))
    assert (max(initial)-min(initial))/max(abs(initial[0]),1)<1e-10
    validation.append(dict(variant=v,T=T,solves=6,reference=ref,initial_objective_relative_range=(max(initial)-min(initial))/max(abs(initial[0]),1),repeat_trajectories_agree=True,physical_and_opendss_checks=True,source_input_checks=True))
    for method in ('ddp','ma97','ma57'):
        selected=[r for r in per if r['variant']==v and r['T']==T and r['method']==method];assert len(selected)==2
        r=dict(variant=v,T=T,method=method)
        for k in ('full_optimizer_s','quality_s','iterations','total_process_wall_including_validation_s'):
            r[k]=statistics.median(x[k] for x in selected)
        r['full_repeat_spread_fraction']=abs(selected[0]['full_optimizer_s']-selected[1]['full_optimizer_s'])/r['full_optimizer_s']
        if method=='ddp':
            for k in ('backward_passes','phase_A_s','scan_wall_s','phase_C_s','scan_fraction_of_full','scan_per_pass_mean_s','scan_per_pass_median_s','scan_per_pass_min_s','scan_per_pass_max_s'):r[k]=statistics.median(x[k] for x in selected)
        medians.append(r)
    own=[r for r in medians if r['variant']==v and r['T']==T];ddp=next(r for r in own if r['method']=='ddp');centr=[r for r in own if r['method']!='ddp']
    f=min(centr,key=lambda r:r['full_optimizer_s']);q=min(centr,key=lambda r:r['quality_s'])
    comparisons.append(dict(variant=v,T=T,batteries=int(variants[v,T]['batteries']),ddp_full_s=ddp['full_optimizer_s'],central_full_backend=f['method'],central_full_s=f['full_optimizer_s'],full_speedup=f['full_optimizer_s']/ddp['full_optimizer_s'],ddp_quality_s=ddp['quality_s'],central_quality_backend=q['method'],central_quality_s=q['quality_s'],quality_speedup=q['quality_s']/ddp['quality_s'],shared_buffer_bytes=int(variants[v,T]['shared_bytes_all_time_steps'])))
assert len(per)==48 and len(medians)==24 and len(comparisons)==8
save('per_run.csv',per);save('medians.csv',medians);save('comparison.csv',comparisons);save('job_memory.csv',mem)
(OUT/'validation.json').write_text(json.dumps(dict(measured_solves=48,smoke_solves=6,groups=validation,all_pass=True,quality_objective_relative_tolerance=.001,original_primal_tolerance=1e-5,quality_times_observed_not_interpolated=True,ddp_clock='reconstructed from timed iteration rows and backward-pass endpoint',central_clock='IPOPT_ITER callback, original primal from printed internal table'),indent=2)+'\n')
for r in comparisons:print(json.dumps(r))
print('COLLECTION_VALIDATED solves=48 smoke=6 groups=8')
