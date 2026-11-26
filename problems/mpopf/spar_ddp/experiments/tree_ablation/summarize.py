import csv,pathlib,re,statistics,sys
root=pathlib.Path(sys.argv[1]);rows=[]
for folder in sorted(root.glob('*_T96_r*')):
    if not (folder/'pair.complete').exists(): raise SystemExit(f'Incomplete pair: {folder}')
    for arm in ('tree','umfpack'):
        text=(folder/f'{arm}.log').read_text()
        head=re.search(r'TREE_ABLATION_RUN system=(\S+) T=(\d+) repeat=(\d+) backend=(\S+)',text)
        solve=re.search(r'solve complete: ([\d.]+) s, iterations=(\d+), status=0\b',text)
        if head is None or solve is None: raise SystemExit(f'Missing validated solve: {folder}/{arm}')
        # Phase durations are elapsed coordinator times, not summed worker times.
        phases=[]
        for line in text.splitlines():
            if line.startswith('FILTERDDP_PROCS_PASS '):
                phases.append({k:float(v) for k,v in re.findall(r'\b(step1_s|rounds_s|step3_s)=([\d.eE+\-]+)',line)})
        if not phases: raise SystemExit('Missing process phase timing')
        rows.append(dict(system=head[1],T=int(head[2]),repeat=int(head[3]),backend=arm,
            full_optimizer_s=float(solve[1]),iterations=int(solve[2]),
            phase_A_s=sum(p['step1_s'] for p in phases),phase_B_s=sum(p['rounds_s'] for p in phases),
            phase_C_s=sum(p['step3_s'] for p in phases)))
if len(rows)!=18: raise SystemExit(f'Expected 18 arms (3 systems x 3 repeats x 2), found {len(rows)}')
writer=csv.writer(sys.stdout);writer.writerow(['system','T','metric','tree_median','umfpack_median','umfpack_over_tree'])
for system in sorted({r['system'] for r in rows}):
    for metric in ('full_optimizer_s','phase_A_s','phase_B_s','phase_C_s','iterations'):
        med={a:statistics.median(r[metric] for r in rows if r['system']==system and r['backend']==a) for a in ('tree','umfpack')}
        writer.writerow([system,96,metric,med['tree'],med['umfpack'],med['umfpack']/med['tree'] if med['tree'] else ''])
with (root/'per_run.csv').open('x',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
