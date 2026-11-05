"""Strict matched scoring; does not run or alter an optimizer."""
import argparse,csv,hashlib,json,math,re,statistics,subprocess,sys
from pathlib import Path
P=argparse.ArgumentParser()
P.add_argument("root",type=Path);P.add_argument("out",type=Path)
P.add_argument("--mode",choices=["t24","core128","core128t24"],required=True)
P.add_argument("--historical-central",type=Path)
a=P.parse_args();root,out=a.root,a.out;out.mkdir(parents=True,exist_ok=True)
historical=a.mode=="t24"
T=96 if a.mode=="core128" else 24
assert (root/(a.mode+".finished")).is_file()
def read(p):
 with p.open(newline="") as f:return list(csv.DictReader(f))
def save(name,data):
 assert data
 with (out/name).open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in data for k in r)));w.writeheader();w.writerows(data)
def num(x):
 x=float(x);assert math.isfinite(x);return x
def physical(folder):
 text=(folder/"physical.log").read_text()
 v={k:num(re.search(k+r"=(\S+)",text)[1]) for k in ("original_max","objective")}
 assert 0<=v["original_max"]<=1e-5
 dss=(folder/"opendss.log").read_text()
 assert "below limit: 0" in dss and "above limit: 0" in dss
 return v
jobs=list((root/"jobs").iterdir());assert len(jobs)==1
job=jobs[0];hardware=(job/"hardware.txt").read_text()
assert ("Gold 6438Y+" if historical else "6760P") in hardware
stderr=(job/"stderr.log").read_text()
# Julia emits benign precompile summaries to stderr, not numerical errors.
unclassified=[s for s in stderr.splitlines() if s.strip() and not re.fullmatch(r"\s*\d+ dependenc(?:y|ies) successfully precompiled in \d+ seconds?(?:\. \d+ already precompiled\.)?\s*",s)]
assert not unclassified,unclassified
for p in [root/"source-check-final.log",root/"harness-check-final.log",root/"inherited-source-check.log",root/"inherited-input-check.log",job/"source-check.log",job/"harness-check.log",job/"source-check-after.log",job/"harness-check-after.log"]:
 text=p.read_text();assert text.strip() and all(s.endswith(": OK") for s in text.splitlines()),p
acct=(root/"accounting.txt").read_text()
assert re.search(r"(?m)^"+re.escape(job.name)+r"\|COMPLETED\|[^|]+\|0:0\|",acct)
subprocess.run([sys.executable,str(Path(__file__).with_name("collect_ddp.py")),str(root),str(out/"ddp"),"--mode",a.mode],check=True)
validation=json.loads((out/"ddp/validation.json").read_text());assert validation["validated"] and not validation["strict_equivalence_failures"]
assert "strict_existing_match_gate = true" in (root/f"comparisons_T{T}/reference_repeats.toml").read_text()
ddp=read(out/"ddp/per_run.csv")
assert len({(r["x0_sha256"],r["ubar_sha256"]) for r in ddp})==1, "Physical start arrays differ across repeats"
central=[];initial=[]
if historical:
 assert a.historical_central, "Use existing published central reference on the usual node"
 h=[r for r in read(a.historical_central) if r["system"]=="large10kC_1ph" and r["T"]=="24" and r["solver"] in ("ma57","ma97") and r["gap_percent"]=="0.1" and r["metric"]=="original"]
 assert len(h)==4 and {r["repeat"] for r in h}=={"1","2"}
 for r in h:
  assert r["status"]=="reached" and num(r["primal"])<=1e-5 and num(r["actual_gap"])<=.001
  central.append(dict(arm=r["solver"],repeat=int(r["repeat"]),full_s=num(r["full_seconds"]),quality_s=num(r["seconds"]),reference_objective=num(r["reference_objective"]),scope="existing October16 central reference",source=r["source"]))
 (out/"historical_reference.json").write_text(json.dumps(dict(path=str(a.historical_central),sha256=hashlib.sha256(a.historical_central.read_bytes()).hexdigest()),indent=2))
 save("reused_central_rows.csv",h)
else:
 for arm in ("ma57","ma97"):
  for rep in (1,2):
   folder=root/f"results/10k_1020_T{T}_{arm}_r{arm}_{rep}"
   assert (folder/"COMPLETE").is_file()
   text=(folder/"solve.log").read_text();s=json.loads((folder/"timing_summary.json").read_text())
   assert s["central"][0]["status"]=="LOCALLY_SOLVED"
   assert s["preprocessing"][0]["matched_physical_start"]=="1"
   assert "IPOPT_WARMUP_CHECK" in text and "initial_objective_difference=0.000000000e+00" in text
   assert "MATCHED_CENTRAL_RUNTIME blas=1" in text and "omp=1" in text
   assert json.loads((folder/"logging_check.json").read_text())["logging_validated"]
   for label in ("Dual infeasibility","Constraint violation","Complementarity"):
    vals=re.findall(r"(?m)^"+label+r"\.*:\s+(\S+)\s+(\S+)",text)
    assert vals and 0<=num(vals[-1][1])<=1e-5,(folder,label,vals[-1:])
   v=physical(folder);trace=read(folder/"quality_trace.csv")
   assert trace
   for x in trace:
    for k in ("seconds","objective","original","dual"):num(x[k])
    assert num(x["original"])>=0
   assert all(num(y["seconds"])>=num(x["seconds"]) for x,y in zip(trace,trace[1:]))
   central.append(dict(arm=arm,repeat=rep,full_s=num(s["central"][0]["solve_time_s"]),objective=v["objective"],independent_original_primal=v["original_max"],trace=trace,folder=folder,scope="fresh same-node reference"))
   initial.append(num(trace[0]["objective"]))
reference=statistics.mean(r["reference_objective"] if historical else r["objective"] for r in central)
if not historical:assert (max(r["objective"] for r in central)-min(r["objective"] for r in central))/max(1,abs(reference))<=1e-5
per=[]
for row in ddp:
 folder=root/f"results/10k_1020_T{T}_ddp_r{row['arm']}_{row['repeat']}"
 trace=read(folder/"quality_trace.csv");initial.append(num(trace[0]["objective"]))
 assert abs(num(row["objective"])-reference)/max(1,abs(reference))<=.001
 for x in trace:
  for k in ("seconds","objective","original","dual"):num(x[k])
  assert num(x["original"])>=0
 assert all(num(y["seconds"])>=num(x["seconds"]) for x,y in zip(trace,trace[1:]))
 hits=[x for x in trace if num(x["original"])<=1e-5 and abs(num(x["objective"])-reference)/max(1,abs(reference))<=.001]
 assert hits and num(trace[-1]["seconds"])<=num(row["full_s"])+.01
 per.append(dict(**row,quality_s=num(hits[0]["seconds"]),quality_iteration=hits[0]["iteration"],quality_original=num(hits[0]["original"]),reference_objective=reference,reference_scope="historical central" if historical else "fresh same node"))
for row in central:
 if not historical:
  trace=row.pop("trace");folder=row.pop("folder")
  hits=[x for x in trace if num(x["original"])<=1e-5 and abs(num(x["objective"])-reference)/max(1,abs(reference))<=.001]
  assert hits and num(trace[-1]["seconds"])<=row["full_s"]+.01
  row["quality_s"]=num(hits[0]["seconds"])
 per.append(row)
assert (max(initial)-min(initial))/max(1,abs(initial[0]))<1e-10,"Physical starting objectives differ"
med=[]
for arm in dict.fromkeys(r["arm"] for r in per):
 selected=[r for r in per if r["arm"]==arm];assert len(selected)==2
 med.append(dict(arm=arm,repeats=2,full_s=statistics.median(num(r["full_s"]) for r in selected),quality_s=statistics.median(num(r["quality_s"]) for r in selected),full_min=min(num(r["full_s"]) for r in selected),full_max=max(num(r["full_s"]) for r in selected)))
centr=[r for r in med if r["arm"] in ("ma57","ma97")]
full=min(centr,key=lambda r:r["full_s"]);quality=min(centr,key=lambda r:r["quality_s"])
cmp=[dict(arm=r["arm"],full_central_backend=full["arm"],full_speedup=full["full_s"]/r["full_s"],quality_central_backend=quality["arm"],quality_speedup=quality["quality_s"]/r["quality_s"]) for r in med if r not in centr]
save("per_run.csv",per);save("medians.csv",med);save("comparisons.csv",cmp)
(out/"validation.json").write_text(json.dumps(dict(validated=True,mode=a.mode,objective_gate=.001,physical_gate=1e-5,strict_dispatch_gate=1e-10,stderr_benign_precompile_only=bool(stderr.strip()),central_scope="existing October16" if historical else "fresh same node"),indent=2))
print(json.dumps(cmp,indent=2))

