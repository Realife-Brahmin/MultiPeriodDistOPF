import csv,json,re,statistics,sys,hashlib
from pathlib import Path

root=Path(sys.argv[1])
out=Path(sys.argv[2])
out.mkdir(parents=True,exist_ok=True)
assert (root/"study.finished").exists()
for line in (root/"source.sha256").read_text().splitlines():
    expected,name=line.split(maxsplit=1)
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected,name
if (root/"scratch-after-construction.log").exists():
    scratch=(root/"scratch-after-construction.log").read_text()
    assert "SCRATCH_COEFFICIENT_BYTES 24977760" in scratch
    assert "SCRATCH_OBJECT_BYTES_AFTER_CONSTRUCTION 24977976" in scratch

assert (root/"unit.log").read_text().strip().endswith("RIGHT_KERNEL_SYNTHETIC_PASS")
assert "LATE_COMPLETE accurate=true" in (root/"late-1.log").read_text()
for name in ("source-check.log","capture-check.log","late-runtime-check.log"):
    text=(root/name).read_text()
    assert "FAILED" not in text and ": OK" in text,name

def kv(line):
    return dict(re.findall(r"(\w+)\s*=\s*([^\s]+)",line))
rows=[]
rounds=[]
diagnostics=[]
errors=[]
for config in ("kernel63","threads31","threads15","combination"):
    folder=root/config
    if not folder.exists():continue
    text=(folder/"replay.log").read_text()
    assert "BALANCED_COMPLETE" in text and "ERROR:" not in text, config
    topology={}
    for line in (folder/"cpu_topology.csv").read_text().splitlines():
        if not line.startswith("#"):
            cpu,core,socket=map(int,line.split(","))
            topology[cpu]=(socket,core)
    widths={}
    for maps in folder.glob("cpu_groups_*.csv"):
        b=int(maps.stem.split("_")[-1])
        groups=[list(map(int,s.split(","))) for s in maps.read_text().splitlines()]
        cpus=[cpu for group in groups for cpu in group]
        assert len(set(cpus))==len(cpus)<=64 and len(groups[0])==1
        assert len({topology[c] for c in cpus})==len(cpus)
        assert all(len({topology[c][0] for c in group})==1 for group in groups)

        assert all(len(g)==b for g in groups[1:])
        widths[b]=len(groups)-1
    modepasses={}
    last_diag=None
    for line in text.splitlines():
        if line.startswith("CONFIG "):
            d=kv(line)
            mode=d["mode"]
            b=int(mode.split("_")[1]) if "_" in mode else (next(iter(widths)))
            assert int(d["blas_threads"].rstrip(",)"))==b,line
        elif line.startswith("ROUND_DETAIL "):
            d=kv(line);last_diag=d
            if d["rep"]=="-1":
                assert float(d["compile_s"])==0
                assert int(d["allocation_bytes"])>=0
                diagnostics.append(dict(config=config,**d))
        elif line.startswith("ROUND "):
            d=kv(line)
            tasks=list(map(int,d["tasks_per_worker"].split(",")))
            b=int(d["mode"].split("_")[1]) if "_" in d["mode"] else next(iter(widths))
            assert len(tasks)==widths[b] and sum(tasks)==int(d["merges"])
            rounds.append(dict(config=config,**d,max_tasks=max(tasks),min_tasks=min(tasks),
                               busiest_s=last_diag["busiest_s"],kernel_sum_s=last_diag["kernel_sum_s"],
                               completion_spread_s=last_diag["completion_spread_s"]))
        elif line.startswith("ERRORS "):
            d=kv(line)
            assert max(float(d["grad"]),float(d["hess"]))<1e-11
            errors.append(dict(config=config,**d))
        elif line.startswith("BALANCED_PASS "):
            d=kv(line);rep=int(d["rep"]);mode=d["scheduling"]
            assert float(d["max_relative_error"])<1e-11
            if rep>0:
                modepasses.setdefault(mode,[]).append((rep,float(d["wall_s"])))
    for mode,values in modepasses.items():
        assert sorted(rep for rep,_ in values)==[1,2,3,4,5]
        b=int(mode.split("_")[1]) if "_" in mode else next(iter(widths))
        selected=[d for d in diagnostics if d["config"]==config and d["mode"]==mode]
        assert len(selected)==8
        rows.append(dict(config=config,kernel=mode.split("_")[0],workers=widths[b],blas_threads=b,
                         median_s=statistics.median(t for _,t in values),min_s=min(t for _,t in values),max_s=max(t for _,t in values),
                         diagnostic_alloc_bytes=sum(int(d["allocation_bytes"]) for d in selected),
                         diagnostic_gc_s=sum(float(d["gc_s"]) for d in selected),
                         mean_merge_ms=1000*statistics.mean(sum(float(d["kernel_sum_s"]) for d in rounds if d["config"]==config and d["mode"]==mode and int(d["rep"])==rep)/256 for rep,_ in values),
                         passes=json.dumps(values)))
    assert len([e for e in errors if e["config"]==config])==8*len(modepasses)

for b in (1,2,4):
    path=root/f"late-{b}.log"
    if not path.exists():continue
    text=path.read_text()
    assert "LATE_COMPLETE accurate=true" in text
    for line in text.splitlines():
        if line.startswith("LATE_DIRECTION ") or line.startswith("LATE_BATTERY "):
            d=kv(line);assert float(d["error"])<1e-7 and float(d["relative"])<1e-5
for name,data in (("summary",rows),("rounds",rounds),("diagnostics",diagnostics),("errors",errors)):
    with (out/f"{name}.csv").open("w",newline="") as f:
        fields=list(dict.fromkeys(k for r in data for k in r))
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(data)
(out/"validation.json").write_text(json.dumps(dict(validated=True,configurations=len(rows),all_firstpass_gate=1e-11,late_scaled_gate=1e-7,late_relative_gate=1e-5,zero_compile_in_allocation_diagnostics=True),indent=2))
print(json.dumps(rows,indent=2))
