from pathlib import Path
import csv, json, statistics
from collections import defaultdict
root = Path(__file__).parent
log = (root/"rahul_merge_profile_nov05_v2/slurm-31574751.out").read_text()
header = next(s.split(",")[1:] for s in log.splitlines() if s.startswith("MERGE_OP_FIELDS,"))
rows = [dict(zip(header,s.split(",")[1:])) for s in log.splitlines() if s.startswith("MERGE_OP,")]
assert "BALANCED_COMPLETE" in log and "ERROR:" not in log
errors = [s for s in log.splitlines() if s.startswith("ERRORS ")]
assert len(errors)==12
for line in errors:
    fields=dict(x.split("=") for x in line.split()[1:])
    assert max(float(fields["grad"]),float(fields["hess"]))<1e-11
rows = [r for r in rows if int(r["rep"])>0]
assert len(rows)==5*(368+256)
assert sum(int(r["cholesky"]) for r in rows)>0
operations = header[8:]
output = {}
for mode in ("balanced","grouped"):
    subset = [r for r in rows if r["mode"]==mode]
    kernel = sum(float(r["kernel_s"]) for r in subset)
    totals = {k:sum(int(r[k]) for r in subset)/1e9 for k in operations}
    assert all(r["fallback"]=="0" for r in subset)
    print(mode, "merges",len(subset), "summed_worker_seconds", kernel)
    print(json.dumps({k:{"worker_s_per_pass":v/5,"percent_kernel":100*v/kernel} for k,v in totals.items()},indent=2))
    for two in ("true","false"):
        merges = [r for r in subset if r["two_sided"]==two]
        print("two_sided",two,"mean_merge_ms",1000*statistics.mean(float(r["kernel_s"]) for r in merges))
    rounds = defaultdict(lambda:defaultdict(list))
    for r in subset:
        rounds[(r["rep"],r["round"])][r["pid"]].append(r)
    critical = dict.fromkeys(operations,0.)
    for workers in rounds.values():
        busiest = max(workers.values(), key=lambda rs:sum(float(r["kernel_s"]) for r in rs))
        for k in operations:
            critical[k] += sum(int(r[k]) for r in busiest)/1e9/5
    print("critical_worker_op_seconds_per_pass",json.dumps(critical))
    output[mode]={"kernel_worker_seconds_per_pass":kernel/5,"operations_worker_seconds_per_pass":{k:v/5 for k,v in totals.items()},"critical_worker_seconds_per_pass":critical}
(root/"merge_operation_summary.json").write_text(json.dumps(output,indent=2))
with (root/"merge_operations.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=header);w.writeheader();w.writerows(rows)
for campaign in ("rahul_ablation_nov05","rahul_fresh_288_nov05"):
    for folder in sorted((root/campaign/"results").iterdir()):
        assert (folder/"COMPLETE").exists()
        j=json.loads((folder/"timing_summary.json").read_text())
        assert "OPENDSS_SUMMARY" in (folder/"opendss.log").read_text()
        print(folder.name,j["full_optimizer_s"],j["scan_wall_s"])
