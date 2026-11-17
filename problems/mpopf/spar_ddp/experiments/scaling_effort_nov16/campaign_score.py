"""Time to near-optimality for the two-repeat runs (kamiak/campaign_job.sh), with
the October 16 campaign's rules (near_opt_completed.py there):

- reference objective: the final objective of the first MA97 repeat of the same
  system and horizon (this campaign's for systems it ran Ipopt on; otherwise the
  October campaign's, results/proper_logging_campaign_oct16/near_opt_completed.csv,
  whose Ipopt rows are reused unchanged for those systems);
- DDP: iterate k is known at the end of pass k's backward sweep, at
  setup + sum(total_s of earlier passes) + backward_s of pass k, setup = reported
  solve time - sum(all total_s); residuals from FILTERDDP_ORIGINAL_PRIMAL
  (original model) and the iteration's primal_inf (internal);
- Ipopt: IPOPT_ITER lines (elapsed_s), original residual from Ipopt's own table;
- a target is reached at the first iterate with relative gap <= gap AND residual
  <= threshold (1e-6 ieee123, 1e-5 otherwise);
- two-repeat medians, and the faster of MA57/MA97 per case.

    python campaign_score.py <job result dir> [...]   (dirs holding summary.txt and logs)
Writes campaign_rows.csv and prints the median table at the 0.1% / original rule.
"""
import csv, gzip, re, statistics, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OCT = HERE.parents[1] / "results" / "proper_logging_campaign_oct16" / "near_opt_completed.csv"
kv = lambda s: dict(re.findall(r"(\w+)=([^\s|]+)", s))


def read(p):
    p = Path(p)
    if not p.exists() and Path(str(p) + ".gz").exists():
        p = Path(str(p) + ".gz")
    return gzip.open(p, "rt", errors="ignore").read() if p.suffix == ".gz" else p.read_text(errors="ignore")


def ddp_records(text):
    if "WARMUP" in text:                                  # the timed solve only
        text = text[text.index("\nWARMUP"):]
    lines = text.splitlines()
    orig = [kv(x) for x in lines if x.startswith("FILTERDDP_ORIGINAL_PRIMAL ")]
    timing = [kv(x) for x in lines if x.startswith("FILTERDDP_ITER_TIMING ")]
    table = {}
    for x in lines:
        m = re.match(r"^\s*(\d+)\s+([-0-9.]+e[-+]\d+)\s+([0-9.]+e[-+]\d+)\s", x)
        if m: table[int(m[1])] = (float(m[2]), float(m[3]))
    solve = float(re.search(r"solve complete: ([0-9.]+) s", text)[1])
    setup = solve - sum(float(t["total_s"]) for t in timing)
    recs, elapsed = [], 0.0
    for i, t in enumerate(timing):
        k = int(t["iteration"])
        o = orig[i] if len(orig) == len(timing) else None
        if o is not None:
            assert (o["iteration"], o["barrier_iteration"]) == (t["iteration"], t["barrier_iteration"])
        obj, internal = (float(o["objective"]), float(o["primal_inf"])) if o else table.get(k, (float("nan"), float("nan")))
        recs.append(dict(iteration=k, objective=obj, internal=internal,
                         original=float(o["original_max"]) if o else float("nan"),
                         seconds=setup + elapsed + float(t["backward_s"])))
        elapsed += float(t["total_s"])
    return recs, solve


def ipopt_records(text, internal_log):
    orig = {}
    for line in read(internal_log).splitlines():
        m = re.match(r"^\s*(\d+)r?\s+([-0-9.]+e[-+]\d+)\s+([0-9.]+e[-+]\d+)\s", line)
        if m: orig[int(m[1])] = float(m[3])
    recs = []
    for line in text.splitlines():
        if line.startswith("IPOPT_ITER "):
            d = kv(line); k = int(d["iteration"])
            recs.append(dict(iteration=k, objective=float(d["objective"]), internal=float(d["inf_pr"]),
                             original=orig.get(k, float("nan")), seconds=float(d["elapsed_s"])))
    solve = float(re.search(r"solve_time_s=([0-9.]+)", text)[1])
    return recs, solve


runs = []                                                  # (system, T, solver, repeat, records, full_s, final objective)
for d in map(Path, sys.argv[1:]):
    for f in sorted(list(d.glob("fddp_*_T*_paired_*.log")) + list(d.glob("fddp_*_T*_paired_*.log.gz"))):
        m = re.match(r"fddp_(.+)_T(\d+)_paired_.*_(ddp)(?:_r(\d))?\.log", f.name.replace(".gz", ""))
        if not m: continue
        text = read(f); recs, full = ddp_records(text)
        runs.append((m[1], int(m[2]), "ddp", m[4] or "1", recs, full,
                     float(re.search(r"FilterDDP objective=([-0-9.eE+]+)", text)[1])))
    sysname = None
    if (d / "hardware.txt").exists():
        m = re.search(r"campaign_job\.sh (\S+)|new_systems_job\.sh (\S+)", read(d / "hardware.txt"))
        sysname = (m[1] or m[2]) if m else None
    for f in sorted(list(d.glob("ipopt_T*_ma*.log")) + list(d.glob("ipopt_T*_ma*.log.gz"))):
        name = f.name.replace(".gz", "")
        m = re.match(r"ipopt_T(\d+)_(ma57|ma97)(?:_r(\d))?\.log$", name)
        if not m: continue
        text = read(f)
        s = re.search(r"CENTRAL_IPOPT system=(\S+)", text)
        if not s: continue
        internal = f.with_name(name.replace(".log", "_internal.log"))
        recs, full = ipopt_records(text, internal)
        runs.append((s[1], int(m[1]), m[2], m[3] or "1", recs, full,
                     float(re.search(r"CENTRAL_IPOPT .* objective=([-0-9.eE+]+)", text)[1])))

refs = {(s, T): obj for (s, T, solver, rep, _, _, obj) in runs if solver == "ma97" and rep == "1"}
octrows = list(csv.DictReader(open(OCT))) if OCT.exists() else []
for r in octrows:
    refs.setdefault((r["system"], int(r["T"])), float(r["reference_objective"]))

rows = []
for (s, T, solver, rep, recs, full, _) in runs:
    if (s, T) not in refs: continue
    ref = refs[(s, T)]; thr = 1e-6 if s == "ieee123C_1ph" else 1e-5
    for gap in (0.005, 0.001):
        for metric in ("original", "internal"):
            hit = next((x for x in recs if x[metric] <= thr and abs(x["objective"] - ref) / abs(ref) <= gap), None)
            rows.append(dict(system=s, T=T, solver=solver, repeat=rep, gap=gap, metric=metric, threshold=thr,
                             reference=ref, reached=hit is not None, iteration=hit["iteration"] if hit else "",
                             seconds=hit["seconds"] if hit else "", full_seconds=full, source="nov17"))
for r in octrows:                                          # October's Ipopt rows for systems not rerun here
    if r["solver"] in ("ma57", "ma97") and not any(x[0] == r["system"] and x[1] == int(r["T"]) and x[2] == r["solver"] for x in runs):
        rows.append(dict(system=r["system"], T=int(r["T"]), solver=r["solver"], repeat=r["repeat"], gap=float(r["gap_percent"]) / 100,
                         metric=r["metric"], threshold=float(r["primal_threshold"]), reference=float(r["reference_objective"]),
                         reached=r["status"] == "reached", iteration=r["iteration"], seconds=r["seconds"],
                         full_seconds=float(r["full_seconds"]), source="oct16"))
with open(HERE / "campaign_rows.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

def med(system, T, solver, gap=0.001, metric="original"):
    v = [float(r["seconds"]) for r in rows if r["system"] == system and r["T"] == T and r["solver"] == solver
         and abs(r["gap"] - gap) < 1e-12 and r["metric"] == metric and r["reached"]]
    n = len([r for r in rows if r["system"] == system and r["T"] == T and r["solver"] == solver and abs(r["gap"] - gap) < 1e-12 and r["metric"] == metric])
    return (statistics.median(v) if v and len(v) == n else None), n

for metric in ("original", "internal"):
    print(f"\nTime to the 0.1% gap with {metric} residual <= threshold, two-repeat medians (s):")
    print(f"{'system':26s} {'T':>4s} {'DDP':>9s} {'MA57':>9s} {'MA97':>9s} {'best Ipopt / DDP':>17s}  repeats")
    for (s, T) in sorted({(r["system"], r["T"]) for r in rows if r["solver"] == "ddp"}):
        d, nd = med(s, T, "ddp", metric=metric); a, na = med(s, T, "ma57", metric=metric); b, nb = med(s, T, "ma97", metric=metric)
        best = min([x for x in (a, b) if x is not None], default=None)
        f = lambda x: f"{x:9.2f}" if x is not None else f"{'-':>9s}"
        print(f"{s:26s} {T:4d} {f(d)} {f(a)} {f(b)} {(f'{best / d:17.2f}' if best and d else f'{chr(45):>17s}')}  {nd}/{na}/{nb}")
