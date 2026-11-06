# Adapted from dense_merge_fullsolve_nov05/collect.py; identical scientific gates.
import argparse
import csv
import hashlib
import json
import math
import re
import statistics
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("root", type=Path)
parser.add_argument("out", type=Path)
parser.add_argument("--partial", action="store_true")
parser.add_argument("--horizon",type=int,choices=(24,96,144,288),required=True)
args = parser.parse_args()
root, out = args.root, args.out
horizon, budget = args.horizon,128
reference_arm = "base"
configs = {arm: (24 if horizon==24 else 48,2) for arm in ("base","cache")}
out.mkdir(parents=True, exist_ok=True)
arms = list(configs)
rows, traces, comparisons, missing = [], {}, [], []

if not args.partial:
    jobs = list((root/"jobs").iterdir())
    assert len(jobs) == 1
    job = jobs[0]
    hardware = (job/"hardware.txt").read_text()
    assert "6760P" in hardware and "NumCPUs=128" in hardware
    stderr = (job/"stderr.log").read_text()
    unexpected = [line for line in stderr.splitlines() if line.strip() and not re.fullmatch(r"\s*\d+ dependenc(?:y|ies) successfully precompiled in \d+ seconds?(?:\. \d+ already precompiled\.)?\s*",line)]
    assert not unexpected, unexpected
    assert re.search(r"(?m)^"+re.escape(job.name)+r"\|COMPLETED\|[^|]+\|0:0\|",(root/"accounting.txt").read_text())
    for path in [root/"inherited-source-check.log",root/"inherited-input-check.log",root/"source-check-final.log",root/"harness-check-final.log",job/"source-check.log",job/"source-check-after.log",job/"harness-check.log",job/"harness-check-after.log"]:
        text = path.read_text()
        assert text.strip() and all(line.endswith(": OK") for line in text.splitlines()), path
    assert "strict_existing_match_gate = true" in (root/f"comparisons_T{horizon}/reference_repeats.toml").read_text()


def read_csv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def read_comparison(path):
    result = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        result[key] = value == "true" if value in ("true", "false") else float(value)
    return result


for rep in (1, 2, 3):
    for arm in arms:
        folder = root / f"results/10k_1020_T{horizon}_ddp_r{arm}_{rep}"
        if not (folder / "COMPLETE").exists():
            if rep < 3:
                missing.append(folder.name)
            continue
        text = (folder / "solve.log").read_text()
        for bad in ("WARMUP_STREAM_MISMATCH", "FILTERDDP_PROCS_FAIL", "ERROR: LoadError", "BR_DEBUG_FAIL", "AssertionError"):
            assert bad not in text, (folder, bad)
        assert "WARMUP full" in text
        assert re.search(r"solve complete: [\d.]+ s, iterations=\d+, status=0\b", text)
        final = re.search(r"final residuals: primal=(\S+) dual=(\S+) complementarity=(\S+)", text)
        assert final and all(0 <= float(v) <= 1e-5 for v in final.groups())
        physical = (folder / "physical.log").read_text()
        residual = float(re.search(r"original_max=(\S+)", physical)[1])
        objective = float(re.search(r"objective=(\S+)", physical)[1])
        assert math.isfinite(objective) and 0 <= residual <= 1e-5
        dss = (folder / "opendss.log").read_text()
        assert "below limit: 0" in dss and "above limit: 0" in dss
        workers, threads = configs[arm]
        maps = list((root / "jobs").glob(f"*/{horizon}_{arm}_{rep}/cpu_groups.csv"))
        assert len(maps) == 1
        groups = [set(map(int, line.split(","))) for line in maps[0].read_text().splitlines()]
        assert len(groups) == workers + 1 and len(groups[0]) == 1
        assert all(len(group) == threads for group in groups[1:])
        allocated = [cpu for group in groups for cpu in group]
        assert len(set(allocated)) == len(allocated) <= budget
        topology = {}
        for line in (maps[0].parent / "cpu_topology.csv").read_text().splitlines():
            if not line.startswith("#"):
                cpu, core, socket = map(int, line.split(","))
                topology[cpu] = (socket, core)
        assert len({topology[cpu] for cpu in allocated}) == len(allocated)
        assert all(len({topology[cpu][0] for cpu in group}) == 1 for group in groups)
        for tag in ("after-setup", "after-solve"):
            audit = (folder / f"worker-audit-{tag}.txt").read_text().splitlines()
            assert len(audit) == workers
            for index, line in enumerate(audit):
                recorded = re.findall(r'allowed\s*=\s*"([^"]+)"', line)
                assert recorded
                for mask in recorded:
                    cpus = set()
                    for span in mask.split(","):
                        ends = list(map(int, span.split("-")))
                        cpus.update(range(ends[0], ends[-1] + 1))
                    assert cpus == groups[index+1], (folder, tag, index, cpus, groups[index+1])
                assert int(re.search(r"blas\s*=\s*(\d+)", line)[1]) == threads
                assert int(re.search(r"blas_default\s*=\s*(\d+)", line)[1]) == threads
                assert "libopenblas" in line
        summary = json.loads((folder / "timing_summary.json").read_text())
        passes = read_csv(folder / "passes.csv")
        timing = read_csv(folder / "iteration_timing.csv")
        quality = read_csv(folder / "quality_trace.csv")
        assert len(passes) == len(timing) == len(quality) > 0
        assert len(re.findall(r"^FULL_STUDY_PHASE_B_THREADS ", text, re.M)) == len(passes)
        assert all(int(p["workers"]) == workers for p in passes)
        for p in quality:
            assert all(math.isfinite(float(p[k])) for k in ("objective", "original", "dual", "complementarity"))
        expected_merges = {24:(68,52),96:(368,256),144:(624,420),288:(1392,912)}[horizon][0 if arm == "original" else 1]
        assert all(int(p["combinations"]) == expected_merges for p in passes)
        initial = re.search(r"FULL_STUDY_START x0_sha256=(\w+) ubar_sha256=(\w+)", text)
        assert initial
        row = dict(arm=arm, repeat=rep, T=horizon, workers=workers, blas_threads=threads,
                   full_s=summary["full_optimizer_s"], A_s=summary["phase_A_s"],
                   B_s=summary["scan_wall_s"], C_s=summary["phase_C_s"],
                   accepted_steps=sum(t["outcome"] == "accepted" for t in timing),
                   solver_iterations=int(re.search(r"solve complete: .*iterations=(\d+),", text)[1]),
                   backward_passes=len(passes),
                   mu_updates=sum(float(b["mu"]) < float(a["mu"]) for a, b in zip(passes, passes[1:])),
                   backtracks=sum(int(t["backtracks"]) for t in timing),
                   objective=objective, independent_original_primal=residual,
                   x0_sha256=initial[1], ubar_sha256=initial[2],
                   coordinator_maxrss_mib=summary["memory"][-1]["coordinator_maxrss_mib"],
                   worker_maxrss_mib_max=summary["memory"][-1]["worker_maxrss_mib_max"],
                   worker_maxrss_mib_sum=summary["memory"][-1]["worker_maxrss_mib_sum"])
        wall = (folder / "wall.txt").read_text()
        row["process_wall_including_setup_warmup_validation_s"] = int(re.search(r"total_process_wall_s=([0-9]+)", wall)[1])
        build = re.search(r"build complete: ([0-9.]+) s", text)
        warm = re.search(r"WARMUP full iterations=[0-9]+ time_s=([0-9.]+)", text)
        assert build and warm
        row["coordinator_model_build_s"] = float(build[1])
        row["warmup_wall_s"] = float(warm[1])
        cache_passes=[]
        for line in re.findall(r"^BR_PASS solve=1 .*",text,re.M):
            cache_passes.append(dict(re.findall(r"(\w+)=(\S+)",line)))
        assert len(cache_passes)==len(passes), (folder,len(cache_passes),len(passes))
        assert all(x["debug"]=="false" for x in cache_passes)
        expected={24:52,96:256,144:420,288:912}[horizon]
        if arm=="base":
            assert all(x["hit"]=="false" and int(x["reused"])==0 and int(x["bytes"])==0 for x in cache_passes)
        else:
            assert cache_passes[0]["hit"]=="false" and int(cache_passes[0]["populated"])==expected
            assert any(x["hit"]=="true" for x in cache_passes)
            assert all(int(x["reused"])==expected for x in cache_passes if x["hit"]=="true")
            assert all(int(x["bytes"])==expected*3*1020*1020*8 for x in cache_passes)
            assert "collected=true" in text
        row["eligible_requests"]=sum(x["requested"]=="true" for x in cache_passes)
        row["hit_passes"]=sum(x["hit"]=="true" for x in cache_passes)
        row["eligible_miss_passes"]=sum(x["requested"]=="true" and x["hit"]=="false" for x in cache_passes)
        row["fresh_passes"]=sum(x["hit"]=="false" for x in cache_passes)
        row["populated_merges"]=sum(int(x["populated"]) for x in cache_passes)
        row["reused_merges"]=sum(int(x["reused"]) for x in cache_passes)
        row["factorization_fallbacks"]=sum(int(x["fallback"]) for x in cache_passes)
        row["cache_numeric_bytes_peak"]=max(int(x["bytes"]) for x in cache_passes)
        row["worker_gc_summed_s"]=sum(float(x["gc_sum_s"]) for x in passes)
        row["worker_gc_max_per_pass_sum_s"]=sum(float(x["gc_max_s"]) for x in passes)
        refs=json.loads(Path(__file__).with_name("quality_reference.json").read_text())
        reference=float(refs[str(horizon)]["objective"])
        hits=[x for x in quality if float(x["original"])<=1e-5 and abs(float(x["objective"])-reference)/max(1,abs(reference))<=.001]
        assert hits
        row["quality_s"]=float(hits[0]["seconds"])
        row["quality_iteration"]=int(hits[0]["iteration"])
        row["reference_objective"]=reference
        row["reference_scope"]=refs[str(horizon)]["scope"]
        rows.append(row)
        traces[(arm, rep)] = quality
        if arm != reference_arm:
            path = root / f"comparisons_T{horizon}/{arm}_r{rep}.toml"
            if path.exists():
                comparisons.append(dict(arm=arm, repeat=rep, **read_comparison(path)))
            else:
                missing.append(str(path.relative_to(root)))

if rows:
    assert len({(r["x0_sha256"],r["ubar_sha256"]) for r in rows}) == 1

for rep in (1, 2, 3):
    reference = next((r for r in rows if r["arm"] == reference_arm and r["repeat"] == rep), None)
    if reference is None:
        continue
    for row in (r for r in rows if r["repeat"] == rep):
        row["objective_gap_to_reference"] = abs(row["objective"] - reference["objective"]) / max(1, abs(reference["objective"]))
        assert row["objective_gap_to_reference"] <= .001
        assert row["x0_sha256"] == reference["x0_sha256"] and row["ubar_sha256"] == reference["ubar_sha256"]

trajectory_differences = []
for (arm, rep), trace in traces.items():
    if arm == reference_arm or (reference_arm, rep) not in traces:
        continue
    ref = {(p["iteration"], p["barrier_iteration"]): p for p in traces[(reference_arm, rep)]}
    candidate = {(p["iteration"], p["barrier_iteration"]): p for p in trace}
    common = set(ref) & set(candidate)
    trajectory_differences.append(dict(arm=arm, repeat=rep, common_pass_keys=len(common),
        missing_reference_keys=len(set(ref)-set(candidate)), extra_candidate_keys=len(set(candidate)-set(ref)),
        max_objective_relative_gap=max(abs(float(candidate[k]["objective"])-float(ref[k]["objective"]))/max(1,abs(float(ref[k]["objective"]))) for k in common),
        max_original_primal_difference=max(abs(float(candidate[k]["original"])-float(ref[k]["original"])) for k in common)))

medians = []
for arm in arms:
    selected = [r for r in rows if r["arm"] == arm]
    if len(selected) < 2:
        continue
    result = dict(arm=arm, repeats=len(selected))
    for key in ("full_s", "A_s", "B_s", "C_s", "accepted_steps", "solver_iterations", "backward_passes", "mu_updates", "backtracks", "quality_s", "hit_passes", "eligible_requests", "cache_numeric_bytes_peak"):
        result[key] = statistics.median(r[key] for r in selected)
    result["full_min_s"] = min(r["full_s"] for r in selected)
    result["full_max_s"] = max(r["full_s"] for r in selected)
    result["strict_solution_match"] = arm == reference_arm or (sum(c["arm"] == arm for c in comparisons) == len(selected) and all(c["strict_existing_match_gate"] for c in comparisons if c["arm"] == arm))
    medians.append(result)

for name, data in (("per_run", rows), ("medians", medians), ("solutions", comparisons), ("trajectory_differences", trajectory_differences)):
    if data:
        fields = list(dict.fromkeys(k for row in data for k in row))
        with (out / f"{name}.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(data)
validated = not missing and (root / f"T{horizon}.finished").exists()
(out / "validation.json").write_text(json.dumps(dict(validated=validated, missing=missing,
    objective_gate=.001, physical_gate=1e-5, strict_dispatch_gate=1e-10,
    strict_equivalence_failures=[c for c in comparisons if not c["strict_existing_match_gate"]]), indent=2))
print(json.dumps(medians, indent=2))
if not args.partial:
    assert validated, missing
