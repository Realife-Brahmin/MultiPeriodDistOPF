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
parser.add_argument("--T", type=int, default=96)
parser.add_argument("--partial", action="store_true")
args = parser.parse_args()
root, out, horizon = args.root, args.out, args.T
out.mkdir(parents=True, exist_ok=True)
arms = ["c1", "c2", "c3", "c4"] if horizon == 96 else ["c1", (root / "selected_candidate.txt").read_text().strip()]
rows, traces, comparisons, missing = [], {}, [], []


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
        for bad in ("WARMUP_STREAM_MISMATCH", "FILTERDDP_PROCS_FAIL", "ERROR: LoadError"):
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
        workers, threads = (31, 2) if arm == "c4" else (63, 1)
        for tag in ("after-setup", "after-solve"):
            audit = (folder / f"worker-audit-{tag}.txt").read_text().splitlines()
            assert len(audit) == workers
            for line in audit:
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
        expected_merges = (1392 if horizon == 288 else 368) if arm == "c1" else (912 if horizon == 288 else 256)
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
        rows.append(row)
        traces[(arm, rep)] = quality
        if arm != "c1":
            path = root / f"comparisons_T{horizon}/{arm}_r{rep}.toml"
            if path.exists():
                comparisons.append(dict(arm=arm, repeat=rep, **read_comparison(path)))
            else:
                missing.append(str(path.relative_to(root)))

for rep in (1, 2, 3):
    reference = next((r for r in rows if r["arm"] == "c1" and r["repeat"] == rep), None)
    if reference is None:
        continue
    for row in (r for r in rows if r["repeat"] == rep):
        row["objective_gap_to_reference"] = abs(row["objective"] - reference["objective"]) / max(1, abs(reference["objective"]))
        assert row["objective_gap_to_reference"] <= .001
        assert row["x0_sha256"] == reference["x0_sha256"] and row["ubar_sha256"] == reference["ubar_sha256"]

trajectory_differences = []
for (arm, rep), trace in traces.items():
    if arm == "c1" or ("c1", rep) not in traces:
        continue
    ref = {(p["iteration"], p["barrier_iteration"]): p for p in traces[("c1", rep)]}
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
    for key in ("full_s", "A_s", "B_s", "C_s", "accepted_steps", "solver_iterations", "backward_passes", "mu_updates", "backtracks"):
        result[key] = statistics.median(r[key] for r in selected)
    result["full_min_s"] = min(r["full_s"] for r in selected)
    result["full_max_s"] = max(r["full_s"] for r in selected)
    result["strict_solution_match"] = arm == "c1" or (sum(c["arm"] == arm for c in comparisons) == len(selected) and all(c["strict_existing_match_gate"] for c in comparisons if c["arm"] == arm))
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
