"""Strict collector: partial direct generation on exactly one immutable window."""
import argparse
import csv
import hashlib
import json
import math
import statistics
import tomllib
from pathlib import Path


def readcsv(p, delimiter=","):
    with p.open(newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def main(root, out):
    root, out = Path(root), Path(out)
    assert not out.exists();out.mkdir(parents=True)
    r = root / "results" / "large_late"
    m = tomllib.loads((r / "metadata.toml").read_text())
    assert (root / "direct_initial.finished").is_file()
    assert m["accurate"] and m["poisoned_Sh0_independence_pass"] and m["negative_control_rejected"]
    assert (m["nB"], m["N"], m["window"]) == (1020, 24, 4)
    assert m["available_stage_records"] == [1, 12, 24] and m["direct_generated_stages"] == [1]
    assert m["unavailable_stage_records"] == [2, 3, 4] and not m["true_future_direct"]
    assert not any(m[k] for k in ("full_direct_pipeline_tested", "optimizer_run", "medium_extended", "new_capture", "reuse_across_iterates", "initial_dense_SVD_for_direct_stage"))
    assert m["initial_full_boundary_matrices_constructed_by_direct_path"] == 0 and m["algorithm_dense_exports"] == 0
    assert all(m[k] for k in ("exact_timed_helpers_warmed", "condition_probe_warmed", "compilation_recorded_every_phase"))
    for key in ("ordinary_error", "saved_probe_error", "consumer_reference_scaled", "value_max_scaled", "direction_scaled"):
        assert m[key] < 1e-7
    assert m["direction_relative"] < 1e-5 and m["negative_scaled"] > 1e-7
    assert m["leaves_sha256"] == "cd7922cf75df0743cf5e3fd38e2e87e4e65a861346e89bc66512cf1ffb90f416"
    assert m["input_sha256"] == "b750a212459d7972a4a72b18ef29d7ee8337c22db970a2e2c0eab5ce40dfce9c"
    for line in (root / "source.sha256").read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        assert hashlib.sha256((root / name.lstrip(" *")).read_bytes()).hexdigest() == expected
    for name in ("source-check.log", "inherited-source-check.log", "input-check.log"):
        log = (root / "logs" / name).read_text()
        assert log.strip() and all(x.endswith(": OK") for x in log.splitlines())
    units = (root / "logs" / "unit.log").read_text()
    assert all(x in units for x in ("COMPACT_UNIT_PASS", "COMPACT_CHAIN_UNIT_PASS", "DIRECT_INITIAL_UNIT_PASS", "DIRECT_VALUE_LOG_REGRESSION_PASS"))
    errors = list(root.glob("slurm-*.err"));assert len(errors) == 1 and not errors[0].read_text().strip()
    log = (root / "logs" / "large_late.log").read_text()
    assert "DIRECT_INITIAL_COMPLETE accurate=true full_direct_pipeline=false" in log and "ERROR:" not in log
    acct = readcsv(root / "accounting.txt", "|")
    assert len(acct) == 3 and all(x["State"] == "COMPLETED" and x["ExitCode"] == "0:0" for x in acct)
    assert all(x["AllocCPUS"] == "64" and x["ReqMem"] == "300Gn" for x in acct)
    hardware = (root / "hardware.txt").read_text();assert "Gold 6438Y+" in hardware and "NumCPUs=64" in hardware
    values = readcsv(r / "values.csv");assert len(values) == 22
    assert all(math.isfinite(float(x["scaled"])) and float(x["scaled"]) < 1e-7 for x in values if x["phase"] != "raw_skew")
    ds = readcsv(r / "directions.csv")
    assert {x["field"] for x in ds} == {"du", "dphi", "dzl", "dzu", "dx", "dp"}
    assert all(float(x["scaled"]) < 1e-7 and float(x["relative_floor_1e_10"]) < 1e-5 for x in ds if x["mode"] != "negative_composed_gradient")
    assert max(float(x["scaled"]) for x in ds if x["mode"] == "negative_composed_gradient" and x["field"] == "du") > 1e-7
    cs = readcsv(r / "repeat_checks.csv");assert len(cs) == 6
    assert {(int(x["repetition"]), x["mode"]) for x in cs} == {(j, k) for j in (1, 2, 3) for k in ("dense", "partial_direct")}
    for x in cs:
        assert x["pass"] == "true" and float(x["direction_scaled"]) < 1e-7 and float(x["direction_relative"]) < 1e-5 and float(x["prefix_scaled"]) < 1e-7
    ts = readcsv(r / "timings.csv")
    phases = {
        "dense": {"dense_leaf1_generation", "three_merges", "closed1_actions", "directions", "cleanup_gc"},
        "partial_direct": {"direct_leaf1_generation", "unavailable_leaf2to4_discovery", "true_future_discovery", "compact_updates", "condition_probe", "closed1_actions", "directions", "cleanup_gc"},
    }
    totals = []
    for rep in (1, 2, 3):
        for mode in phases:
            xs = [x for x in ts if int(x["repetition"]) == rep and x["mode"] == mode]
            assert len(xs) == len(phases[mode]) and {x["phase"] for x in xs} == phases[mode]
            pos = 1 if (mode == "dense") == bool(rep % 2) else 2
            assert all(int(x["order"]) == pos for x in xs)
            assert all(math.isfinite(float(x[k])) and float(x[k]) >= 0 for x in xs for k in ("seconds", "alloc_bytes", "gc_s"))
            assert all(float(x["compile_s"]) == 0 and float(x["recompile_s"]) == 0 for x in xs), "mixed compilation/execution is not a warmed benchmark"
            by = {x["phase"]: float(x["seconds"]) for x in xs}
            total = sum(by.values());discovery = by.get("unavailable_leaf2to4_discovery", 0)+by.get("true_future_discovery", 0)
            setup = discovery+by.get("direct_leaf1_generation", by.get("dense_leaf1_generation", 0))
            totals.append(dict(repetition=rep, mode=mode, fresh_charged_s=total, setup_s=setup,
                               recurring_s=total-setup, without_cleanup_s=total-by["cleanup_gc"],
                               alloc_bytes=sum(int(x["alloc_bytes"]) for x in xs), gc_s=sum(float(x["gc_s"]) for x in xs)))
    medians = {mode: {k: statistics.median(x[k] for x in totals if x["mode"] == mode) for k in totals[0] if k not in ("repetition", "mode")} for mode in phases}
    with (out / "totals.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(totals[0]));w.writeheader();w.writerows(totals)
    ps = []
    for mode in phases:
        for phase in sorted(phases[mode]):
            xs = [x for x in ts if x["mode"] == mode and x["phase"] == phase]
            ps.append(dict(mode=mode, phase=phase, seconds=statistics.median(float(x["seconds"]) for x in xs),
                           alloc_bytes=statistics.median(int(x["alloc_bytes"]) for x in xs), gc_s=statistics.median(float(x["gc_s"]) for x in xs)))
    with (out / "phases.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(ps[0]));w.writeheader();w.writerows(ps)
    summary = dict(complete=True, points=1, repetitions=3, direct_stage_coverage="period1 only; chain2:4 and true future still discovered",
                   metadata=m, medians=medians, partial_direct_over_dense=medians["partial_direct"]["fresh_charged_s"]/medians["dense"]["fresh_charged_s"],
                   full_direct_pipeline_tested=False, optimizer_tested=False, actual_resident_memory_saving_demonstrated=False)
    (out / "validation.json").write_text(json.dumps(summary, indent=2)+"\n")
    with (out / "evidence.sha256").open("w") as f:
        for p in sorted(root.rglob("*")):
            if p.is_file():
                f.write(hashlib.sha256(p.read_bytes()).hexdigest()+"  "+p.relative_to(root).as_posix()+"\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "metadata"}, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser();p.add_argument("root");p.add_argument("out");a = p.parse_args();main(a.root, a.out)
