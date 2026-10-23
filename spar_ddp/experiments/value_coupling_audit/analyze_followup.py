"""Strict posthoc collection of the three-point October23 frozen follow-up."""
import argparse
import csv
import hashlib
import json
import math
import statistics
import tomllib
from pathlib import Path

CASES = ("ieee2522C_1ph_late", "large10kC_1ph_late", "ieee2522C_1ph_early")


def read_csv(p):
    with p.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(p, rows):
    if not rows:
        return
    with p.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def analyze(root, out, partial=False):
    out.mkdir(parents=True, exist_ok=True)
    if not partial:
        assert (root / "followup.finished").is_file()
        errs = list((root / "logs").glob("slurm-*.err"))
        assert len(errs) == 1 and not errs[0].read_text().strip()
        assert "BLOCK_CONTROL_UNIT_PASS" in (root / "logs/unit.log").read_text()
        for name in ("source-check.log", "inherited-source-check.log", "input-check.log"):
            lines = (root / "logs" / name).read_text().splitlines()
            assert lines and all(x.endswith(": OK") for x in lines)
        for line in (root / "source.sha256").read_text().splitlines():
            expected, source = line.split(maxsplit=1)
            archived = root / source.rsplit("/", 1)[-1]
            assert hashlib.sha256(archived.read_bytes()).hexdigest() == expected
        assert "DIRECTION_MEASUREMENT_REGRESSION_PASS" in (root / "logs/unit.log").read_text()
        if (root / "followup_equivalence.jl").is_file():
            unit = (root / "logs/unit.log").read_text()
            assert "FIXED_NODE_EQUIVALENCE_PASS" in unit
            assert "FIXED_NODE_REAL_LEAF_EQUIVALENCE_PASS" in unit
            legacy = (root / "logs/legacy-source-check.log").read_text().splitlines()
            assert legacy and all(x.endswith(": OK") for x in legacy)
        hardware = (root / "hardware.txt").read_text()
        assert "Gold 6438Y+" in hardware and "NumCPUs=64" in hardware
        assert "mem=300G" in hardware
    summary, warmed, components, missing = [], [], [], []
    field_errors, rank_stats = [], []
    for case in CASES:
        path = root / "results" / case
        if not (path / "metadata.toml").is_file():
            missing.append(case)
            continue
        meta = tomllib.loads((path / "metadata.toml").read_text())
        assert meta["probe_error"] < 1e-7 and meta["ordinary_error"] < 1e-7
        assert meta["direction_scaled_gate"] == 1e-7
        assert meta["direction_relative_gate"] == 1e-5 and meta["relative_floor"] == 1e-10
        assert meta["blas_threads"] == meta["julia_threads"] == 1
        assert not meta["optimizer_run"] and not meta["compact_kernel"]
        assert not meta["resident_memory_saving_demonstrated"]
        assert not meta["component_io_in_timing"]
        assert "FOLLOWUP_COMPLETE" in (root / "logs" / (case + ".log")).read_text()
        dirs = read_csv(path / "directions.csv")
        costs = read_csv(path / "costs.csv")
        inert = read_csv(path / "inertia.csv")
        ranks = read_csv(path / "ranks.csv")
        reps = read_csv(path / "warmed.csv")
        assert dirs and inert and ranks and costs
        assert all(math.isfinite(float(x[k])) for x in dirs for k in ("absolute", "reference", "scaled", "relative_floor_1e_10"))
        for trial in meta["trials"]:
            label = trial["mode"]
            rows = [r for r in dirs if r["mode"] == label]
            assert rows and {"dx", "dp", "du", "dphi", "dzl", "dzu"} <= {r["field"] for r in rows}
            scaled = max(float(r["scaled"]) for r in rows)
            rel = max(float(r["relative_floor_1e_10"]) for r in rows)
            assert scaled == trial["scaled"] and rel == trial["relative"]
            assert trial["pass"] == (scaled < 1e-7 and rel < 1e-5)
            ir = [r for r in inert if r["mode"] == label]
            assert ir and trial["inertia_ok"] == all(r["same"] == "true" for r in ir)
            summary.append(dict(case=case, mode=label, scaled=scaled, relative_floor_1e_10=rel,
                                direction_pass=trial["pass"], inertia_pass=trial["inertia_ok"],
                                payload_ratio=trial["payload_ratio"], rank_max=trial["rank_max"],
                                raw_fallbacks=trial["raw_fallbacks"], chosen=label == meta["chosen"]))
            for field in sorted({r["field"] for r in rows}):
                fr = [r for r in rows if r["field"] == field]
                field_errors.append(dict(case=case, mode=label, field=field,
                    max_absolute=max(float(r["absolute"]) for r in fr),
                    max_scaled=max(float(r["scaled"]) for r in fr),
                    max_relative_floor_1e_10=max(float(r["relative_floor_1e_10"]) for r in fr)))
            for depth in sorted({int(r["scan_depth"]) for r in ranks if r["mode"] == label}):
                rr = [r for r in ranks if r["mode"] == label and int(r["scan_depth"]) == depth]
                rank_stats.append(dict(case=case, mode=label, scan_depth=depth,
                    cross_blocks=len(rr), dense_blocks=sum(r["dense"] == "true" for r in rr),
                    rank_max=max(int(r["rank"]) for r in rr)))
        if meta["chosen"] == "none_raw_fallback":
            assert not reps
            continue
        chosen = next(t for t in meta["trials"] if t["mode"] == meta["chosen"])
        assert chosen["pass"] and chosen["inertia_ok"] and chosen["payload_ratio"] < .99
        assert len(reps) == 6
        for repeat in (1, 2, 3):
            rs = [r for r in reps if int(r["repetition"]) == repeat]
            expected = ["dense", "candidate"] if repeat % 2 else ["candidate", "dense"]
            assert [r["mode"] for r in rs] == expected
            cr = [r for r in costs if r["mode"] == f"repeat{repeat}"]
            assert cr
            ct = sum(float(r["seconds"]) for r in cr)
            candidate = next(r for r in rs if r["mode"] == "candidate")
            assert ct <= float(candidate["total_s"]) * 1.01
            for component in ("condition_control", "construction", "materialization", "combination"):
                rr = [r for r in cr if r["component"] == component]
                assert rr
                components.append(dict(case=case, repetition=repeat, component=component,
                    seconds=sum(float(r["seconds"]) for r in rr),
                    alloc_bytes=sum(int(r["alloc_bytes"]) for r in rr), gc_s=sum(float(r["gc_s"]) for r in rr)))
            for r in rs:
                assert float(r["scaled"]) < 1e-7 and float(r["relative"]) < 1e-5
                assert float(r["total_s"]) > 0 and int(r["alloc_bytes"]) >= 0
                warmed.append(dict(case=case, **r))
    if not partial:
        assert not missing
    medians = []
    for case in CASES:
        dense = [float(r["total_s"]) for r in warmed if r["case"] == case and r["mode"] == "dense"]
        candidate = [float(r["total_s"]) for r in warmed if r["case"] == case and r["mode"] == "candidate"]
        if dense:
            medians.append(dict(case=case, dense_median_s=statistics.median(dense),
                candidate_median_s=statistics.median(candidate),
                candidate_over_dense=statistics.median(candidate)/statistics.median(dense),
                candidate_min_s=min(candidate), candidate_max_s=max(candidate)))
    write_csv(out / "accuracy.csv", summary)
    write_csv(out / "warmed.csv", warmed)
    write_csv(out / "components.csv", components)
    write_csv(out / "medians.csv", medians)
    write_csv(out / "field_errors.csv", field_errors)
    write_csv(out / "rank_stats.csv", rank_stats)
    validation = dict(complete=not partial and not missing, missing=missing, points=len(CASES)-len(missing),
        warmed_points=len(medians), optimizer_tested=False, compact_kernel=False,
        resident_memory_saving_demonstrated=False, phase_b_timing=False,
        coefficient_storage_only=True, offline_direction_oracle=True,
        mixed_cold_accuracy_times_excluded=True, minimum_repetitions=3)
    (out / "validation.json").write_text(json.dumps(validation, indent=2)+"\n")
    files = [p for p in root.rglob("*") if p.is_file() and "julia-depot" not in p.parts]
    with (out / "evidence.sha256").open("w") as f:
        for p in sorted(files):
            f.write(hashlib.sha256(p.read_bytes()).hexdigest()+"  "+str(p.relative_to(root)).replace("\\", "/")+"\n")
    print(json.dumps(validation))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("root", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--partial", action="store_true")
    a = p.parse_args()
    analyze(a.root, a.out, a.partial)
