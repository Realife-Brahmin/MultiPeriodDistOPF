"""Strict collection of the one-operation compact Schur experiment."""
import argparse
import csv
import hashlib
import json
import math
import statistics
import tomllib
from pathlib import Path


def rows(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write(path, data):
    if data:
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)


def analyze(root, out):
    assert (root / "compact.finished").is_file()
    assert "COMPACT_UNIT_PASS" in (root / "logs/unit.log").read_text()
    log = (root / "logs/large_late.log").read_text()
    assert "COMPACT_REFERENCE_PASS" in log and "COMPACT_COMPLETE" in log
    assert "ERROR:" not in log
    err = list((root / "logs").glob("slurm-*.err"))
    assert len(err) == 1 and not err[0].read_text().strip()
    for file in ("source-check.log", "inherited-source-check.log", "input-check.log"):
        lines = (root / "logs" / file).read_text().splitlines()
        assert lines and all(x.endswith(": OK") for x in lines)
    for line in (root / "source.sha256").read_text().splitlines():
        sha, name = line.split(maxsplit=1)
        assert hashlib.sha256((root / name.rsplit("/", 1)[-1]).read_bytes()).hexdigest() == sha
    hw = (root / "hardware.txt").read_text()
    assert "Gold 6438Y+" in hw and "NumCPUs=64" in hw and "mem=300G" in hw
    case = root / "results/large_late"
    m = tomllib.loads((case / "metadata.toml").read_text())
    assert m["nB"] == 1020
    assert m["leaves_sha256"] == "cd7922cf75df0743cf5e3fd38e2e87e4e65a861346e89bc66512cf1ffb90f416"
    assert m["input_sha256"] == "b750a212459d7972a4a72b18ef29d7ee8337c22db970a2e2c0eab5ce40dfce9c"
    out.mkdir(parents=True, exist_ok=True)
    accepted = m["outcome"] == "compact_accuracy_pass"
    totals, phase_medians, repeats = [], [], []
    if m["outcome"] == "compact_rejected_dense_fallback":
        assert m["reason"] and m["fallback_s"] > 0 and not m["optimizer_run"]
    else:
        assert m["compact_component"] and not m["compact_full_scan"]
        assert not m["full_optimizer_run"] and not m["medium_extended"]
        assert m["legacy_export_required"]
        assert m["tested_merge_left_period"] == 3
        assert m["direction_audit_negative_control_rejected"]
        assert m["negative_control_full_direction_scaled"] > 1e-7
        assert m["dense_actions_exact_path_warmed"]
        assert m["saved_probe_error"] < 1e-7
        assert m["julia_threads"] == m["blas_threads"] == 1
        assert m["stats"]["max_local_dense_dim"] <= 64
        directions = rows(case / "directions.csv")
        values = rows(case / "value_errors.csv")
        assert {"Maa", "Mab", "Mbb", "ma", "mb", "V", "v"} <= {r["field"] for r in values}
        hybrid = [r for r in directions if r["mode"] == "compact_hybrid"]
        assert {"dx", "dp", "du", "dphi", "dzl", "dzu"} <= {r["field"] for r in hybrid}
        assert all(math.isfinite(float(r["scaled"])) for r in values + directions)
        negative = [r for r in directions if r["mode"] == "sensitivity_negative_control" and r["field"] in
                    {"du", "dphi", "dzl", "dzu"}]
        assert negative and max(float(r["scaled"]) for r in negative) > 1e-7
        actual = (max(float(r["scaled"]) for r in values) < 1e-7
                  and max(float(r["scaled"]) for r in hybrid) < 1e-7
                  and max(float(r["relative_floor_1e_10"]) for r in hybrid) < 1e-5
                  and m["action_max_scaled"] < 1e-7 and m["solve_residual_scaled"] < 1e-7)
        assert accepted == actual
        timings = rows(case / "timings.csv")
        checks = rows(case / "repeat_checks.csv")
        if not accepted:
            assert not timings and not checks
        else:
            assert len(checks) == 3 and all(r["pass"] == "true" for r in checks)
            assert len(timings) == 36
            for check in checks:
                rep = int(check["repetition"])
                assert rep in (1, 2, 3)
                assert float(check["values_scaled"]) < 1e-7
                assert float(check["directions_scaled"]) < 1e-7
                assert float(check["directions_relative"]) < 1e-5
                rr = [r for r in directions if r["mode"] == f"repeat{rep}"]
                assert rr and {"dx", "dp", "du", "dphi", "dzl", "dzu"} <= {r["field"] for r in rr}
                assert max(float(r["scaled"]) for r in rr) < 1e-7
                assert max(float(r["relative_floor_1e_10"]) for r in rr) < 1e-5
            for rep in (1, 2, 3):
                order = ("dense", "compact") if rep % 2 else ("compact", "dense")
                for pos, mode in enumerate(order, 1):
                    rs = [r for r in timings if int(r["repetition"]) == rep and r["mode"] == mode]
                    expected = ({"workspace", "merge", "actions", "legacy_downstream", "cleanup_gc"}
                                if mode == "dense" else
                                {"fresh_input_factor_discovery", "merge_factor_updates", "condition_probe", "actions",
                                 "legacy_dense_export", "legacy_downstream", "cleanup_gc"})
                    assert len(rs) == len(expected) and {r["phase"] for r in rs} == expected
                    assert all(int(r["order"]) == pos for r in rs)
                    assert all(math.isfinite(float(r["seconds"])) and float(r["seconds"]) >= 0 for r in rs)
                    total = sum(float(r["seconds"]) for r in rs)
                    recurring = sum(float(r["seconds"]) for r in rs if r["phase"] in
                                    ({"merge", "actions"} if mode == "dense" else
                                     {"merge_factor_updates", "condition_probe", "actions"}))
                    cached_pipeline = sum(float(r["seconds"]) for r in rs if r["phase"] not in
                                          {"fresh_input_factor_discovery", "workspace"})
                    totals.append(dict(repetition=rep, mode=mode, all_phases_s=total,
                        observed_cold_topology_plus_warm_phases_s_not_benchmark=total+(m["topology_setup_s"] if mode == "compact" else 0),
                        execution_without_explicit_cleanup_s=sum(float(r["seconds"]) for r in rs if r["phase"] != "cleanup_gc"),
                        recurring_component_s=recurring,
                        cached_inputs_legacy_pipeline_s=cached_pipeline,
                        alloc_bytes=sum(int(r["alloc_bytes"]) for r in rs),
                        inside_phase_gc_s=sum(float(r["gc_s"]) for r in rs)))
                    repeats.extend(dict(r) for r in rs)
            for mode in ("dense", "compact"):
                for phase in sorted({r["phase"] for r in timings if r["mode"] == mode}):
                    rs = [r for r in timings if r["mode"] == mode and r["phase"] == phase]
                    phase_medians.append(dict(mode=mode, phase=phase,
                        median_s=statistics.median(float(r["seconds"]) for r in rs),
                        min_s=min(float(r["seconds"]) for r in rs),max_s=max(float(r["seconds"]) for r in rs),
                        median_alloc_bytes=statistics.median(int(r["alloc_bytes"]) for r in rs)))
    write(out / "totals.csv", totals)
    write(out / "phases.csv", phase_medians)
    summary = dict(complete=True, outcome=m["outcome"], accuracy_pass=accepted,
                   optimizer_tested=False, fully_compact_scan=False, isolated_resident_memory_comparison=False,
                   phase_b_timing=False, repeats_per_mode=3 if accepted else 0)
    if accepted:
        dense = [r["all_phases_s"] for r in totals if r["mode"] == "dense"]
        compact = [r["all_phases_s"] for r in totals if r["mode"] == "compact"]
        summary.update(dense_all_phases_median_s=statistics.median(dense),
                       compact_all_phases_median_s=statistics.median(compact),
                       compact_over_dense=statistics.median(compact)/statistics.median(dense))
    (out / "validation.json").write_text(json.dumps(summary, indent=2)+"\n")
    with (out / "evidence.sha256").open("w") as f:
        for p in sorted(root.rglob("*")):
            if p.is_file() and "julia-depot" not in p.parts:
                f.write(hashlib.sha256(p.read_bytes()).hexdigest()+"  "+p.relative_to(root).as_posix()+"\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("root", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    analyze(a.root, a.out)
