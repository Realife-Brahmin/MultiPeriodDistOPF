"""Strict complete-input four-stage campaign point collector, not a solve report."""
import argparse
import csv
import hashlib
import json
import math
import statistics
import tomllib
from pathlib import Path


def rows(p, delimiter=","):
    with p.open(newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write_csv(p, xs):
    with p.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(xs[0]));w.writeheader();w.writerows(xs)


def main(root, out):
    root, out = Path(root), Path(out)
    assert not out.exists();out.mkdir(parents=True)
    p = root / "results" / "large_late"
    m = tomllib.loads((p / "metadata.toml").read_text())
    cv = tomllib.loads((root / "capture_validation.toml").read_text())
    assert cv["validated"] and cv["selected_stages"] == [1, 2, 3, 4]
    assert cv["classification"] == "bit_exact_original_Newton_state"
    assert cv["leaf_bit_exact"] and cv["nominal_bit_exact"] and cv["scalar_state_matches"]
    assert (m["nB"], m["N"], m["window"]) == (1020, 24, 4)
    assert m["generated_stages"] == [1, 2, 3, 4] and m["all_four_stage_policies_consumed"]
    assert m["boundary_construction_charged_both_modes"] and not m["true_future_given_free"]
    assert m["incoming_dense_tail_stages"] == 20 and not m["full_compact_horizon"] and not m["optimizer_run"]
    assert not m["initial_factors_reused_across_iterates"] and m["algorithm_dense_boundary_exports"] == 0
    assert m["poisoned_Sh0_independence_pass"] and m["negative_control_rejected"]
    assert m["negative_scaled"] > 1e-7
    for k in ("ordinary_error", "saved_probe_error", "independent_consumer_scaled", "dense_window_scaled"):
        assert m[k] < 1e-7
    assert m["leaves_sha256"] == "cd7922cf75df0743cf5e3fd38e2e87e4e65a861346e89bc66512cf1ffb90f416"
    assert m["input_sha256"] == "b750a212459d7972a4a72b18ef29d7ee8337c22db970a2e2c0eab5ce40dfce9c"
    assert (root / "campaign_replay.finished").is_file()
    for line in (root / "source.sha256").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        assert hashlib.sha256((root / name.lstrip(" *")).read_bytes()).hexdigest() == digest
    for name in ("source-check.log", "capture-source-check.log", "runtime-check.log", "input-check.log"):
        text = (root / "logs" / name).read_text()
        assert text.strip() and all(x.endswith(": OK") for x in text.splitlines())
    unit = (root / "logs" / "unit.log").read_text()
    assert all(x in unit for x in ("COMPACT_UNIT_PASS", "COMPACT_CHAIN_UNIT_PASS", "DIRECT_INITIAL_UNIT_PASS", "CAMPAIGN_COMPACT_UNIT_PASS"))
    log = (root / "logs" / "large_late.log").read_text()
    assert "CAMPAIGN_PROBE_COMPLETE" in log and "ERROR:" not in log
    errors = list(root.glob("slurm-*.err"));assert len(errors) == 1 and not errors[0].read_text().strip()
    acct = rows(root / "accounting.txt", "|")
    assert len(acct) == 3 and all(x["State"] == "COMPLETED" and x["ExitCode"] == "0:0" for x in acct)
    assert all(x["AllocCPUS"] == "64" and x["ReqMem"] == "300Gn" for x in acct)
    hardware = (root / "hardware.txt").read_text();assert "Gold 6438Y+" in hardware and "NumCPUs=64" in hardware
    vs = rows(p / "values.csv");assert len(vs) == 44
    ds = rows(p / "directions.csv")
    inertia = rows(p / "inertia.csv")
    assert len(inertia) == 4 and {int(x["stage"]) for x in inertia} == {1, 2, 3, 4}
    assert all(x["same"] == "true" and x["negative"] == x["reference_negative"] and x["near_zero"] == x["reference_near_zero"] for x in inertia)
    assert {int(x["stage"]) for x in ds} == {1, 2, 3, 4}
    assert {x["field"] for x in ds} == {"du", "dphi", "dzl", "dzu", "dx", "dp"}
    assert max(float(x["scaled"]) for x in ds if x["mode"] == "negative_composed_gradient" and x["field"] == "du") > 1e-7
    ts = rows(p / "timings.csv")
    totals = [];phases = []
    summary = dict(point_complete=True, accuracy_pass=m["accurate"], optimizer_tested=False,
                   full_compact_horizon=False, metadata=m, capture_validation=cv,
                   actual_resident_memory_saving_demonstrated=False)
    if m["accurate"]:
        assert m["value_max_scaled"] < 1e-7 and m["direction_scaled"] < 1e-7 and m["direction_relative"] < 1e-5 and m["inertia_ok"]
        assert all(float(x["scaled"]) < 1e-7 for x in vs if x["phase"] != "raw_skew")
        assert all(float(x["scaled"]) < 1e-7 and float(x["relative_floor_1e_10"]) < 1e-5 for x in ds if x["mode"] != "negative_composed_gradient")
        checks = rows(p / "repeat_checks.csv");assert len(checks) == 6
        assert {(int(x["repetition"]), x["mode"]) for x in checks} == {(i, j) for i in (1, 2, 3) for j in ("dense", "direct_compact")}
        assert all(x["pass"] == "true" and float(x["scaled"]) < 1e-7 and float(x["relative"]) < 1e-5 and float(x["value_actions_scaled"]) < 1e-7 for x in checks)
        expected = {
            "dense": {"incoming_dense_future_construction", "initial_dense_generation", "dense_values", "all_value_actions", "all_stage_consumers", "cleanup_gc"},
            "direct_compact": {"incoming_dense_future_construction", "initial_direct_generation", "incoming_factor_discovery", "compact_values", "condition_probes", "all_value_actions", "all_stage_consumers", "cleanup_gc"},
        }
        for rep in (1, 2, 3):
            for mode in expected:
                xs = [x for x in ts if int(x["repetition"]) == rep and x["mode"] == mode]
                assert len(xs) == len(expected[mode]) and {x["phase"] for x in xs} == expected[mode]
                pos = 1 if (mode == "dense") == bool(rep % 2) else 2
                assert all(int(x["order"]) == pos for x in xs)
                assert all(math.isfinite(float(x[k])) and float(x[k]) >= 0 for x in xs for k in ("seconds", "alloc_bytes", "gc_s"))
                assert all(float(x["compile_s"]) == 0 and float(x["recompile_s"]) == 0 for x in xs), "mixed compilation clocks are not warmed evidence"
                t = {x["phase"]: float(x["seconds"]) for x in xs}
                total = sum(t.values());tail = t["incoming_dense_future_construction"]
                setup = tail+t.get("initial_direct_generation", t.get("initial_dense_generation", 0))+t.get("incoming_factor_discovery", 0)
                totals.append(dict(repetition=rep, mode=mode, whole_charged_s=total, incoming_construction_s=tail,
                                   selected_window_fresh_s=total-tail, setup_s=setup, recurring_s=total-setup,
                                   without_cleanup_s=total-t["cleanup_gc"], alloc_bytes=sum(int(x["alloc_bytes"]) for x in xs),
                                   gc_s=sum(float(x["gc_s"]) for x in xs)))
        for mode, names in expected.items():
            for name in sorted(names):
                xs = [x for x in ts if x["mode"] == mode and x["phase"] == name]
                phases.append(dict(mode=mode, phase=name, seconds=statistics.median(float(x["seconds"]) for x in xs),
                                   alloc_bytes=statistics.median(int(x["alloc_bytes"]) for x in xs), gc_s=statistics.median(float(x["gc_s"]) for x in xs)))
        med = {mode: {k: statistics.median(x[k] for x in totals if x["mode"] == mode) for k in totals[0] if k not in ("mode", "repetition")} for mode in expected}
        paired = all(next(x["whole_charged_s"] for x in totals if x["mode"] == "direct_compact" and x["repetition"] == rep) < next(x["whole_charged_s"] for x in totals if x["mode"] == "dense" and x["repetition"] == rep) for rep in (1, 2, 3))
        w_ratio = med["direct_compact"]["selected_window_fresh_s"]/med["dense"]["selected_window_fresh_s"]
        whole_ratio = med["direct_compact"]["whole_charged_s"]/med["dense"]["whole_charged_s"]
        summary.update(fresh_cost_tested=True, zero_measured_compilation_verified=True, medians=med,
                       selected_window_ratio=w_ratio, whole_charged_ratio=whole_ratio, all_paired_whole_improve=paired,
                       advancement_cost_gate=w_ratio <= .9 and whole_ratio < 1 and paired)
        write_csv(out / "totals.csv", totals);write_csv(out / "phases.csv", phases)
    else:
        assert not ts and not rows(p / "repeat_checks.csv")
        summary.update(fresh_cost_tested=False, advancement_cost_gate=False,
                       reason="unchanged scientific gate rejected candidate; no warmed performance evidence")
    (out / "validation.json").write_text(json.dumps(summary, indent=2)+"\n")
    with (out / "evidence.sha256").open("w") as f:
        for file in sorted(root.rglob("*")):
            if file.is_file():
                f.write(hashlib.sha256(file.read_bytes()).hexdigest()+"  "+file.relative_to(root).as_posix()+"\n")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("metadata", "capture_validation")}, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser();p.add_argument("root");p.add_argument("out");a = p.parse_args();main(a.root, a.out)
