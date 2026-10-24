"""Strict posthoc collector for the ONE large-late compact chain experiment."""
import argparse
import csv
import hashlib
import json
import math
import statistics
import tomllib
from pathlib import Path


def rows(p):
    with p.open(newline="") as f:
        return list(csv.DictReader(f))


def main(root, out):
    root, out = Path(root), Path(out)
    assert not out.exists(), "fresh analysis directory required"
    out.mkdir(parents=True)
    r = root / "results" / "large_late"
    m = tomllib.loads((r / "metadata.toml").read_text())
    assert (root / "compact_chain.finished").is_file()
    assert m["nB"] == 1020 and m["window"] == 4 and m["N"] == 24
    assert m["accurate"] and m["inertia_ok"] and m["negative_control_rejected"]
    assert m["algorithm_dense_exports"] == 0
    assert not m["optimizer_run"] and not m["medium_extended"]
    assert m["leaves_sha256"] == "cd7922cf75df0743cf5e3fd38e2e87e4e65a861346e89bc66512cf1ffb90f416"
    assert m["input_sha256"] == "b750a212459d7972a4a72b18ef29d7ee8337c22db970a2e2c0eab5ce40dfce9c"
    assert m["saved_probe_error"] < 1e-7 and m["ordinary_error"] < 1e-7
    assert m["max_value_scaled"] < 1e-7
    assert m["direction_scaled"] < 1e-7 and m["direction_relative"] < 1e-5
    assert m["consumer_reference_scaled"] < 1e-7
    assert m["negative_scaled"] >= 1e-7
    unit = (root / "logs" / "unit.log").read_text()
    assert "COMPACT_UNIT_PASS" in unit and "COMPACT_CHAIN_UNIT_PASS" in unit
    for name in ("source-check.log", "inherited-source-check.log", "input-check.log"):
        log = (root / "logs" / name).read_text()
        assert log.strip() and all(x.endswith(": OK") for x in log.splitlines())
    # Recompute every overlay checksum from the archived actual sources.
    for line in (root / "source.sha256").read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        assert hashlib.sha256((root / name.lstrip(" *")).read_bytes()).hexdigest() == expected
    stderr = list(root.glob("slurm-*.err"))
    assert len(stderr) == 1 and not stderr[0].read_text().strip()
    with (root / "accounting.txt").open(newline="") as f:
        acct = list(csv.DictReader(f, delimiter="|"))
    assert len(acct) == 3 and all(x["State"] == "COMPLETED" and x["ExitCode"] == "0:0" for x in acct)
    assert all(x["AllocCPUS"] == "64" and x["ReqMem"] == "300Gn" for x in acct)
    assert (root / "overlay.commit").read_text().strip() == "eb4f0439"
    log = (root / "logs" / "large_late.log").read_text()
    assert "COMPACT_CHAIN_COMPLETE accurate=true" in log and "ERROR:" not in log
    hardware = (root / "hardware.txt").read_text()
    assert "Gold 6438Y+" in hardware and "NumCPUs=64" in hardware
    assert "mem=300G" in hardware or "mem=307200M" in hardware
    checks = rows(r / "repeat_checks.csv")
    assert len(checks) == 6
    assert {(int(x["repetition"]), x["mode"]) for x in checks} == {(j, k) for j in range(1, 4) for k in ("dense", "compact")}
    for x in checks:
        assert x["pass"] == "true" and float(x["scaled"]) < 1e-7 and float(x["relative"]) < 1e-5
    values = rows(r / "values.csv")
    assert len(values) == 14 and all(float(x["scaled"]) < 1e-7 for x in values)
    directions = rows(r / "directions.csv")
    assert {x["field"] for x in directions} == {"du", "dphi", "dzl", "dzu", "dx", "dp"}
    assert all(float(x["scaled"]) < 1e-7 and float(x["relative_floor_1e_10"]) < 1e-5
               for x in directions if x["mode"] != "corrupt_composed_gradient")
    assert max(float(x["scaled"]) for x in directions
               if x["mode"] == "corrupt_composed_gradient" and x["field"] in ("du", "dphi", "dzl", "dzu")) > 1e-7
    ts = rows(r / "timings.csv")
    phases = {
        "dense": {"three_merges", "downstream_directions", "cleanup_gc"},
        "compact": {"fresh_input_discovery", "compact_merges_and_close", "condition_probe", "downstream_directions", "cleanup_gc"},
    }
    totals = []
    for rep in range(1, 4):
        for mode in ("dense", "compact"):
            selected = [x for x in ts if int(x["repetition"]) == rep and x["mode"] == mode]
            assert len(selected) == len(phases[mode]) and {x["phase"] for x in selected} == phases[mode]
            expected_order = 1 if (mode == "dense") == bool(rep % 2) else 2
            assert all(int(x["order"]) == expected_order for x in selected)
            for x in selected:
                assert all(math.isfinite(float(x[k])) and float(x[k]) >= 0 for k in ("seconds", "alloc_bytes", "gc_s"))
            time = {x["phase"]: float(x["seconds"]) for x in selected}
            total = sum(time.values())
            discovery = time.get("fresh_input_discovery", 0)
            totals.append(dict(repetition=rep, mode=mode, fresh_charged_s=total,
                               same_iterate_factors_available_s=total-discovery,
                               fresh_without_cleanup_s=total-time["cleanup_gc"],
                               cached_without_cleanup_s=total-discovery-time["cleanup_gc"],
                               alloc_bytes=sum(int(x["alloc_bytes"]) for x in selected),
                               gc_s=sum(float(x["gc_s"]) for x in selected)))
    with (out / "totals.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(totals[0]));w.writeheader();w.writerows(totals)
    medians = {}
    for mode in phases:
        medians[mode] = {k: statistics.median(x[k] for x in totals if x["mode"] == mode)
                         for k in totals[0] if k not in ("repetition", "mode")}
    phase_summary = []
    for mode, names in phases.items():
        for name in sorted(names):
            xs = [x for x in ts if x["mode"] == mode and x["phase"] == name]
            phase_summary.append(dict(mode=mode, phase=name, seconds=statistics.median(float(x["seconds"]) for x in xs),
                                      alloc_bytes=statistics.median(int(x["alloc_bytes"]) for x in xs),
                                      gc_s=statistics.median(float(x["gc_s"]) for x in xs)))
    with (out / "phases.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(phase_summary[0]));w.writeheader();w.writerows(phase_summary)
    rank_summary = []
    rs = rows(r / "ranks.csv")
    for label in sorted({x["representation"] for x in rs}):
        xs = [x for x in rs if x["representation"] == label]
        rootrow = next(x for x in xs if x["depth"] == "0")
        rank_summary.append(dict(representation=label, rows=int(rootrow["node_rows"]),
                                 coefficient_bytes=8*int(rootrow["entries"]),
                                 max_stored_factor_width_including_dense_identity=max(int(x[k]) for x in xs for k in ("lr_rank", "rl_rank")),
                                 max_compressed_rank=max([int(x[k]) for x in xs for k, flag in (("lr_rank", "lr_dense"), ("rl_rank", "rl_dense")) if x[flag] == "false"] or [0]),
                                 dense_cross_nodes=sum(x["lr_dense"] == "true" for x in xs)))
    with (out / "rank_growth.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rank_summary[0]));w.writeheader();w.writerows(rank_summary)
    summary = dict(complete=True, points=1, repetitions=3, metadata=m, medians=medians,
                   input_coefficient_bytes_including_terminal_v=m["input_coefficient_bytes"]+8*m["nB"],
                   coefficient_accounting_note="raw metadata omitted terminal.v only; collector adds its 8*nB bytes, without changing raw evidence",
                   fresh_compact_over_dense=medians["compact"]["fresh_charged_s"]/medians["dense"]["fresh_charged_s"],
                   conditional_cached_compact_over_dense=medians["compact"]["same_iterate_factors_available_s"]/medians["dense"]["fresh_charged_s"],
                   production_defaults_changed=False, full_optimizer_tested=False,
                   actual_resident_memory_savings_demonstrated=False)
    (out / "validation.json").write_text(json.dumps(summary, indent=2)+"\n")
    with (out / "evidence.sha256").open("w") as f:
        for p in sorted(root.rglob("*")):
            if p.is_file():
                f.write(hashlib.sha256(p.read_bytes()).hexdigest()+"  "+p.relative_to(root).as_posix()+"\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "metadata"}, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser();p.add_argument("root");p.add_argument("out")
    a = p.parse_args();main(a.root, a.out)
