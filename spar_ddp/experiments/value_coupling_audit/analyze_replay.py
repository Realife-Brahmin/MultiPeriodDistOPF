"""Reproduce frozen-replay evidence summaries; stdlib, Python >=3.9.

No solver invocation. Clocks are serial/cold diagnostic observations, NOT
parallel wall or benchmark timings. --partial preserves incomplete outcomes.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import statistics
from pathlib import Path

MODES = ("dense", "exact", "drop", "1e-8", "1e-5")
TAGS = ("early", "barrier_before", "barrier_after", "middle", "late")
CASES = ("ieee123C_1ph_T3", "ieee2522C_1ph_T24", "large10kC_1ph_T24")


def rows(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def scalar(text, key):
    value = re.search(r"^" + re.escape(key) + r"\s*=\s*(.+)$", text, re.M).group(1)
    return json.loads(value)


def peak(data, key, default=0.0):
    return max((float(r[key]) for r in data), default=default)


def write_csv(path, data):
    if not data:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def analyze(root, out, partial=False):
    out.mkdir(parents=True, exist_ok=True)
    summaries, blocks, depths, spectral, missing, validation = [], [], [], [], [], []
    for case in CASES:
        for tag in TAGS:
            directory = root / "results" / case / tag
            log = root / "results" / case / (tag + ".log")
            text = log.read_text(errors="replace") if log.exists() else ""
            if not (directory / "metadata.toml").exists():
                missing.append({"case": case, "tag": tag, "log_exists": log.exists(),
                                "runtime_error": "ERROR:" in text})
                continue
            assert "FROZEN_REPLAY_COMPLETE" in text, (case, tag, "completion marker")
            assert "ERROR:" not in text and "UNHANDLED TASK ERROR" not in text
            meta = (directory / "metadata.toml").read_text()
            n = scalar(meta, "nB")
            N = scalar(meta, "N")
            window = scalar(meta, "window")
            assert window == (4 if case.startswith("large") else N)
            assert scalar(meta, "compact_kernel") is False
            assert scalar(meta, "timing_is_benchmark") is False
            if not partial:
                assert scalar(meta, "reference_path") == "production_inplace"
                assert scalar(meta, "ordinary_value_max_error") < 1e-7
            for key in ("leaves_sha256", "input_sha256"):
                assert re.fullmatch("[0-9a-f]{64}", scalar(meta, key))
            err = rows(directory / "errors.csv")
            direction = rows(directory / "directions.csv")
            probe = rows(directory / "probes.csv")
            cost = rows(directory / "costs.csv")
            rank = rows(directory / "ranks.csv")
            inertia = rows(directory / "inertia.csv")
            sep = rows(directory / "separators.csv")
            assert set(r["mode"] for r in err) == set(MODES)
            assert set(r["mode"] for r in direction) == set(MODES)
            assert len(probe) == (4 if window < N else 12)
            assert peak(probe, "scaled_error") < 1e-7
            for dataset in (err, direction):
                for row in dataset:
                    for key in ("absolute_fro", "reference_fro", "unit_floor_scaled", "relative_fro"):
                        assert math.isfinite(float(row[key])), (case, tag, row)
            references = {(r["depth"], r["lo"], r["hi"]): r
                          for r in inertia if r["mode"] == "dense"}
            identity = {"case": case, "tag": tag, "iteration": scalar(meta, "iteration"),
                        "barrier_iteration": scalar(meta, "barrier_iteration"),
                        "mu": scalar(meta, "mu"), "nB": n, "window": window}
            # Numerical rank at two relative singular-value cutoffs is only
            # descriptive, NOT proof of exact algebraic rank or a fixed bound.
            sg = {}
            for r in rows(directory / "spectra.csv"):
                key = tuple(r[k] for k in ("scan_depth", "lo", "hi", "hierarchy_depth", "node"))
                sg.setdefault(key, []).append(float(r["singular_value"]))
            for key, values in sg.items():
                top = max(values)
                spectral.append(dict(identity, scan_depth=int(key[0]), lo=int(key[1]), hi=int(key[2]),
                    hierarchy_depth=int(key[3]), node=int(key[4]), singular_count=len(values),
                    sigma_max=top, sigma_min=min(values),
                    numerical_rank_rel_1e12=sum(x > top * 1e-12 for x in values) if top else 0,
                    numerical_rank_rel_1e8=sum(x > top * 1e-8 for x in values) if top else 0))
            for mode in MODES:
                e = [r for r in err if r["mode"] == mode]
                d = [r for r in direction if r["mode"] == mode]
                c = [r for r in cost if r["mode"] == mode]
                rnk = [r for r in rank if r["mode"] == mode]
                ii = [r for r in inertia if r["mode"] == mode]
                assert len(ii) == len(references)
                assert {r["block"] for r in e} == {"Maa", "Mab", "Mbb", "ma", "mb", "V", "v"}
                assert {r["field"] for r in d} == {"dx", "dp", "du", "dphi", "dzl", "dzu"}
                assert len(d) == 2 * window + 8 * (1 if window < N else 3)
                if mode in ("dense", "exact"):
                    assert peak(e, "unit_floor_scaled") < 1e-7
                if mode in ("1e-8", "1e-5"):
                    assert peak(c, "reconstruction_relative_fro") <= float(mode) + 1e-12
                changed = 0
                for r in ii:
                    ref = references[r["depth"], r["lo"], r["hi"]]
                    changed += (r["negative"], r["near_zero"]) != (ref["negative"], ref["near_zero"])
                # Exclude already-closed dense PintVal records from the matrix
                # representation storage summary; their storage is unchanged.
                cc = [r for r in c if int(r["dense_three_blocks_bytes"]) == 8 * (3*n*n + 2*n)]
                ratios = [float(r["payload_bytes"]) / float(r["dense_three_blocks_bytes"]) for r in cc]
                full_dirs = [r for r in d if r["field"] in ("du", "dphi", "dzl", "dzu")]
                summaries.append(dict(identity, mode=mode,
                    probe_gate_max=peak(probe, "scaled_error"),
                    block_max_scaled=peak(e, "unit_floor_scaled"),
                    direction_max_scaled=peak(d, "unit_floor_scaled"),
                    full_direction_max_scaled=peak(full_dirs, "unit_floor_scaled"),
                    direction_max_relative=peak(d, "relative_fro"),
                    inertia_changes=changed, negative_K=max(int(r["negative"]) for r in ii),
                    K_condition_max=peak(ii, "abs_condition"),
                    leaf_and_interval_payload_ratio_median=statistics.median(ratios),
                    payload_ratio_max=max(ratios),
                    rank_max=max((int(r["rank"]) for r in rnk), default=0),
                    dense_fallback_fraction=(sum(r["dense_fallback"] == "true" for r in rnk) / len(rnk)) if rnk else 0,
                    serial_combine_s=sum(float(r["combine_s"]) for r in c),
                    construction_s=sum(float(r["construction_s"]) for r in c),
                    materialization_s=sum(float(r["materialization_s"]) for r in c),
                    combine_gc_s=sum(float(r["combine_gc_s"]) for r in c),
                    combine_alloc_bytes=sum(int(r["combine_alloc_bytes"]) for r in c)))
                for source, label in ((e, "value"), (d, "direction")):
                    field = "block" if label == "value" else "field"
                    for name in sorted({r[field] for r in source}):
                        group = [r for r in source if r[field] == name]
                        worst = max(group, key=lambda r: float(r["unit_floor_scaled"]))
                        blocks.append(dict(identity, mode=mode, kind=label, field=name,
                            max_scaled=float(worst["unit_floor_scaled"]),
                            absolute_at_max=float(worst["absolute_fro"]),
                            reference_at_max=float(worst["reference_fro"]),
                            max_relative=peak(group, "relative_fro")))
                for depth in sorted({int(r["scan_depth"]) for r in rnk}):
                    rs = [r for r in rnk if int(r["scan_depth"]) == depth]
                    depths.append(dict(identity, mode=mode, scan_depth=depth,
                        rank_max=max(int(r["rank"]) for r in rs),
                        rank_median=statistics.median(int(r["rank"]) for r in rs),
                        nonzero_cross_blocks=sum(float(r["cross_fro"]) > 0 for r in rs),
                        blocks=len(rs), dense_fallbacks=sum(r["dense_fallback"] == "true" for r in rs)))
            validation.append(dict(identity, full_policy_probe_gate=True,
                ordinary_scan_gate=True, all_modes_complete=True,
                max_parent_interface=max(int(r["parent_interface_width"]) for r in sep)))
    if not partial:
        assert not missing, missing
        assert (root / "replay.finished").exists()
        hardware = (root / "hardware.txt").read_text()
        assert "Gold 6438Y+" in hardware
        for name in ("source-check.log", "input-check.log"):
            content = (root / "logs" / name).read_text()
            assert "FAILED" not in content and ": OK" in content
        for path in (root / "logs").glob("slurm-*.err"):
            assert not path.read_text().strip(), str(path)
    write_csv(out / "summary.csv", summaries)
    write_csv(out / "block_errors.csv", blocks)
    write_csv(out / "rank_by_depth.csv", depths)
    write_csv(out / "spectral_ranks.csv", spectral)
    report = {"complete_points": len(validation), "expected_points": 15,
              "complete": not missing and not partial, "strict_campaign_checks": not partial,
              "missing_or_incomplete": missing, "validated": validation,
              "timings_are_benchmarks": False, "compact_kernel": False}
    (out / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    manifest = []
    for path in sorted((root / "results").rglob("*")):
        if path.is_file():
            manifest.append(hashlib.sha256(path.read_bytes()).hexdigest() + "  " + str(path.relative_to(root)).replace("\\", "/"))
    (out / "replay_files.sha256").write_text("\n".join(manifest) + "\n")
    print(json.dumps({"complete_points": len(validation), "missing": len(missing), "complete": report["complete"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--partial", action="store_true")
    args = parser.parse_args()
    analyze(args.root, args.output, args.partial)
