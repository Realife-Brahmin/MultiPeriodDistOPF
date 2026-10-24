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


def analyze(root, out, partial=False, continuation=None):
    out.mkdir(parents=True, exist_ok=True)
    summaries, blocks, depths, spectral, triangles, missing, validation = [], [], [], [], [], [], []
    for case in CASES:
        for tag in TAGS:
            point_root = root
            directory = root / "results" / case / tag
            if continuation and (continuation / "results" / case / tag / "metadata.toml").exists():
                assert not (directory / "metadata.toml").exists(), (case, tag, "duplicate completed point")
                point_root = continuation
                directory = point_root / "results" / case / tag
            log = point_root / "results" / case / (tag + ".log")
            series_log = False
            if not log.exists():
                log = point_root / "logs" / (case + "_series.log")
                series_log = log.exists()
            text = log.read_text(errors="replace") if log.exists() else ""
            if series_log:
                segment = re.search(r"SERIES_POINT_BEGIN [^\n]*tag=" + re.escape(tag) +
                                    r" snapshot=[^\n]*\n(.*?)(?=SERIES_POINT_BEGIN|\Z)", text, re.S)
                text = segment.group(1) if segment else ""
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
            identity = {"case": case, "tag": tag,
                        "source_attempt": "continuation" if point_root == continuation else "primary",
                        "iteration": scalar(meta, "iteration"),
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
                triangle_uncertified = 0
                for r in ii:
                    ref = references[r["depth"], r["lo"], r["hi"]]
                    changed += (r["negative"], r["near_zero"]) != (ref["negative"], ref["near_zero"])
                    # Recorded eigenvalues use (K+K')/2; production uses the
                    # symmetric UPPER triangle. Bound their difference instead
                    # of silently equating these two floating-point matrices.
                    alpha = float(r["symmetry_error"])
                    amax = max(abs(float(r["min_eig"])), abs(float(r["max_eig"])))
                    condition = float(r["abs_condition"])
                    gap = amax / condition if math.isfinite(condition) and condition > 0 else 0.0
                    denom = 1.0 - alpha * alpha / 4.0
                    perturbation = (0.5 * alpha * max(1.0, math.sqrt(n) * amax / math.sqrt(denom))) if denom > 0 else math.inf
                    certified = gap > perturbation + float(r["zero_threshold"])
                    triangle_uncertified += not certified
                    triangles.append(dict(identity, mode=mode, depth=int(r["depth"]),
                        lo=int(r["lo"]), hi=int(r["hi"]), eigen_gap=gap,
                        upper_vs_average_norm_bound=perturbation,
                        inertia_signs_certified=certified))
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
                    upper_triangle_inertia_uncertified=triangle_uncertified,
                    K_condition_max=peak(ii, "abs_condition"),
                    leaf_and_interval_payload_ratio_median=statistics.median(ratios),
                    payload_ratio_max=max(ratios),
                    rank_max=max((int(r["rank"]) for r in rnk), default=0),
                    dense_fallback_fraction=(sum(r["dense_fallback"] == "true" for r in rnk) / len(rnk)) if rnk else 0,
                    raw_dense_fallbacks=sum(r.get("dense_raw_fallback") == "true" for r in c),
                    symmetry_defect_max=peak(c, "symmetry_defect_relative") if "symmetry_defect_relative" in c[0] else 0.0,
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
            validation.append(dict(identity, leaves_sha256=scalar(meta, "leaves_sha256"),
                input_sha256=scalar(meta, "input_sha256"), full_policy_probe_gate=True,
                ordinary_scan_gate=True, all_modes_complete=True,
                max_parent_interface=max(int(r["parent_interface_width"]) for r in sep)))
    if not partial:
        assert not missing, missing
        if continuation:
            original = [r for r in validation if r["source_attempt"] == "primary"]
            remaining = [r for r in validation if r["source_attempt"] == "continuation"]
            assert len(original) == 13 and len(remaining) == 2
            assert {(r["case"], r["tag"]) for r in remaining} == {
                ("large10kC_1ph_T24", "middle"), ("large10kC_1ph_T24", "late")}
            assert (root / "smoke.complete").exists() and not (root / "replay.finished").exists()
            assert (continuation / "remaining.finished").exists()
            # Explicitly retain the KNOWN harness timeout; never silently
            # require clean stderr from a batch which actually failed.
            stderr = list((root / "logs").glob("slurm-*.err"))
            assert len(stderr) == 1
            lines = [s for s in stderr[0].read_text().splitlines() if s.strip()]
            assert len(lines) == 1 and re.search(r"Killed\s+timeout .*12m julia .*replay_series\.jl", lines[0]), lines
            interrupted = (root / "logs" / "large10kC_1ph_T24_series.log").read_text()
            assert re.search(r"signal (?:\(15\)|15): Terminated", interrupted) and "ERROR:" not in interrupted
            assert sum("FROZEN_SERIES_COMPLETE" in (root / "logs" / (c + "_series.log")).read_text()
                       for c in CASES[:2]) == 2
            assert len(list((continuation / "results").rglob("metadata.toml"))) == 2
        else:
            assert (root / "replay.finished").exists()
        for attempt in ([root, continuation] if continuation else [root]):
            hardware = (attempt / "hardware.txt").read_text()
            assert "Gold 6438Y+" in hardware
            checks = ["source-check.log", "input-check.log"]
            if attempt == continuation:
                checks.append("inherited-source-check.log")
            for name in checks:
                content = (attempt / "logs" / name).read_text()
                assert "FAILED" not in content and ": OK" in content
            if not continuation or attempt == continuation:
                for path in (attempt / "logs").glob("slurm-*.err"):
                    assert not path.read_text().strip(), str(path)
    write_csv(out / "summary.csv", summaries)
    write_csv(out / "block_errors.csv", blocks)
    write_csv(out / "rank_by_depth.csv", depths)
    write_csv(out / "spectral_ranks.csv", spectral)
    write_csv(out / "triangle_certificates.csv", triangles)
    write_csv(out / "point_sources.csv", validation)
    report = {"complete_points": len(validation), "expected_points": 15,
              "complete": not missing and not partial, "strict_campaign_checks": not partial,
              "completion_kind": "13_plus_2_continuation" if continuation else "single_attempt",
              "primary_attempt_timeout_preserved": bool(continuation),
              "primary_timeout_classification_validated": bool(continuation) and not partial,
              "missing_or_incomplete": missing, "validated": validation,
              "timings_are_benchmarks": False, "compact_kernel": False}
    (out / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    manifest = []
    for attempt in ([root, continuation] if continuation else [root]):
        prefix = "continuation/" if attempt == continuation else "primary/"
        paths = list((attempt / "results").rglob("*")) + list((attempt / "logs").rglob("*"))
        paths += [attempt / "hardware.txt", attempt / "source.sha256"]
        for path in sorted(paths):
            if path.is_file():
                manifest.append(hashlib.sha256(path.read_bytes()).hexdigest() + "  " + prefix + str(path.relative_to(attempt)).replace("\\", "/"))
    (out / "replay_files.sha256").write_text("\n".join(manifest) + "\n")
    print(json.dumps({"complete_points": len(validation), "missing": len(missing), "complete": report["complete"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--partial", action="store_true")
    parser.add_argument("--continuation", type=Path)
    args = parser.parse_args()
    analyze(args.root, args.output, args.partial, args.continuation)
