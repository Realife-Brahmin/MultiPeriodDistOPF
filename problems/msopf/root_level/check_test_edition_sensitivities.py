"""Validate completed local sensitivity runs without solving or changing their raw evidence.

Usage: python problems/msopf/root_level/check_test_edition_sensitivities.py <run-directory> [...]
Writes validation.json only after every check passes. These are predictions, not measured regions.
"""
import csv
import json
import math
import sys
from pathlib import Path


def rows(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def validate(root):
    metadata = (root / "metadata.txt").read_text()
    assert "exit_code=0" in metadata, "Run is unfinished or exited with an error"
    replay_file, = root.glob("multi_replay_*.csv")
    sensitivity_file, = root.glob("multi_sens_*.csv")
    prediction_file, = root.glob("multi_predicted_*.csv")
    replay, sens, pred = rows(replay_file), rows(sensitivity_file), rows(prediction_file)
    assert len(replay) == 7, f"Expected baseline and six perturbations, found {len(replay)}"
    assert len({r["angles_deg"] for r in replay}) == 7, "Repeated angle configurations"
    for r in replay:
        assert r["status"] in ("LOCALLY_SOLVED", "OPTIMAL")
        assert all(r[k] == "true" for k in ("converged", "pq_band", "passed"))
        for k, limit in (("dP_kW", 1e-3), ("dQ_kvar", 1e-3), ("dV_pu", 1e-7), ("dB_kW", 1e-3)):
            assert math.isfinite(float(r[k])) and 0 <= float(r[k]) <= limit, (k, r[k])
    names = sorted({r["substation"] for r in sens})
    assert len(names) == 4 and len(sens) == 96
    for name in names:
        assert sorted(int(r["hour"]) for r in sens if r["substation"] == name) == list(range(1, 25))
    keys = [k for k in sens[0] if k.startswith("K_")]
    assert len(keys) == 3
    assert all(math.isfinite(float(r[k])) for r in sens for k in ["import_kW", *keys])
    assert len(pred) == 156
    assert {(r["case"], int(r["direction"])) for r in pred} == {(c, i) for c in ("nobatt", "cyclic") for i in range(1, 79)}
    for r in pred:
        u = [float(r[k]) for k in ("u1", "u2", "u3")]
        assert abs(sum(v*v for v in u) - 1) < 1e-12
        reaches = {}
        for name in names:
            hourly = [q for q in sens if q["substation"] == name]
            p = [float(q["import_kW"]) for q in hourly]
            slope = [sum(float(q[k])*v for k, v in zip(keys, u)) for q in hourly]
            if r["case"] == "cyclic":
                reaches[name] = sum(p)/-sum(slope) if sum(slope) < 0 else math.inf
            else:
                reaches[name] = min((a/-b if b < 0 else math.inf) for a, b in zip(p, slope))
        expected = min(reaches.values())
        observed = float(r["predicted_deg"])
        assert observed > 0, "Equal-angle operating point is outside a predicted import half-space"
        assert (math.isinf(expected) and math.isinf(observed)) or abs(expected-observed) <= 2e-5*max(1, abs(expected)), (expected, observed)
        assert r["limiting_substation"] in names or math.isinf(observed)
    summary = {
        "validated": True,
        "kind": "sensitivity-based estimates, not measured feasibility boundaries",
        "solves": len(replay), "hours_per_solve": 24, "sensitivity_rows": len(sens), "predicted_rows": len(pred),
        "worst_replay": {k: max(float(r[k]) for r in replay) for k in ("dP_kW", "dQ_kvar", "dV_pu", "dB_kW")},
        "mean_import_kW": {n: sum(float(r["import_kW"]) for r in sens if r["substation"] == n)/24 for n in names},
        "mean_sensitivity_kW_per_deg": {n: {k: sum(float(r[k]) for r in sens if r["substation"] == n)/24 for k in keys} for n in names},
    }
    (root / "validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"PASS {root.name}: 7 replays, 96 hourly sensitivity rows, 156 predicted reaches")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    for directory in sys.argv[1:]:
        validate(Path(directory))
