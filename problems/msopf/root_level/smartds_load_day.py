"""Hourly loads of sf6937 from SMART-DS's own customer profiles (user, 2026-11-03).

Each deck load (the customers of one service transformer, LoadMembers.csv of convert_to_1ph.jl) gets, for each
hour h = 1..24 (hour ending h:00), the sum over its customers of nominal kW x the mean of the four 15-minute
multipliers of that hour, and the same for kvar with the kvar multipliers, divided by 3 (the per-phase
equivalent). The day is the dataset's peak day, 148 (28 May 2018), fetched by tools/fetch_smartds_day.py.
Writes test_networks/smartds_p1u_xf_1ph/LoadDay148.csv and checks the total at the dataset's own peak instant
(14:00) against its analysis/Summary_data.csv.

    python -I problems/msopf/root_level/smartds_load_day.py
"""
import csv
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
SRC = os.path.join(ROOT, "test_networks", "sources", "smartds_sfo_p1u")
DECK = os.path.join(ROOT, "test_networks", "smartds_p1u_xf_1ph")
DAY = 148


def main():
    files = {}                                         # shape -> (kW multiplier file, kvar multiplier file)
    for m in re.finditer(r"New Loadshape\.(\S+).*?mult = \(file=(?:\.\./)+profiles/([^)\s]+)\)\s+qmult = \(file=(?:\.\./)+profiles/([^)\s]+)\)",
                         open(os.path.join(SRC, "LoadShapes.dss")).read()):
        files[m.group(1).lower()] = (m.group(2), m.group(3))
    prof = {}
    with open(os.path.join(SRC, f"profiles_day{DAY}.csv")) as f:
        for r in list(csv.reader(f))[1:]:
            prof[r[0]] = [float(x) for x in r[1:]]
    loads = {}                                         # deck load -> 24 kW, 24 kvar (per-phase equivalent)
    at14 = 0.0                                         # 3-phase customer kW at 14:00, as the dataset's analysis
    with open(os.path.join(DECK, "LoadMembers.csv")) as f:
        for r in csv.DictReader(f):
            kp, kq = files[r["shape"]]
            p, q = prof[kp], prof[kq]
            P, Q = float(r["kW"]), float(r["kvar"])
            v = loads.setdefault(r["load"], [0.0] * 48)
            for h in range(24):
                v[h] += P * sum(p[4 * h:4 * h + 4]) / 4 / 3
                v[24 + h] += Q * sum(q[4 * h:4 * h + 4]) / 4 / 3
            at14 += P * p[56]
    out = os.path.join(DECK, f"LoadDay{DAY}.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["load"] + [f"kW_{h}" for h in range(1, 25)] + [f"kvar_{h}" for h in range(1, 25)])
        for n in sorted(loads):
            w.writerow([n] + [f"{x:.6g}" for x in loads[n]])
    tot = [3 * sum(v[h] for v in loads.values()) for h in range(48)]
    with open(os.path.join(SRC, "analysis", "Summary_data.csv")) as f:
        ref = float(list(csv.DictReader(f))[0]["Total peak time real customer power (kW)"])
    print(f"wrote {out}: {len(loads)} loads")
    print(f"customer kW at 14:00 on day {DAY}: {at14:.1f}, dataset's analysis {ref:.1f} (ratio {at14 / ref:.5f})")
    print("hourly totals, 3-phase MW (hour ending): " + ", ".join(f"{h + 1}:{tot[h] / 1000:.1f}" for h in range(24)))
    print(f"peak hour {max(range(24), key=lambda h: tot[h]) + 1}, min/max {min(tot[:24]) / max(tot[:24]):.3f}")


if __name__ == "__main__":
    main()
