"""Iteration counts of the barrier-schedule and hybrid-Hessian runs (run_barrier_sweep.sh).

    python ddp/examples/power_system/barrier_sweep_summary.py [--dir=ddp/results/barrier_sweep/logs]

Per case and variant: iterations and backward passes to the solver's own
tolerance (passes that reused the factors in brackets), and the first
iteration at which the objective is within 0.5% / 1e-4 of the reference (the
MA57 run of the clean batch) with the primal infeasibility below the system
threshold, for the internal infeasibility (int) and for the violation of the
constraints as posed (orig). Untimed runs: counts only.
"""
import glob
import os
import re
import sys

KV = re.compile(r"(\w+)=([-+0-9.eE]+)")
ROW = re.compile(r"^\s+(\d+)\s+([-0-9.]+e[-+]\d+)\s+([0-9.]+e[-+]\d+)\s")
THR = {"ieee123C_1ph": 1e-6, "ieee2522C_1ph": 1e-5, "large10kC_1ph": 1e-5}
NAME = {"ieee123C_1ph": "ieee123", "ieee2522C_1ph": "med2522", "large10kC_1ph": "large10k"}
H = "ddp/results/ipopt_hsl/logs"


def reference(S, T):
    tag = "_screen_lf_pv_hp" if S == "large10kC_1ph" else "_screen_lf_hp"
    path = f"{H}/ipopt_ma57{tag}_{S}_T{T}.log"
    if not os.path.exists(path):
        return None
    for line in open(path, encoding="utf-8", errors="ignore"):
        m = re.search(r"CENTRAL_IPOPT .* objective=([-0-9.eE+]+)", line)
        if m:
            return float(m.group(1))
    return None


def read(path):
    rows, orig = {}, {}
    passes = reused = 0
    done = None
    switch = None
    for line in open(path, encoding="utf-8", errors="ignore"):
        m = ROW.match(line)
        if m:
            rows.setdefault(int(m.group(1)), (float(m.group(2)), float(m.group(3))))
        elif line.startswith("FILTERDDP_ORIGINAL_PRIMAL"):
            d = {k: float(v) for k, v in KV.findall(line)}
            orig.setdefault(int(d["iteration"]), d["original_max"])
        elif line.startswith("FILTERDDP_ITER_TIMING"):
            passes += 1
        elif line.startswith("FILTERDDP_REUSE_PASS"):
            reused += 1
        elif line.startswith("FILTERDDP_HYBRID_SWITCH"):
            switch = int(re.search(r"iteration=(\d+)", line).group(1))
        else:
            m = re.match(r"solve complete: ([0-9.]+) s, iterations=(\d+), status=(\S+)", line)
            if m:
                done = (int(m.group(2)), m.group(3))
    return rows, orig, passes, reused, done, switch


def first(rows, orig, ref, gap, thr, original):
    for k in sorted(rows):
        obj, pr = rows[k]
        if original:
            pr = orig.get(k, float("inf"))
        if abs(obj - ref) / abs(ref) <= gap and pr <= thr:
            return k
    return None


d = "ddp/results/barrier_sweep/logs"
for a in sys.argv[1:]:
    if a.startswith("--dir="):
        d = a.split("=", 1)[1]
cases = {}
for path in sorted(glob.glob(f"{d}/bar_*.log")):
    m = re.match(r"bar_(\w+?_1ph)_T(\d+)_(diag|exact|hybrid)_(.+)\.log", os.path.basename(path))
    if m:
        cases.setdefault((m.group(1), int(m.group(2))), []).append((m.group(4), m.group(3), path))
fmt = lambda k: "--" if k is None else str(k)
for (S, T), runs in sorted(cases.items(), key=lambda c: (list(NAME).index(c[0][0]), c[0][1])):
    ref = reference(S, T)
    print(f"\n{NAME[S]} T={T}  (reference {ref})")
    print(f"  {'variant':12s} {'arm':6s} {'it':>4s} {'passes':>10s} {'st':>2s} | {'0.5% int':>8s} {'orig':>5s} | {'1e-4 int':>8s} {'orig':>5s} | {'final gap':>9s} {'switch':>6s}")
    for name, arm, path in sorted(runs, key=lambda r: (r[1], r[0] != "base", r[0])):
        rows, orig, passes, reused, done, switch = read(path)
        if done is None or ref is None or not rows:
            print(f"  {name:12s} {arm:6s}  (not finished)")
            continue
        last = rows[max(rows)][0]
        cols = [first(rows, orig, ref, g, THR[S], o) for g in (5e-3, 1e-4) for o in (False, True)]
        print(f"  {name:12s} {arm:6s} {done[0]:4d} {passes:5d} ({reused:2d}) {done[1]:>2s} | {fmt(cols[0]):>8s} {fmt(cols[1]):>5s} |"
              f" {fmt(cols[2]):>8s} {fmt(cols[3]):>5s} | {abs(last - ref) / abs(ref):9.1e} {fmt(switch):>6s}")
