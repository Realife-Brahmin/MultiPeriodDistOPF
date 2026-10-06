"""Do two FilterDDP logs follow the same iterates?

For changes that keep the mathematics but not the order of floating-point
operations (a different elimination order, say), identical output is too much
to ask. This compares, pass by pass: the outcome, step size, backtracks and
barrier parameter of every iteration (FILTERDDP_ITER_TIMING), and the
equality-residual statistics (FILTERDDP_FEASIBILITY), and reports the number
of passes and the largest relative difference. Final objectives are compared
too.

    python compare_runs.py <reference_log> <log>
"""
import re
import sys

KV = re.compile(r"(\w+)=(\S+)")


def lines(path, tag):
    out = []
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith(tag):
            out.append(dict(KV.findall(line)))
    return out


def reldiff(a, b):
    a, b = float(a), float(b)
    return abs(a - b) / max(abs(a), abs(b), 1e-300) if (a != 0 or b != 0) else 0.0


def objective(path):
    m = None
    for line in open(path, encoding="utf-8", errors="replace"):
        g = re.search(r"FilterDDP objective=([-0-9.eE+]+)", line)
        m = g.group(1) if g else m
    return m


ref, new = sys.argv[1], sys.argv[2]
ia, ib = lines(ref, "FILTERDDP_ITER_TIMING"), lines(new, "FILTERDDP_ITER_TIMING")
fa, fb = lines(ref, "FILTERDDP_FEASIBILITY"), lines(new, "FILTERDDP_FEASIBILITY")
same_trace = len(ia) == len(ib) and all(
    x["iteration"] == y["iteration"] and x["outcome"] == y["outcome"] and x["backtracks"] == y["backtracks"]
    for x, y in zip(ia, ib))
step = max((reldiff(x["step_size"], y["step_size"]) for x, y in zip(ia, ib)), default=0.0)
mu = max((reldiff(x["mu"], y["mu"]) for x, y in zip(ia, ib)), default=0.0)
res = max((reldiff(x[k], y[k]) for x, y in zip(fa, fb) for k in ("equality_rms", "equality_max")), default=0.0)
oa, ob = objective(ref), objective(new)
print(f"passes {len(ia)} vs {len(ib)}; same outcomes and backtracks: {same_trace}; "
      f"largest relative difference: step {step:.1e}, mu {mu:.1e}, equality residuals {res:.1e}; "
      f"objective {oa} vs {ob} ({reldiff(oa, ob):.1e})")
