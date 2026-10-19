"""Which bound stopped the step, per iteration, from a FilterDDP log.

Reads the FILTERDDP_FTB lines (FILTERDDP_FTB_DIAGNOSTIC=1): one per trial step
rejected by the fraction-to-boundary rule, with the stage, the test that
failed and the worst entry. Maps the control index to its variable type using
the driver's control layout
    ps, qs, P[L], Q[L], v[N], ell[L], pb[B], qnorm[D], soc_slack[L], energy_slack[B]
and prints, per iteration, the accepted step and what rejected each larger
trial; then a summary of the last rejection before the accepted step (the
bound that set the step size).

    python step_blockers.py <N buses> <L lines> <B batteries> <D ders> <log> [--all]
"""
import re
import sys
from collections import Counter, defaultdict

N, L, B, D = (int(a) for a in sys.argv[1:5])
log = sys.argv[5]
show_all = "--all" in sys.argv

edges = [("P_Subs", 1), ("Q_Subs", 1), ("P", L), ("Q", L), ("v", N), ("ell", L),
         ("P_B", B), ("qnorm", D), ("soc_slack", L), ("energy_slack", B)]


def var(i):
    for name, n in edges:
        if i <= n:
            return name, i
        i -= n
    return "?", i


rej = defaultdict(list)       # iteration -> [(step, stage, kind, var, pos, failed, before, after)]
acc = {}                      # iteration -> (step, filter backtracks, outcome)
for line in open(log, encoding="utf-8", errors="ignore"):
    if line.startswith("FILTERDDP_FTB "):
        d = dict(kv.split("=") for kv in line.split()[1:])
        name, pos = var(int(d["worst_index"]))
        rej[int(d["iteration"])].append((float(d["step"]), int(d["stage"]), d["kind"], name, pos,
                                         int(d["failed"]), float(d["before"]), float(d["after"])))
    elif line.startswith("FILTERDDP_ITER_TIMING"):
        d = dict(kv.split("=") for kv in line.split()[1:] if "=" in kv)
        if "step_size" in d and d.get("outcome") in ("accepted", "forward_failed"):
            acc[int(d["iteration"])] = (float(d["step_size"]), int(d["backtracks"]), d["outcome"])

setter = Counter()
steps = Counter()
for k in sorted(acc):
    step, bt, outcome = acc[k]
    steps[step] += 1
    r = rej.get(k, [])
    last = r[-1] if r else None
    label = "none (full step)" if last is None else f"{last[2]}:{last[3]}"
    if last is not None and last[3] == "v" and last[4] == 1:
        label += "[substation]"
    setter[label] += 1
    if show_all:
        what = "; ".join(f"{s:.3g}->{kind}:{name}#{pos}@t{t} ({n} entries, {b:.1e}->{a:.1e})"
                          for s, t, kind, name, pos, n, b, a in r)
        print(f"it {k:3d} step {step:.3g} filter_backtracks {bt}  {what}")
total = sum(setter.values())
print(f"STEP_BLOCKERS {log.replace(chr(92), '/').split('/')[-1]}: {total} iterations")
print("  accepted step sizes: " + ", ".join(f"{s:.3g} x{n}" for s, n in sorted(steps.items(), reverse=True)))
for label, n in setter.most_common():
    print(f"  last bound to reject a larger step: {label:42s} {n:4d} iterations ({100 * n / total:.0f}%)")
