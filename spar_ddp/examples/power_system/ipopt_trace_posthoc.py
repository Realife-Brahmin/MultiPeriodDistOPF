"""When an Ipopt run first met a looser criterion, from its IPOPT_TRACE=1 log.

The counterpart of near_opt_posthoc.py, with the same rule: the first
iteration whose objective is within GAP of the final objective and whose
primal infeasibility is below the system threshold, and the time at which
that iterate was known (the callback of that iteration).

Two primal infeasibilities are reported, because they differ a great deal:
  internal  the callback's inf_pr: an inequality counts through its slack,
            d(x) - s. FilterDDP's primal_inf is of this kind (its cone rows
            carry an explicit slack), so this is the like-for-like one.
  original  the inf_pr column of Ipopt's own iteration table (read from
            <log>_ipoptlog.txt): the violation of the constraints as posed,
            zero whenever an inequality holds.

    python ipopt_trace_posthoc.py [--gaps=5e-3,1e-4,1e-5] <primal_threshold> <ipopt_log> [...]
"""
import os
import re
import sys

KV = re.compile(r"(\w+)=([-+0-9.eE]+)")
ROW = re.compile(r"^\s*(\d+)r?\s+([-0-9.]+e[-+]\d+)\s+([0-9.]+e[-+]\d+)\s")


def analyse(path, primal, gaps):
    trace = []
    final = None
    for line in open(path, encoding="utf-8", errors="ignore"):
        if line.startswith("IPOPT_ITER "):
            trace.append({k: float(v) for k, v in KV.findall(line)})
        else:
            m = re.search(r"CENTRAL_IPOPT .* status=(\S+) iterations=(\d+) objective=([-0-9.eE+]+) solve_time_s=([0-9.]+)", line)
            if m:
                final = (m.group(1), int(m.group(2)), float(m.group(3)), float(m.group(4)))
    if not trace or final is None:
        print(f"{path}: no IPOPT_ITER lines or no CENTRAL_IPOPT line")
        return
    status, iters, obj, solve_s = final
    orig = {}
    table = path[:-4] + "_ipoptlog.txt"
    if os.path.exists(table):
        for line in open(table, encoding="utf-8", errors="ignore"):
            m = ROW.match(line)
            if m:
                orig[int(m.group(1))] = float(m.group(3))
    name = path.replace(chr(92), "/").split("/")[-1]
    print(f"{name}: status={status} iterations={iters} solve_s={solve_s:.3f} first_iterate_s={trace[0]['elapsed_s']:.3f}"
          f" objective_ratio={obj / trace[-1]['objective']:.9f}")
    for gap in gaps:
        parts = []
        for label, pr in (("internal", lambda d: d["inf_pr"]), ("original", lambda d: orig.get(int(d["iteration"]), float("inf")))):
            hit = next((d for d in trace if abs(d["objective"] - obj) / abs(obj) <= gap and pr(d) <= primal), None)
            parts.append(f"{label}: " + ("not reached" if hit is None else
                         f"iteration {int(hit['iteration'])}, {hit['elapsed_s']:.3f} s ({hit['elapsed_s'] / solve_s:.0%})"))
        print(f"  gap<={gap:g}, primal<={primal:g}  " + "   ".join(parts))


args = sys.argv[1:]
gaps = [5e-3]
if args and args[0].startswith("--gaps="):
    gaps = [float(g) for g in args.pop(0).split("=", 1)[1].split(",")]
primal = float(args.pop(0))
for p in args:
    analyse(p, primal, gaps)
