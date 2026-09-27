"""Markdown tables from ddp/results/kkt_ordering/kkt_ordering_benchmark.csv.

Per (system, arm): every configuration's median ordering, factorization and
(n_x+1)-column solve time, factor-plus-wide total relative to FilterDDP's
UMFPACK baseline, fill as each solver stores it, ordering-only fill, and the
wide-solve residual. Also a one-line "best by factor+wide" and "least
ordering-only fill" per case.

    python ddp/examples/power_system/summarize_kkt_ordering.py [csv] > summary.md
"""
import csv
import sys
from collections import defaultdict

path = sys.argv[1] if len(sys.argv) > 1 else "ddp/results/kkt_ordering/kkt_ordering_benchmark.csv"
rows = list(csv.DictReader(open(path)))
cases = defaultdict(list)
for r in rows:
    cases[(r["system"], r["arm"])].append(r)

f = lambda r, k: float(r[k]) if r.get(k) not in (None, "") else float("nan")
order = {"ieee123C_1ph": 0, "ieee2522C_1ph": 1, "large10kC_1ph": 2}
for (system, arm) in sorted(cases, key=lambda c: (order.get(c[0], 9), c[1])):
    rs = cases[(system, arm)]
    ok = [r for r in rs if r["status"] == "ok"]
    base = next(r for r in ok if r["config"] == "UMFPACK_default")
    b = f(base, "factor_s") + f(base, "solvewide_s")
    r0 = rs[0]
    print(f"\n### {system}, {arm} Hessian: n = {int(r0['n']):,}, nnz(K) = {int(r0['nnz_K']):,}, "
          f"RHS = {r0['rhs_cols']} columns\n")
    print("| config | ordering used | ordering s | factor s | 1-RHS solve s | wide solve s "
          "| factor+wide s | vs baseline | fill (stored) | fill (ordering only) | wide residual | pivots |")
    print("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for r in rs:
        if r["status"] != "ok":
            print(f"| {r['config']} | failed | | | | | | | | | | {r['note']} |")
            continue
        fw = f(r, "factor_s") + f(r, "solvewide_s")
        print(f"| {r['config']} | {r['ordering_used']} ({r['strategy_used']}) | {f(r,'ordering_s'):.4g} "
              f"| {f(r,'factor_s'):.4g} | {f(r,'solve1_s'):.3g} | {f(r,'solvewide_s'):.4g} | {fw:.4g} "
              f"| {fw / b:.2f}x | {f(r,'fill_ratio'):.2f} | {f(r,'struct_fill_sym_order'):.2f} "
              f"| {f(r,'relres_wide'):.1e} | {r['pivot_stats']} |")
    best = min(ok, key=lambda r: f(r, "factor_s") + f(r, "solvewide_s"))
    leanest = min(ok, key=lambda r: f(r, "struct_fill_sym_order"))
    wide_share = f(base, "solvewide_s") / (f(base, "ordering_s") + b)
    print(f"\nBaseline wide-solve share of ordering+factor+wide: {100 * wide_share:.0f}%. "
          f"Fastest factor+wide: {best['config']} "
          f"({(f(best,'factor_s') + f(best,'solvewide_s')) / b:.2f}x baseline). "
          f"Least ordering-only fill: {leanest['config']} ({f(leanest,'struct_fill_sym_order'):.2f}).")
