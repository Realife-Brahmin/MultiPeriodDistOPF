"""Summary of the PINT_SCAN lines of run_pint_scan.sh logs (pint_scan.jl).

    python ddp/examples/power_system/pint_scan_summary.py ddp/results/pint_scan/logs/pint_*_scan.log

One row per log and BLAS thread count: medians over the scanned backward
passes. Times in ms per backward pass.
  sweep     the value recursion in sequence, from the stage data
  count     the same by pairwise combination, one worker per period: slowest
            combination of each round, plus the slowest stage function and the
            slowest final battery solve (max: slowest; med: median of each)
  solver    the solver's own sequential part for the same passes
"""
import re
import statistics
import sys

num = re.compile(r"(\w+)=([-+0-9.eE]+)")


def med(v):
    return statistics.median(v) if v else float("nan")


print(f"{'log':34s} {'thr':>3s} {'T':>3s} {'nB':>5s} {'scans':>5s} {'depth':>5s} {'err':>8s} {'cholF':>5s} "
      f"{'sweep':>9s} {'cnt_max':>9s} {'cnt_med':>9s} {'x_max':>6s} {'x_med':>6s} {'merge':>8s} {'step':>8s} {'solver':>9s}")
for path in sys.argv[1:]:
    rows = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        if not line.startswith("PINT_SCAN "):
            continue
        d = {k: float(v) for k, v in num.findall(line)}
        rows.setdefault(int(d["blas_threads"]), []).append(d)
    name = path.replace("\\", "/").split("/")[-1].replace("pint_", "").replace(".log", "")
    for thr in sorted(rows):
        r = rows[thr]
        err = max(max(d["scan_err_Vxx"], d["scan_err_Vx"], d["sweep_err_Vxx"], d["sweep_err_Vx"]) for d in r)
        chol = int(sum(d["chol_fail_sweep"] + d["chol_fail_scan"] for d in r))
        print(f"{name:34s} {thr:3d} {int(r[0]['T']):3d} {int(r[0]['nB']):5d} {len(r):5d} {int(r[0]['depth']):5d} {err:8.1e} {chol:5d} "
              f"{1e3 * med([d['sweep_sum_s'] for d in r]):9.3f} {1e3 * med([d['count_max_s'] for d in r]):9.3f} "
              f"{1e3 * med([d['count_med_s'] for d in r]):9.3f} {med([d['speedup_max'] for d in r]):6.2f} "
              f"{med([d['speedup_med'] for d in r]):6.2f} {1e3 * med([d['merge_med_s'] for d in r]):8.3f} "
              f"{1e3 * med([d['sweep_step_med_s'] for d in r]):8.3f} {1e3 * med([d['solver_seq_sum_s'] for d in r]):9.3f}")
