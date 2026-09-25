"""One row per full FilterDDP run in ddp/results/kkt_ordering/fullrun_blocked/.

Reads each fddp_*.log: iterations, time to near-optimality, summed per-stage
factorization and (n_x+1)-column solve time (iteration 0 excluded, it carries
compilation), objective, final equality residual, the blocked-solve setting and
the median background load logged during the run.

    python ddp/examples/power_system/summarize_blocked_fullrun.py
"""
import csv
import glob
import re
import statistics

OUT = "ddp/results/kkt_ordering/fullrun_blocked"
rows = []
for log in sorted(glob.glob(f"{OUT}/fddp_*.log")):
    txt = open(log).read()
    factor = solve = 0.0
    for m in re.finditer(r"FILTERDDP_TIMING iteration=(\d+) .*?factor_s=([\d.]+) solve_s=([\d.]+)", txt):
        if int(m.group(1)) == 0:
            continue
        factor += float(m.group(2)); solve += float(m.group(3))
    wall = re.search(r"solve complete: ([\d.]+) s, iterations=(\d+)", txt)
    near = re.search(r"FILTERDDP_NEAR_OPT iteration=(\d+) elapsed_s=([\d.]+)", txt)
    if not (wall and near):
        print(f"skipping unfinished run {log}")
        continue
    env = re.search(r"PIPELINE_ENV (.*)", txt).group(1)
    blocked = re.search(r"blocked_solve=(\S+)", env).group(1)
    rewrites = " ".join(re.findall(r"(direct_diag=\S+|triplet=\S+|kkt_pattern_cache=\S+)", env)) or "not recorded (off)"
    try:
        load = [float(r[1]) for r in csv.reader(open(log.replace("fddp_", "load_").replace(".log", ".csv")))
                if r and r[0][:1].isdigit()]
        load = f"{statistics.median(load):.2f}" if load else ""
    except FileNotFoundError:
        load = ""
    rows.append([log.split("fddp_")[1][:-4], blocked, rewrites, wall.group(2), wall.group(1),
                 near.group(1), near.group(2), f"{factor:.3f}", f"{solve:.3f}",
                 re.search(r"FilterDDP objective=([-\d.eE+]+)", txt).group(1),
                 re.search(r"max_equality_residual=([\d.eE+-]+)", txt).group(1), load])

with open(f"{OUT}/fullrun_blocked_summary.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["run", "blocked_solve", "assembly_rewrites", "iterations", "wall_s",
                "near_opt_iteration", "near_opt_elapsed_s", "factor_s_sum_excl_iter0",
                "solve_s_sum_excl_iter0", "objective", "max_equality_residual", "median_background_cores"])
    w.writerows(rows)
for r in rows:
    print(f"{r[0]:48s} it={r[3]:>3s} near_opt_s={float(r[6]):8.1f} solve={float(r[8]):8.1f} obj={r[9]} load={r[11]}")
