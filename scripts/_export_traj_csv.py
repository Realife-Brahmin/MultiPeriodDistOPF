#!/usr/bin/env python3
"""Export BF + tADMM convergence trajectories to tidy CSVs for PGFPlots.
BF: parsed from ipopt_bf.log, split into infeasible/feasible, truncated at BF's
first feasible-and-within-0.5% iterate. tADMM: trajectory from convergence_data.csv,
truncated at the authoritative near-optimal milestone in near_optimal_summary.csv
(which uses the correct benchmark: BF when converged, final tADMM when BF fails).
Lines emitted in HOURS. Prints scalars to bake into the PGF TeX.
"""
import re, os, statistics

FEAS_TOL = 1e-4
GAP_TOL  = 0.005
PROC = "envs/tadmm/processedData"
OUT  = "../IAS-Trans-2025-Scaling-MPOPF-Computation-via-Temporal-Decomposition/figures"

CASES = [
    dict(key="med2522_T144", bf="ieee2552C_1ph_T144",
         tcsv="ieee2552C_1ph_T144/convergence_data.csv"),
    dict(key="large10k_T48", bf="large10kC_1ph_T48",
         tcsv="large10kC_1ph_T48/rho_sweep/rho_30000/convergence_data.csv"),
]

def parse_ipopt(path):
    rx = re.compile(r"^\s*(\d+)r?\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)")
    out = []
    for line in open(path):
        m = rx.match(line)
        if m: out.append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
    return out

def parse_total(path):
    for line in open(path):
        m = re.search(r"Wall-clock time:\s*([\d.]+)\s*seconds", line)
        if m: return float(m.group(1))
    return float("nan")

def parse_csv(path):
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    hdr = lines[0].split(",")
    rows = []
    for l in lines[1:]:
        v = l.split(",")
        d = {}
        for i, c in enumerate(hdr):
            try: d[c] = float(v[i])
            except: d[c] = v[i] if i < len(v) else None
        rows.append(d)
    return hdr, rows

def downsample(pts, target=160):
    n = len(pts)
    if n <= target: return pts
    stride = max(1, round(n/target))
    idx = list(range(0, n, stride))
    if n-1 not in idx: idx.append(n-1)
    return [pts[i] for i in sorted(set(idx))]

def write_csv(name, pts):
    with open(os.path.join(OUT, name), "w", newline="") as f:
        f.write("t,obj\n")
        for t, o in pts: f.write(f"{t:.6f},{o:.4f}\n")
    return len(pts)

for c in CASES:
    bf = parse_ipopt(os.path.join(PROC, c["bf"], "ipopt_bf.log"))
    total = parse_total(os.path.join(PROC, c["bf"], "results_socp_bf.txt"))
    n = len(bf); dt = total/n
    t   = [i*dt for i in range(n)]
    obj = [b[1] for b in bf]
    feas = [b[2] <= FEAS_TOL for b in bf]
    bf_ref = obj[-1]
    first_feas = next((i for i in range(n) if feas[i]), None)
    first_near = next((i for i in range(n) if feas[i] and abs(obj[i]-bf_ref) <= GAP_TOL*bf_ref), None)
    bf_end = first_near if first_near is not None else n-1

    # tADMM trajectory + authoritative NO summary
    _, trows = parse_csv(os.path.join(PROC, c["tcsv"]))
    sdir = os.path.dirname(os.path.join(PROC, c["tcsv"]))
    _, srows = parse_csv(os.path.join(sdir, "near_optimal_summary.csv"))
    s = srows[0]
    no_iter = int(s["near_opt_iter"]); no_t = s["near_opt_eff_time"]
    no_obj = s["near_opt_obj"]; ref = s["ref_obj"]; ref_src = s["ref_source"]

    S = 1/3600.0
    # tADMM curve up to the NO iteration, dropping wild leading consensus-forming
    # iterates (objective far below the benchmark) that would draw as a spike.
    ref_half = 0.5*ref
    start = next((i for i in range(len(trows)) if trows[i]["objective"] >= ref_half), 0)
    tser = [(trows[i]["cum_eff_time"]*S, trows[i]["objective"])
            for i in range(start, len(trows)) if int(trows[i]["iteration"]) <= no_iter]
    # BF split (overlap at first_feas for a continuous line)
    if first_feas is None or first_feas > bf_end:
        inf_idx = list(range(0, bf_end+1)); fea_idx = []
    else:
        inf_idx = list(range(0, first_feas+1)); fea_idx = list(range(first_feas, bf_end+1))
    bf_inf = downsample([(t[i]*S, obj[i]) for i in inf_idx])
    bf_fea = downsample([(t[i]*S, obj[i]) for i in fea_idx]) if fea_idx else []

    n1 = write_csv(f"traj_{c['key']}_bf_inf.csv", bf_inf)
    n2 = write_csv(f"traj_{c['key']}_bf_fea.csv", bf_fea) if bf_fea else 0
    n3 = write_csv(f"traj_{c['key']}_tadmm.csv", downsample(tser))

    bf_no_t = t[first_near]*S; bf_no_o = obj[first_near]
    t_no_t = no_t*S
    spd = bf_no_t / t_no_t
    obj_all  = [obj[i] for i in range(bf_end+1)] + [p[1] for p in tser]
    obj_conv = [obj[i] for i in fea_idx] + [p[1] for p in tser]
    ycap = statistics.quantiles(obj_all, n=100)[94]
    yhi = max(ycap, max(obj_conv))*1.14
    ylo = min(obj_all)*0.96
    xmax = max(bf_no_t, t_no_t)*1.05
    print(f"\n=== {c['key']} === rows bf_inf={n1} bf_fea={n2} tadmm={n3}")
    print(f"  ref_obj={ref:.2f} ({ref_src})")
    print(f"  BF NO:    t={bf_no_t:.4f}h obj={bf_no_o:.2f} (k={first_near}, first_feas={first_feas}, ntot={n})")
    print(f"  tADMM NO: t={t_no_t:.4f}h obj={no_obj:.2f} (iter={no_iter})")
    print(f"  speedup={spd:.2f}x  ylo={ylo:.1f} yhi={yhi:.1f} xmax={xmax:.4f}")
