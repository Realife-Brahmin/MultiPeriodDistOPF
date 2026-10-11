"""Tables of the clean batch of 2026-10-11 (run_clean_hp_batch.sh).

For every case, both solvers under one stopping rule: the first iterate whose
objective is within GAP of the reference optimum (the MA57 run's final
objective) and whose primal infeasibility is below the system threshold.

    python ddp/examples/power_system/clean_hp_table.py [--gaps=5e-3,1e-4,1e-5] [--primal=internal|original] [--threshold=1e-5]

Ipopt: time of that iterate from its trace (IPOPT_TRACE), per linear solver,
and its full solve time (tol 1e-6).
FilterDDP (tolerance 1e-5, diagonal and exact Hessian): time of that iterate
  seq     as measured (sequential run, compilation excluded by the warm-up)
  count   one worker per period (parsim_from_log.py's accounting, up to that
          iterate): slowest period for the independent parts, the sequential
          sweep as measured
  paired  exact Hessian only, a projection: the sweep replaced by the pairwise
          combination (rounds from the pint_scan2 logs of the same system and
          T, slowest combination of each round) and every stage's own
          sequential work counted as the slowest stage
Primal infeasibility: `internal` is each solver's own (an inequality counted
through its slack); `original` is the violation of the constraints as posed.
"""
import glob
import os
import re
import statistics
import sys

ROOT = "ddp/results"
KV = re.compile(r"(\w+)=([-+0-9.eE]+)")
THR = {"ieee123C_1ph": 1e-6, "ieee2522C_1ph": 1e-5, "large10kC_1ph": 1e-5}
CASES = [("ieee123C_1ph", 6), ("ieee123C_1ph", 24), ("ieee123C_1ph", 96), ("ieee123C_1ph", 144),
         ("ieee2522C_1ph", 6), ("ieee2522C_1ph", 24), ("ieee2522C_1ph", 96), ("ieee2522C_1ph", 144),
         ("large10kC_1ph", 6), ("large10kC_1ph", 24), ("large10kC_1ph", 48), ("large10kC_1ph", 96),
         ("large10kC_1ph", 144)]
NAME = {"ieee123C_1ph": "ieee123", "ieee2522C_1ph": "med2522", "large10kC_1ph": "large10k"}


def fields(line):
    return {k: float(v) for k, v in KV.findall(line)}


def ipopt(path, gap, primal, original, ref=None):
    """(time at the rule, iteration, full solve time, iterations, objective) or None"""
    if not os.path.exists(path):
        return None
    trace, final = [], None
    for line in open(path, encoding="utf-8", errors="ignore"):
        if line.startswith("IPOPT_ITER "):
            trace.append(fields(line))
        else:
            m = re.search(r"CENTRAL_IPOPT .* status=(\S+) iterations=(\d+) objective=([-0-9.eE+]+) solve_time_s=([0-9.]+)", line)
            if m:
                final = (m.group(1), int(m.group(2)), float(m.group(3)), float(m.group(4)))
    if final is None or not trace:
        return None
    orig = {}
    table = path[:-4] + "_ipoptlog.txt"
    if os.path.exists(table):
        row = re.compile(r"^\s*(\d+)r?\s+([-0-9.]+e[-+]\d+)\s+([0-9.]+e[-+]\d+)\s")
        for line in open(table, encoding="utf-8", errors="ignore"):
            m = row.match(line)
            if m:
                orig[int(m.group(1))] = float(m.group(3))
    obj = final[2] if ref is None else ref
    for d in trace:
        pr = orig.get(int(d["iteration"]), float("inf")) if original else d["inf_pr"]
        if abs(d["objective"] - obj) / abs(obj) <= gap and pr <= primal:
            return d["elapsed_s"], int(d["iteration"]), final[3], final[1], final[2]
    return None, None, final[3], final[1], final[2]


def ddp(path):
    """per-pass records of a FILTERDDP_PARSIM=1 log"""
    if not os.path.exists(path):
        return None
    rows, orig, passes = {}, {}, []
    back, fwd, stage_seq = [], [], []
    solve_s = None
    status = None
    row = re.compile(r"^\s+(\d+)\s+([-0-9.]+e[-+]\d+)\s+([0-9.]+e[-+]\d+)\s")
    for line in open(path, encoding="utf-8", errors="ignore"):
        m = row.match(line)
        if m:
            rows.setdefault(int(m.group(1)), (float(m.group(2)), float(m.group(3))))
        elif line.startswith("FILTERDDP_ORIGINAL_PRIMAL"):
            d = fields(line)
            orig.setdefault(int(d["iteration"]), d["original_max"])
        elif line.startswith("FILTERDDP_PARSIM_STAGE"):
            stage_seq.append(fields(line)["seq_s"])
        elif line.startswith("FILTERDDP_PARSIM_BACKWARD"):
            d = fields(line)
            d["seq_max_s"] = max(stage_seq) if stage_seq else 0.0
            back.append(d)
            stage_seq = []
        elif line.startswith("FILTERDDP_PARSIM_FORWARD"):
            fwd.append(fields(line))
        elif line.startswith("FILTERDDP_ITER_TIMING"):
            f = fields(line)
            # pair_*: the pairwise combinations of a paired pass (paired_pass.jl); absent otherwise
            b = {k: sum(x.get(k, 0.0) for x in back) for k in ("pre_sum_s", "pre_max_s", "seq_sum_s", "post_sum_s", "post_max_s", "seq_max_s", "pair_sum_s", "pair_max_s")} if back else None
            w = {k: sum(x[k] for x in fwd) for k in ("dir_sum_s", "dir_max_s", "trial_sum_s", "trial_max_s")} if fwd else {k: 0.0 for k in ("dir_sum_s", "dir_max_s", "trial_sum_s", "trial_max_s")}
            passes.append((int(f["iteration"]), f, b, w, len(back)))
            back, fwd = [], []
        else:
            m = re.match(r"solve complete: ([0-9.]+) s, iterations=(\d+), status=(\S+)", line)
            if m:
                solve_s, status = float(m.group(1)), m.group(3)
    if solve_s is None or not passes or passes[0][2] is None:
        return None
    return dict(rows=rows, orig=orig, passes=passes, solve_s=solve_s, status=status,
                setup=solve_s - sum(p[1]["total_s"] for p in passes), iters=passes[-1][0])


def ddp_at(run, ref, gap, primal, original, pair=None):
    """(seq, count, paired, iteration) at the rule; pair = (leaf_max, rounds_max) per pass or None"""
    el_seq = el_cnt = el_pair = run["setup"]
    seen = set()
    for k, f, b, w, nb in run["passes"]:
        save_f = (w["dir_sum_s"] - w["dir_max_s"]) + (w["trial_sum_s"] - w["trial_max_s"])
        back_cnt = (f["backward_s"] - (b["pre_sum_s"] - b["pre_max_s"]) - (b["post_sum_s"] - b["post_max_s"])
                    - (b["pair_sum_s"] - b["pair_max_s"]))
        # paired: the sweep's sum over stages becomes the rounds plus the slowest stage
        back_pair = back_cnt - b["seq_sum_s"] + b["seq_max_s"] + nb * (pair[0] + pair[1]) if pair else float("nan")
        if k not in seen and k in run["rows"]:
            seen.add(k)
            obj, pr = run["rows"][k]
            if original:
                pr = run["orig"].get(k, float("inf"))
            if abs(obj - ref) / abs(ref) <= gap and pr <= primal:
                return el_seq + f["backward_s"], el_cnt + back_cnt, el_pair + back_pair, k
        el_seq += f["total_s"]
        el_cnt += back_cnt + f["forward_s"] - save_f + (f["total_s"] - f["backward_s"] - f["forward_s"])
        el_pair += back_pair + f["forward_s"] - save_f + (f["total_s"] - f["backward_s"] - f["forward_s"])
    return None, None, None, None


def pairing(S, T):
    """(leaf_max, rounds_max) seconds per pass, one BLAS thread, medians over the scanned passes"""
    path = f"{ROOT}/pint_scan/logs/pint_{S}_T{T}_scan2.log"
    if not os.path.exists(path):
        return None
    r = [fields(l) for l in open(path, encoding="utf-8", errors="ignore") if l.startswith("PINT_SCAN ") and "blas_threads=1 " in l]
    if not r:
        return None
    return statistics.median(d["leaf_max_s"] for d in r), statistics.median(d["rounds_max_s"] for d in r)


def fmt(x, w=7):
    return " " * (w - 2) + "--" if x is None or x != x else (f"{x:{w}.2f}" if x < 100 else f"{x:{w}.1f}")


args = sys.argv[1:]
gaps = [5e-3, 1e-4, 1e-5]
original = False
ddp_tag = "hp_tol1e5"                # --tag=<tag>: FilterDDP logs fddp_<system>_T<T>_<arm>_<tag>.log
arms = ("diag", "exact")             # --arms=a,b: the two FilterDDP arms shown (diag, exact, hybrid)
for a in args:
    if a.startswith("--gaps="):
        gaps = [float(g) for g in a.split("=", 1)[1].split(",")]
    elif a == "--primal=original":
        original = True
    elif a.startswith("--tag="):
        ddp_tag = a.split("=", 1)[1]
    elif a.startswith("--arms="):
        arms = tuple(a.split("=", 1)[1].split(","))
    elif a.startswith("--threshold="):          # one primal threshold for every system
        THR = {S: float(a.split("=", 1)[1]) for S in THR}
H = f"{ROOT}/ipopt_hsl/logs"
P = f"{ROOT}/parallel_in_time/logs"


def ddp_log(S, T, arm):
    # --tag=a,b: the first of the tags that has a log for this case
    for tag in ddp_tag.split(","):
        path = f"{P}/fddp_{S}_T{T}_{arm}_{tag}.log"
        if os.path.exists(path):
            return path
    return f"{P}/fddp_{S}_T{T}_{arm}_{ddp_tag.split(',')[0]}.log"


print(f"FilterDDP logs: {ddp_tag}; arms: {arms[0]}, {arms[1]}")
print(f"primal infeasibility: {'original (constraints as posed)' if original else 'internal (each solver`s own)'}; times in seconds")
for gap in gaps:
    print(f"\n=== objective within {gap:g} of the optimum, primal below the system threshold ===")
    print(f"{'case':14s} | {'MA57':>7s} {'MA97':>7s} {'MUMPS':>7s} {'MA57ad':>7s} it | {arms[0] + ' seq':>8s} {'count':>7s} it | {arms[1] + ' seq':>9s} {'count':>7s} {'paired':>7s} it")
    for S, T in CASES:
        tag = "_screen_lf_pv_hp" if S == "large10kC_1ph" else "_screen_lf_hp"
        ref_run = ipopt(f"{H}/ipopt_ma57{tag}_{S}_T{T}.log", gap, THR[S], original)
        if ref_run is None:
            continue
        ref = ref_run[4]
        # MA57ad: MA57 with Ipopt's adaptive barrier update (mu_strategy=adaptive), same reference objective
        ip = [ipopt(f"{H}/ipopt_{ls}{tag}{sfx}_{S}_T{T}.log", gap, THR[S], original, ref)
              for ls, sfx in (("ma57", ""), ("ma97", ""), ("mumps", ""), ("ma57", "_adaptive"))]
        cells = [fmt(x[0]) if x else "     --" for x in ip]
        out = f"{NAME[S] + ' T=' + str(T):14s} | {cells[0]} {cells[1]} {cells[2]} {cells[3]} {str(ip[0][1] if ip[0] else '--'):>2s} |"
        for pos, arm in enumerate(arms):
            run = ddp(ddp_log(S, T, arm))
            if run is None:
                out += f" {'--':>8s} {'--':>7s} --" + (" |" if pos == 0 else "")
                continue
            pr = pairing(S, T) if arm == "exact" else None
            seq, cnt, par, k = ddp_at(run, ref, gap, THR[S], original, pr)
            if pos == 0:
                out += f" {fmt(seq, 8)} {fmt(cnt)} {str(k if k is not None else '--'):>2s} |"
            else:
                out += f" {fmt(seq, 9)} {fmt(cnt)} {fmt(par)} {str(k if k is not None else '--'):>2s}"
        print(out)
print("\n=== full runs: Ipopt to tol 1e-6, FilterDDP to its tolerance 1e-5 ===")
print(f"{'case':14s} | {'MA57':>7s} {'MA97':>7s} {'MUMPS':>7s} {'MA57ad':>7s} it | {arms[0] + ' seq':>8s} it st | {arms[1] + ' seq':>9s} it st")
for S, T in CASES:
    tag = "_screen_lf_pv_hp" if S == "large10kC_1ph" else "_screen_lf_hp"
    ip = [ipopt(f"{H}/ipopt_{ls}{tag}{sfx}_{S}_T{T}.log", 1.0, float("inf"), False)
          for ls, sfx in (("ma57", ""), ("ma97", ""), ("mumps", ""), ("ma57", "_adaptive"))]
    if ip[0] is None:
        continue
    out = f"{NAME[S] + ' T=' + str(T):14s} | " + " ".join(fmt(x[2]) if x else "     --" for x in ip) + f" {ip[0][3]:2d} |"
    for pos, arm in enumerate(arms):
        run = ddp(ddp_log(S, T, arm))
        out += f" {fmt(run['solve_s'], 8 if pos == 0 else 9)} {run['iters']:2d} {run['status']:>2s}" + (" |" if pos == 0 else "") if run else f" {'--':>8s} -- --" + (" |" if pos == 0 else "")
    print(out)
