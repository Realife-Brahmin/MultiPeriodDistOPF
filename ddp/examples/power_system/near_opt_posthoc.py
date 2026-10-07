"""Near-optimality point of a FilterDDP run, recovered from its log.

For runs that went to strict tolerance (no FILTERDDP_NEAR_OPT_REFERENCE at run
time), find the first iteration whose objective is within GAP of the matched
Ipopt objective and whose primal infeasibility is below the system threshold,
and the time at which that iterate was known. Same rule as
near_opt_from_logs.jl: iterate k's residuals are known at the end of pass k's
backward sweep, so

    t_k = setup + sum(total_s of earlier passes) + backward_s of pass k,
    setup = reported solve time - sum(all total_s).

    python near_opt_posthoc.py [--gaps=5e-3,1e-4] <primal_threshold> <ipopt_log> <filterddp_log> [...]
    (pairs of ipopt_log filterddp_log may be repeated; --gaps lists the
    objective gaps to report, default 5e-3)
"""
import re
import sys

GAP = 0.005


def ipopt_objective(path):
    for line in open(path, encoding="utf-8", errors="ignore"):
        m = re.search(r"CENTRAL_IPOPT .* objective=([-0-9.eE+]+) solve_time_s=([0-9.]+)", line)
        if m:
            return float(m.group(1)), float(m.group(2))
    raise SystemExit(f"no CENTRAL_IPOPT line in {path}")


def near_opt(path, ref, primal, gap=GAP):
    rows = {}            # iteration -> (objective, primal_inf)
    passes = []          # (iteration, backward_s, total_s)
    solve_s = None
    row = re.compile(r"^\s+(\d+)\s+([-0-9.]+e[-+]\d+)\s+([0-9.]+e[-+]\d+)\s")
    for line in open(path, encoding="utf-8", errors="ignore"):
        m = row.match(line)
        if m:
            rows.setdefault(int(m.group(1)), (float(m.group(2)), float(m.group(3))))
            continue
        if line.startswith("FILTERDDP_ITER_TIMING"):
            d = dict(kv.split("=") for kv in line.split()[1:] if "=" in kv)
            passes.append((int(d["iteration"]), float(d["backward_s"]), float(d["total_s"])))
            continue
        m = re.match(r"solve complete: ([0-9.]+) s, iterations=(\d+), status=(\S+)", line)
        if m:
            solve_s, final_it, status = float(m.group(1)), int(m.group(2)), m.group(3)
    setup = solve_s - sum(p[2] for p in passes)
    elapsed = 0.0
    seen = set()
    for k, backward, total in passes:
        if k not in seen and k in rows:
            seen.add(k)
            obj, pr = rows[k]
            if abs(obj - ref) / abs(ref) <= gap and pr <= primal:
                return k, setup + elapsed + backward, final_it, solve_s, status
        elapsed += total
    return None, None, final_it, solve_s, status


if __name__ == "__main__":
    argv = sys.argv[1:]
    gaps = [GAP]
    if argv and argv[0].startswith("--gaps="):
        gaps = [float(g) for g in argv.pop(0).split("=", 1)[1].split(",")]
    primal = float(argv[0])
    args = argv[1:]
    for ip, fd in zip(args[0::2], args[1::2]):
        ref, ip_s = ipopt_objective(ip)
        for gap in gaps:
            k, t, final_it, solve_s, status = near_opt(fd, ref, primal, gap)
            name = fd.replace("\\", "/").split("/")[-1] + (f" gap<={gap:g}" if len(gaps) > 1 else "")
            if k is None:
                print(f"NEAR_OPT_POSTHOC {name}: not reached (strict: {final_it} iterations, {solve_s:.1f} s, status {status})")
            else:
                print(f"NEAR_OPT_POSTHOC {name}: iteration={k} elapsed_s={t:.1f} "
                      f"(strict: {final_it} iterations, {solve_s:.1f} s) ipopt_s={ip_s:.1f} ratio={t / ip_s:.1f}")
