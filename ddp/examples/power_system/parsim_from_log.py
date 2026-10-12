"""Time of a FilterDDP run with one worker per period, from a FILTERDDP_PARSIM=1 log.

The run itself is sequential and unchanged. Each stage's backward work is
timed in three parts (backward_pass.jl): `pre` needs nothing from stage t+1,
`seq` needs its value function, `post` needs only this stage's battery
solution. With one worker per period the sweep takes

    max_t(pre) + sum_t(seq) + max_t(post)

instead of the sum of all three. In the forward pass the coordinator's part is
a recursion on the states (rec); each period's policy solve, directions and
ratio test (dir), and each trial's evaluations, then run at once and count as
their maximum over periods. Everything else in an iteration (filter, barrier
update, logging) is kept as measured.

This is the accounting of the tADMM paper: per-period work is run in sequence
and the slowest period is counted. It assumes one core per period and ignores
communication and memory contention.

    python parsim_from_log.py <filterddp_log> [...]
"""
import re
import sys

KV = re.compile(r"(\w+)=([-+0-9.eE]+)")


def fields(line):
    return {k: float(v) for k, v in KV.findall(line)}


def analyse(path):
    back, fwd = [], []
    rows = []
    reported = None
    near = None
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith("FILTERDDP_PARSIM_BACKWARD"):
            d = fields(line)
            # paired_pass.jl: stage times without the collections that fell in them;
            # a worker per period collects 1/stages of the pass's garbage
            d["gc_par_s"] = (d.get("pre_gc_s", 0.0) + d.get("post_gc_s", 0.0)) * (1.0 - 1.0 / max(d.get("stages", 1.0), 1.0))
            back.append(d)
        elif line.startswith("FILTERDDP_PARSIM_FORWARD"):
            fwd.append(fields(line))
        elif line.startswith("FILTERDDP_ITER_TIMING"):
            f = fields(line)
            b = {k: sum(x.get(k, 0.0) for x in back) for k in
                 ("pre_sum_s", "pre_max_s", "seq_sum_s", "post_sum_s", "post_max_s", "stages",
                  "pair_sum_s", "pair_max_s", "pre_gc_s", "post_gc_s", "pair_gc_s", "gc_par_s")} if back else None
            w = {k: sum(x[k] for x in fwd) for k in
                 ("dir_sum_s", "dir_max_s", "rec_s", "trial_sum_s", "trial_max_s", "trials")} if fwd else None
            err = max((x["rec_err"] for x in fwd), default=0.0)
            rows.append((f, b, w, err))
            back, fwd = [], []
        elif "solve complete:" in line:
            reported = float(re.search(r"solve complete:\s*([0-9.]+)", line).group(1))
        elif line.startswith("FILTERDDP_NEAR_OPT"):
            near = fields(line)
    if not rows or rows[0][1] is None:
        print(f"{path}: no FILTERDDP_PARSIM lines")
        return
    tot = sum(r[0]["total_s"] for r in rows)
    bsum = sum(r[0]["backward_s"] for r in rows)
    fsum = sum(r[0]["forward_s"] for r in rows)
    B = {k: sum(r[1][k] for r in rows if r[1]) for k in rows[0][1]}
    W = {k: sum(r[2][k] for r in rows if r[2]) for k in
         ("dir_sum_s", "dir_max_s", "rec_s", "trial_sum_s", "trial_max_s", "trials")}
    # pair_*: pairwise combinations of a paired pass (paired_pass.jl), one worker each
    save_b = ((B["pre_sum_s"] - B["pre_max_s"]) + (B["post_sum_s"] - B["post_max_s"])
              + (B["pair_sum_s"] - B["pair_max_s"]) + B["gc_par_s"])
    save_f = (W["dir_sum_s"] - W["dir_max_s"]) + (W["trial_sum_s"] - W["trial_max_s"])
    setup = (reported - tot) if reported is not None else float("nan")
    n = len(rows)
    stages = B["stages"]
    err = max(r[3] for r in rows)
    print(f"PARSIM {path}")
    print(f"  passes={n} stage_solves={int(stages)} reported_solve_s={reported} setup_s={setup:.2f}")
    print(f"  backward: measured {bsum:.2f} s = pre {B['pre_sum_s']:.2f} + seq {B['seq_sum_s']:.2f} + post {B['post_sum_s']:.2f}"
          f" (+ {bsum - B['pre_sum_s'] - B['seq_sum_s'] - B['post_sum_s']:.2f} outside the stages)")
    print(f"            per stage: pre {1e3 * B['pre_sum_s'] / stages:.1f} ms, seq {1e3 * B['seq_sum_s'] / stages:.1f} ms,"
          f" post {1e3 * B['post_sum_s'] / stages:.1f} ms")
    print(f"            one worker per period: {bsum - save_b:.2f} s"
          f" (max pre {B['pre_max_s']:.2f} + seq {B['seq_sum_s']:.2f} + max post {B['post_max_s']:.2f})")
    if B["pre_gc_s"] + B["post_gc_s"] + B["pair_gc_s"] > 0:
        print(f"            collections: {B['pre_gc_s']:.2f} s in step 1, {B['pair_gc_s']:.2f} s in the rounds, {B['post_gc_s']:.2f} s in step 3;"
              f" counted {B['pre_gc_s'] + B['post_gc_s'] - B['gc_par_s'] + B['pair_gc_s']:.2f} s (a worker collects 1/stages of steps 1 and 3),"
              f" all in sequence would add {B['gc_par_s']:.2f} s")
    print(f"  forward:  measured {fsum:.2f} s; directions {W['dir_sum_s']:.2f} (max {W['dir_max_s']:.2f}), state recursion {W['rec_s']:.3f},"
          f" {int(W['trials'])} trials {W['trial_sum_s']:.2f} (max {W['trial_max_s']:.2f}); recursion error {err:.1e}")
    print(f"            one worker per period: {fsum - save_f:.2f} s")
    sim = (reported if reported is not None else tot) - save_b - save_f
    print(f"  solve: measured {reported if reported is not None else tot:.1f} s -> one worker per period {sim:.1f} s"
          f" ({(reported if reported is not None else tot) / sim:.1f}x); per pass {tot / n:.2f} -> {(tot - save_b - save_f) / n:.2f} s")
    if near:
        print(f"  near-optimality at iteration {int(near['iteration'])}: {near['elapsed_s']:.1f} s -> {near['elapsed_s'] - save_b - save_f:.1f} s")


for p in sys.argv[1:]:
    analyse(p)
