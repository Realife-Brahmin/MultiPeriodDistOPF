#!/usr/bin/env python3
"""Memory of a process tree on Linux, sampled from /proc/<pid>/smaps_rollup.

  sample_memory.py sample ROOT_PID OUT.csv [INTERVAL_S]   until ROOT_PID exits
  sample_memory.py summarize OUT.csv SOLVE.log            peak per SC_PHASE interval

Rss counts a page shared by n processes n times; Pss counts it once in total,
so the sum of Pss over the tree is what the tree really occupies.
"""
import os
import sys
import time


def tree(root):
    parent = {}
    for name in os.listdir('/proc'):
        if name.isdigit():
            try:
                with open(f'/proc/{name}/stat') as f:
                    fields = f.read().rsplit(')', 1)[1].split()
                parent[int(name)] = int(fields[1])
            except OSError:
                pass
    found, frontier = {root}, [root]
    while frontier:
        p = frontier.pop()
        for child, par in parent.items():
            if par == p and child not in found:
                found.add(child)
                frontier.append(child)
    return found


def usage(pid):
    rss = pss = 0
    try:
        with open(f'/proc/{pid}/smaps_rollup') as f:
            for line in f:
                if line.startswith('Rss:'):
                    rss = int(line.split()[1])
                elif line.startswith('Pss:'):
                    pss = int(line.split()[1])
    except OSError:
        pass
    return rss, pss


def sample(root, out, interval):
    with open(out, 'w') as f:
        f.write('epoch_s,processes,sum_rss_kb,sum_pss_kb,root_rss_kb,root_pss_kb\n')
        while os.path.exists(f'/proc/{root}'):
            pids = tree(root)
            values = {p: usage(p) for p in pids}
            f.write('%.3f,%d,%d,%d,%d,%d\n' % (time.time(), len(pids), sum(v[0] for v in values.values()),
                                               sum(v[1] for v in values.values()), *values.get(root, (0, 0))))
            f.flush()
            time.sleep(interval)


def julia_pids(uid):
    """Every julia process of this user on this node (workers started by srun are not in the tree)."""
    found = []
    for name in os.listdir('/proc'):
        if name.isdigit():
            try:
                if os.stat(f'/proc/{name}').st_uid != uid:
                    continue
                with open(f'/proc/{name}/comm') as f:
                    if 'julia' in f.read():
                        found.append(int(name))
            except OSError:
                pass
    return found


def sample_node(out, interval, stop):
    """Sum over all of this user's julia processes on the node, until the file `stop` exists."""
    uid = os.getuid()
    with open(out, 'w') as f:
        f.write('epoch_s,processes,sum_rss_kb,sum_pss_kb,root_rss_kb,root_pss_kb\n')
        while not os.path.exists(stop):
            values = [usage(p) for p in julia_pids(uid)]
            f.write('%.3f,%d,%d,%d,0,0\n' % (time.time(), len(values), sum(v[0] for v in values), sum(v[1] for v in values)))
            f.flush()
            time.sleep(interval)


def summarize(samples, log):
    marks = {}
    with open(log) as f:
        for line in f:
            if line.startswith('SC_PHASE '):
                _, name, stamp = line.split()
                marks[name] = float(stamp)
    rows = []
    with open(samples) as f:
        next(f)
        for line in f:
            rows.append([float(x) for x in line.split(',')])
    if rows:
        print('MEMORY %-28s samples=%3d processes=%3d peak_pss_gb=%.2f peak_rss_sum_gb=%.2f' % (
            'whole job', len(rows), max(r[1] for r in rows), max(r[3] for r in rows) / 2 ** 20, max(r[2] for r in rows) / 2 ** 20))
    phases = [('centralized, timed solve', 'central_timed_begin', 'central_timed_end'),
              ('centralized, build to end', 'central_build_begin', 'central_timed_end'),
              ('tADMM, warm-up pass', 'tadmm_warmup_begin', 'tadmm_warmup_end'),
              ('tADMM, timed pass', 'tadmm_timed_begin', 'tadmm_timed_end')]
    for label, a, b in phases:
        if a in marks and b in marks:
            inside = [r for r in rows if marks[a] <= r[0] <= marks[b]]
            if inside:
                print('MEMORY %-28s samples=%3d processes=%3d peak_pss_gb=%.2f peak_rss_sum_gb=%.2f' % (
                    label, len(inside), max(r[1] for r in inside), max(r[3] for r in inside) / 2 ** 20,
                    max(r[2] for r in inside) / 2 ** 20))
            else:
                print('MEMORY %-28s no samples (phase shorter than the sampling interval)' % label)


if __name__ == '__main__':
    if sys.argv[1] == 'sample':
        sample(int(sys.argv[2]), sys.argv[3], float(sys.argv[4]) if len(sys.argv) > 4 else 2.0)
    elif sys.argv[1] == 'sample-node':
        # sample_memory.py sample-node OUT.csv INTERVAL_S STOPFILE
        sample_node(sys.argv[2], float(sys.argv[3]), sys.argv[4])
    else:
        summarize(sys.argv[2], sys.argv[3])
