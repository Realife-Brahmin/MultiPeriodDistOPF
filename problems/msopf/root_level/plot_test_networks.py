"""Topology figures of the MSOPF test decks (2026-10-24).

    python problems/msopf/root_level/plot_test_networks.py [deck ...]

For each deck in test_networks/<deck>/ (written by convert_to_1ph.jl), draws every line from the published bus
coordinates of its original system, coloured by the substation that feeds it (breadth-first from each source over
the deck's lines, ties open), the substations (stars, labelled), the normally-open ties (dashed; ties between
substations heavier), and the PV units and batteries. Secondary (customer) buses have no published coordinates of
their own; they sit at their service transformer, and the title counts them. Writes
problems/msopf/results/test_networks/<deck>.png and a gallery of all decks drawn, gallery.png.
"""
import collections
import math
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
NET = os.path.join(ROOT, "test_networks")
SRC = os.path.join(NET, "sources")
OUT = os.path.join(ROOT, "problems", "msopf", "results", "test_networks")

P1U = os.path.join(SRC, "smartds_sfo_p1u", "Buscoords.dss")
DECKS = {  # deck -> (coordinates file, coordinates are longitude/latitude, title)
    "epri_k1xf_1ph": (os.path.join(SRC, "epri_test_circuits", "epri_dpv", "K1", "Buscoords.dss"), False, "EPRI K1"),
    "epri_ckt24xf_1ph": (os.path.join(SRC, "epri_test_circuits", "ckt24", "buscoords_ckt24.dss"), False, "EPRI Ckt24"),
    "ieee8500xf_1ph": (os.path.join(SRC, "ieee8500", "Buscoords.dss"), False, "IEEE 8500"),
    "ieee9500xf_1ph": (os.path.join(NET, "ieee9500_source", "base", "BusCoords.dss"), False, "IEEE 9500"),
    "smartds_p1u_xf_1ph": (P1U, True, "SMART-DS SFO P1U (region)"),
    "smartds_p1u_hs0_1ph": (P1U, True, "SMART-DS P1U substation hs0"),
    "smartds_p1u_hs1_1ph": (P1U, True, "SMART-DS P1U substation hs1"),
    "smartds_p1u_hs2_1ph": (P1U, True, "SMART-DS P1U substation hs2"),
    "smartds_p1u_hs3_1ph": (P1U, True, "SMART-DS P1U substation hs3"),
}
GALLERY = ["epri_k1xf_1ph", "epri_ckt24xf_1ph", "ieee8500xf_1ph", "ieee9500xf_1ph", "smartds_p1u_xf_1ph"]
COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#17becf", "#bcbd22"]

LINE = re.compile(r"^New\s+Line\.(\S+)\s.*?Bus1=(\S+?)(?:\.\d+)*\s.*?Bus2=(\S+?)(?:\.\d+)*\s.*?r1=(\S+)\s+x1=(\S+)", re.I)
BUS1 = re.compile(r"\bBus1?=(\S+?)(?:\.\d+)*(?:\s|$)", re.I)


def read_coords(path):
    xy = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.split("!")[0].split("//")[0].strip()
        if not line:
            continue
        f = re.split(r"[,\s]+", line)
        if len(f) >= 3:
            try:
                xy[f[0].lower()] = (float(f[1]), float(f[2]))
            except ValueError:
                pass
    return xy


def read_deck(deck):
    d = os.path.join(NET, deck)
    src = []
    for line in open(os.path.join(d, "Master.dss"), encoding="utf-8"):
        m = re.match(r"New\s+Vsource\.(\S+)\s.*?Bus1=(\S+?)(?:\.\d+)*\s", line, re.I)
        if m:
            src.append((m.group(1), m.group(2).lower()))
    branches = []
    for line in open(os.path.join(d, "BranchData.dss"), encoding="utf-8"):
        m = LINE.match(line)
        if m:
            branches.append((m.group(1).lower(), m.group(2).lower(), m.group(3).lower()))
    ties = []
    for line in open(os.path.join(d, "Ties.dss"), encoding="utf-8"):
        m = LINE.match(line)
        if m:
            ties.append((m.group(1).lower(), m.group(2).lower(), m.group(3).lower()))
    def buses(fname, prefix):
        out = []
        for line in open(os.path.join(d, fname), encoding="utf-8"):
            if re.match(r"New\s+" + prefix + r"\.", line, re.I):
                m = BUS1.search(line)
                kw = re.search(r"(?:kW|Pmpp|kWrated)=(\S+)", line, re.I)
                out.append((m.group(1).lower(), float(kw.group(1)) if kw else 0.0))
        return out
    loads = buses("Loads.dss", "Load")
    pv = buses("PVSystem.dss", "PVSystem")
    bt = buses("Storage.dss", "Storage")
    return src, branches, ties, loads, pv, bt


def feeders(src, branches):
    adj = collections.defaultdict(list)
    for _, a, b in branches:
        adj[a].append(b)
        adj[b].append(a)
    lab = {}
    for i, (_, b) in enumerate(src):
        if b in lab:
            continue
        lab[b] = i
        q = collections.deque([b])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v not in lab:
                    lab[v] = i
                    q.append(v)
    return lab, adj


def draw(ax, deck, compact=False):
    cpath, lonlat, title = DECKS[deck]
    xy = read_coords(cpath)
    src, branches, ties, loads, pv, bt = read_deck(deck)
    lab, adj = feeders(src, branches)
    svc = {b: a for n, a, b in branches if n.startswith("svc_")}         # secondary bus -> its primary bus
    # a bus without coordinates of its own: its service transformer's primary bus, else the nearest bus (in lines)
    # that has some, within 10 lines
    def pos(b):
        if b in xy:
            return xy[b]
        if b in svc and svc[b] in xy:
            return xy[svc[b]]
        seen, q = {b}, collections.deque([(b, 0)])
        while q:
            u, k = q.popleft()
            if u in xy:
                return xy[u]
            if k < 10:
                for v in adj.get(u, []):
                    if v not in seen:
                        seen.add(v)
                        q.append((v, k + 1))
        return None
    # a substation: its source bus, except that one substation with several transformers (P1U hs0) is one
    def station(i):
        m = re.search(r"p1uhs\d", src[i][1])
        return m.group(0) if m else src[i][1]
    if lonlat:
        lat0 = sum(y for _, y in xy.values()) / len(xy)
        kx, ky = 111.32 * math.cos(math.radians(lat0)), 110.57             # km per degree
    else:
        kx = ky = 1.0
    def p(b):
        q = pos(b)
        return None if q is None else (q[0] * kx, q[1] * ky)
    segs = collections.defaultdict(list)
    drawn = 0
    for n, a, b in branches:
        if n.startswith("svc_") or n.startswith("xfmr_"):
            continue
        pa, pb = p(a), p(b)
        if pa and pb:
            segs[lab.get(a, -1)].append([pa, pb])
            drawn += 1
    nprim = sum(1 for n, a, b in branches if not n.startswith("svc_"))
    for i, ss in segs.items():
        ax.add_collection(LineCollection(ss, colors=COLORS[i % len(COLORS)] if i >= 0 else "0.6", linewidths=0.6 if compact else 0.8))
    inter = 0
    shown = set()
    for n, a, b in ties:
        pa, pb = p(a), p(b)
        if not (pa and pb):
            continue
        between = a in lab and b in lab and station(lab[a]) != station(lab[b])
        inter += between
        lbl = "open tie between substations" if between else "open tie within a substation"
        ax.plot([pa[0], pb[0]], [pa[1], pb[1]], ls="--", color="k" if between else "0.45", lw=1.6 if between else 0.7, zorder=4,
                label=None if lbl in shown else lbl)
        shown.add(lbl)
    for kind, pts, mk, col in (("PV", pv, "^", "#ff9f1c"), ("battery", bt, "s", "#2a9d8f")):
        q = [p(b) for b, _ in pts]
        q = [x for x in q if x]
        if q:
            ax.scatter([x for x, _ in q], [y for _, y in q], s=6 if compact else 9, marker=mk, color=col, zorder=5, linewidths=0, label=kind)
    placed = []
    for i, (name, b) in enumerate(src):
        # the source bus is the substation transformer's high-voltage side; draw it at its low-voltage bus
        lv = next((y for n, x, y in branches if n.startswith("xfmr_") and x == b), b)
        q = p(lv) or p(b)
        if q:
            ax.scatter([q[0]], [q[1]], s=170 if compact else 220, marker="*", color=COLORS[i % len(COLORS)], edgecolors="k", linewidths=0.8, zorder=6)
            if not compact:
                near = sum(1 for r in placed if math.hypot(r[0] - q[0], r[1] - q[1]) < 1e-3 * (1 + abs(q[0])))   # stacked labels
                tag = name + (f" ({station(i)})" if station(i) != src[i][1] else "")
                ax.annotate(tag, q, xytext=(6, 5 + 11 * near), textcoords="offset points", fontsize=8, weight="bold")
            placed.append(q)
    nsec = len(svc)
    load_mw = sum(k for _, k in loads) / 1000
    nb = len({x for _, a, b in branches for x in (a, b)})
    nsrcbus = len({b for _, b in src})
    stats = (f"{nb - nsrcbus:,} buses ({nb - nsrcbus - nsec:,} primary, {nsec:,} secondary) | {len(src)} source(s) | "
             f"{len(ties)} open ties ({inter} between substations)\n"
             f"load {load_mw:.1f} MW per phase | {len(pv)} PV, {len(bt)} batteries | lines drawn {drawn}/{nprim - len(src)}")
    ax.set_title(f"{title}  ({deck})\n{stats}", fontsize=7 if compact else 9)
    ax.set_aspect("equal")
    ax.autoscale_view()
    ax.set_xticks([])
    ax.set_yticks([])
    if lonlat:
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        L = 1.0 if (x1 - x0) > 3 else 0.5
        ax.plot([x0 + 0.05 * (x1 - x0), x0 + 0.05 * (x1 - x0) + L], [y0 + 0.05 * (y1 - y0)] * 2, color="k", lw=2)
        ax.annotate(f"{L:g} km", (x0 + 0.05 * (x1 - x0), y0 + 0.07 * (y1 - y0)), fontsize=7)
    return drawn, nprim - len(src)


def main(decks):
    os.makedirs(OUT, exist_ok=True)
    for deck in decks:
        fig, ax = plt.subplots(figsize=(8, 8))
        drawn, total = draw(ax, deck)
        ax.legend(loc="upper right", fontsize=8, markerscale=2)
        fig.tight_layout()
        f = os.path.join(OUT, deck + ".png")
        fig.savefig(f, dpi=200)
        plt.close(fig)
        print(f"{deck}: {drawn}/{total} lines drawn -> {os.path.relpath(f, ROOT)}")
    gal = [d for d in GALLERY if d in decks]
    if len(gal) > 1:
        n = len(gal)
        cols = 3 if n > 4 else n
        rows = math.ceil(n / cols)
        fig, axs = plt.subplots(rows, cols, figsize=(5.2 * cols, 5.4 * rows))
        axs = list(axs.flat) if hasattr(axs, "flat") else [axs]
        for ax, deck in zip(axs, gal):
            draw(ax, deck, compact=True)
        for ax in axs[len(gal):]:
            ax.axis("off")
        fig.tight_layout()
        f = os.path.join(OUT, "gallery.png")
        fig.savefig(f, dpi=170)
        plt.close(fig)
        print("gallery ->", os.path.relpath(f, ROOT))


if __name__ == "__main__":
    main(sys.argv[1:] or list(DECKS))
