"""Two-level network figure for feeders built by replicating one prototype.

large10kC is a substation plus 102 near-identical feeders (100 of them share an
identical rooted shape: 101 buses, 10 PV, 10 BESS, 40 levels). Drawing all
10,321 buses just renders the same feeder 102 times, which reads as manufactured
rather than informative. This instead draws:

  (a) the macro view -- substation with one stub per feeder, the representative
      one highlighted;
  (b) that representative feeder in full detail, every bus, PV and BESS.

Marker language is shared with make_network_svg.py so all the paper's network
figures read as one family.

Usage:
    python scripts/make_network_svg_two_level.py large10kC_1ph
"""

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_network_svg import (  # noqa: E402
    REPO, PV_FILL, PV_STROKE, BESS_FILL, BESS_STROKE, LINE_STROKE,
    SUBSTATION_FILL, SUBSTATION_STROKE, SVG_NS, namedview,
    parse_edges, parse_buses, build_tree, force_layout,
)

BUS_FILL, BUS_STROKE = "#ffffff", "#3a3a3a"

# Values the feeder renderer computes internally that render() reports back, so
# a figure that is technically correct but too small to read cannot ship quietly.
_diag = {}

ACCENT = "#1f6f78"  # highlight for the representative feeder
FONT = "Times New Roman, Times, serif"

# Alternating stub tones. Deliberately two neighbours of the family grey
# (#7b7b7b) rather than two distinct hues: strong alternation would read as two
# *categories* of feeder, which would be a lie -- 100 of the 102 are identical.
# These give a light-dark rhythm that survives grayscale printing.
STUB_DARK = "#666666"
STUB_LIGHT = "#9a9a9a"


def subtree_nodes(root, children):
    out, stack = [], [root]
    while stack:
        u = stack.pop()
        out.append(u)
        stack.extend(children[u])
    return out


def canonical_shape(root, children):
    """Order-independent encoding of a rooted subtree, for grouping feeders."""
    stack = [(root, False)]
    done = {}
    while stack:
        u, expanded = stack.pop()
        if expanded:
            done[u] = "(" + "".join(sorted(done[k] for k in children[u])) + ")"
        else:
            stack.append((u, True))
            for k in children[u]:
                stack.append((k, False))
    return done[root]


def pick_representative(feeders, children):
    """The feeder whose shape the most other feeders share."""
    groups = defaultdict(list)
    for f in feeders:
        groups[canonical_shape(f, children)].append(f)
    best = max(groups.values(), key=len)
    return sorted(best)[0], len(best)


def tidy_tree(root, children, depth):
    """Classic layered tree: y = depth, x = leaf slot (internal = mean of kids)."""
    order, stack = [], [root]
    while stack:
        u = stack.pop()
        order.append(u)
        stack.extend(children[u])
    order.sort(key=lambda u: depth[u])

    slot = {}
    nxt = [0.0]

    def assign(u):
        kids = children[u]
        if not kids:
            slot[u] = nxt[0]
            nxt[0] += 1.0
            return slot[u]
        vals = [assign(k) for k in kids]
        slot[u] = sum(vals) / len(vals)
        return slot[u]

    sys.setrecursionlimit(10000)
    assign(root)
    d0 = depth[root]
    return {u: (slot[u], depth[u] - d0) for u in slot}


def _render_feeder_detailed(p, rep, rep_nodes, rep_pv, rep_bess, all_edges,
                            parent, children, depth, x0, w, top_pad, plot_h,
                            aspect=2.5, bus_labels=False):
    """One-line diagram of a single feeder: every bus drawn and numbered.

    The earlier schematic did draw all 101 buses, but at r=2.1 against BESS rings
    of r=7.5 -- accurate, yet it read as a ~14-node toy feeder because the buses
    were invisible beside the markers sitting on them. Here buses are real
    circles, numbered, and placed by the spring solver so the feeder sprawls like
    ieee123 instead of marching down a rigid dendrogram.

    Numbers are local (1..n along the feeder). The global ids are 100001..100101
    -- six digits, unreadable at this size, and meaningless outside the feeder.
    """
    sub = set(rep_nodes)
    edges = [(a, b) for a, b in all_edges if a in sub and b in sub]

    # Seed the spring layout with the tree's own shape so it settles quickly.
    lay = tidy_tree(rep, children, depth)
    seed = {u: (lay[u][0] * 2.2, lay[u][1]) for u in rep_nodes}
    pos = force_layout(rep_nodes, edges, seed, iterations=420, seed=7)

    # The spring solver relaxes toward a roughly square blob, which forces the
    # panel into a tall box; at half the text width the bus numbers then land at
    # under 3 pt. Stretching to a wide strip lets the panel run full width, where
    # the same numbers clear 6 pt. Nothing physical is distorted -- these
    # coordinates are synthetic, only the topology carries meaning.
    # Rotate the feeder so its own longest axis runs horizontally. This buys most
    # of the width for free, without distortion. Squashing to the target aspect
    # is then only a mild correction -- doing it by scaling alone crushes any
    # lateral that happens to run across the stretch axis into a stack of
    # overlapping buses.
    cx_ = sum(q[0] for q in pos.values()) / len(pos)
    cy_ = sum(q[1] for q in pos.values()) / len(pos)
    sxx = sum((q[0] - cx_) ** 2 for q in pos.values())
    syy = sum((q[1] - cy_) ** 2 for q in pos.values())
    sxy = sum((q[0] - cx_) * (q[1] - cy_) for q in pos.values())
    theta = 0.5 * math.atan2(2 * sxy, sxx - syy)
    ct, st = math.cos(-theta), math.sin(-theta)
    pos = {u: ((q[0] - cx_) * ct - (q[1] - cy_) * st,
               (q[0] - cx_) * st + (q[1] - cy_) * ct)
           for u, q in pos.items()}

    xs = [q[0] for q in pos.values()]
    ys = [q[1] for q in pos.values()]
    spanx, spany = max(xs) - min(xs), max(ys) - min(ys)
    natural = spanx / max(spany, 1e-9)
    if natural < aspect:  # only ever squash, never stretch further
        sy_ = natural / aspect
        pos = {u: (q[0], q[1] * sy_) for u, q in pos.items()}

    xs = [q[0] for q in pos.values()]
    ys = [q[1] for q in pos.values()]
    spanx, spany = max(xs) - min(xs), max(ys) - min(ys)
    pad = 22.0
    avail_w, avail_h = w - 2 * pad, plot_h - 2 * pad
    k = min(avail_w / max(spanx, 1e-9), avail_h / max(spany, 1e-9))
    offx = x0 + pad + (avail_w - spanx * k) / 2
    offy = top_pad + pad + (avail_h - spany * k) / 2

    def fx(u):
        return offx + (pos[u][0] - min(xs)) * k

    def fy(u):
        return offy + (pos[u][1] - min(ys)) * k

    local = {u: i + 1 for i, u in enumerate(sorted(rep_nodes, key=lambda b: int(b)))}

    # With ieee123's convention the bus is a small dot and the label beside it is
    # the big element, so spacing is governed by label width, not dot size. Cap
    # the type so a 3-digit number still fits between adjacent buses, then take
    # the dot radius from ieee123's own ratio (label 9.88px against r 3.14).
    lens = sorted(math.hypot(fx(u) - fx(parent[u]), fy(u) - fy(parent[u]))
                  for u in rep_nodes if parent[u] in sub)
    med = lens[len(lens) // 2] if lens else 12.0
    widest = len(str(max(local.values())))
    fs = max(4.0, min(11.0, med / (widest * 0.55)))
    bus_r = fs / 3.15
    _diag["bus_label_units"] = fs
    _diag["bus_r"] = bus_r

    # ieee123C ratios, measured from that file rather than eyeballed. Everything
    # is expressed against the bus-dot radius so the whole figure scales as one.
    #   branch stroke 1.50 / bus_r 3.142 = 0.48   (black, opacity 0.83)
    #   PV square side 10.12 / bus_r     = 3.22   (concentric, behind the dot)
    #   BESS ring r     6.29 / bus_r     = 2.00   (concentric, behind the dot)
    #   label size      9.88 / bus_r     = 3.14
    br_w = bus_r * 0.48
    p.append(f'<g id="feeder-branches" fill="none" stroke="#000000" '
             f'stroke-opacity="0.83" stroke-width="{br_w:.2f}" '
             f'stroke-linecap="round" stroke-linejoin="round">')
    seg = [f"M{fx(parent[u]):.2f},{fy(parent[u]):.2f}L{fx(u):.2f},{fy(u):.2f}"
           for u in rep_nodes if parent[u] in sub]
    p.append(f'<path d="{"".join(seg)}"/>')
    p.append("</g>")

    def perp(u):
        """Unit vector across the bus's branch, used to park its label clear."""
        q = parent[u]
        if q not in sub:
            return 0.0, -1.0
        dx, dy = fx(u) - fx(q), fy(u) - fy(q)
        n = math.hypot(dx, dy) or 1.0
        return -dy / n, dx / n

    # PV and BESS sit concentric behind the bus dot, as in ieee123C -- not offset
    # on a leader line, which was my own invention and read as a foreign
    # convention. Largest marker first so each stays visible under the next.
    p.append(f'<g id="feeder-pv" fill="{PV_FILL}" fill-opacity="0.802" '
             f'stroke="{PV_STROKE}" stroke-width="{bus_r * 0.48:.2f}">')
    for u in rep_pv:
        s = bus_r * 3.22
        p.append(f'<rect x="{fx(u) - s / 2:.2f}" y="{fy(u) - s / 2:.2f}" '
                 f'width="{s:.2f}" height="{s:.2f}"/>')
    p.append("</g>")
    p.append(f'<g id="feeder-bess" fill="{BESS_FILL}" fill-opacity="0.55" '
             f'stroke="{BESS_STROKE}" stroke-width="{bus_r * 0.64:.2f}">')
    for u in rep_bess:
        p.append(f'<circle cx="{fx(u):.2f}" cy="{fy(u):.2f}" '
                 f'r="{bus_r * 2.0:.2f}"/>')
    p.append("</g>")

    # Buses and labels follow ieee123C exactly: a solid black dot at opacity 0.83,
    # and the number as separate Times New Roman text beside it -- not a number
    # inside a white circle, which is what made this figure read as a different
    # drawing convention from the ieee123 one-line diagram.
    p.append(f'<g id="feeder-buses" fill="#000000" fill-opacity="0.83" '
             f'stroke="#000000" stroke-width="{bus_r * 0.1:.2f}">')
    for u in rep_nodes:
        p.append(f'<circle cx="{fx(u):.2f}" cy="{fy(u):.2f}" r="{bus_r:.2f}"/>')
    p.append("</g>")

    # Kept in the file even when hidden: display:none leaves the group intact and
    # togglable from Inkscape's Objects panel, so the numbering can be brought
    # back without regenerating. At the widths this panel gets in a two-column
    # figure the labels fall to ~3.5 pt, below anything printable.
    vis = "" if bus_labels else 'style="display:none" '
    p.append(f'<g id="feeder-bus-labels" {vis}font-family="{FONT}" font-size="{fs:.2f}" '
             f'fill="#000000" stroke="{LINE_STROKE}" '
             f'stroke-width="{fs * 0.02:.2f}" stroke-opacity="0.63" '
             f'text-anchor="middle">')
    for u in rep_nodes:
        ox, oy = perp(u)
        # clear the widest concentric marker (the PV square at 3.22 x bus_r)
        off = bus_r * 1.9 + fs * 0.52
        lx = fx(u) - ox * off
        ly = fy(u) - oy * off + fs * 0.34
        p.append(f'<text x="{lx:.2f}" y="{ly:.2f}">{local[u]}</text>')
    p.append("</g>")

    # Feeder head keeps the accent tie back to the highlighted spoke in panel (a).
    p.append(f'<circle cx="{fx(rep):.2f}" cy="{fy(rep):.2f}" r="{bus_r * 3.6:.2f}" '
             f'fill="none" stroke="{ACCENT}" stroke-width="{bus_r * 0.6:.2f}"/>')


def render(system, out_path, size_mm, seed_angle, head_r=3.2, substation_style="black",
           panel="both", alternate=None, labels=None, label_every=1, label_size=None,
           bare=False, feeder_style="detailed", feeder_aspect=2.5, bus_labels=False):
    # Labels are worth it on the standalone star, where there is room; on the
    # two-panel figure the star is under half the width and they would be noise.
    if labels is None:
        labels = panel == "a"
    if alternate is None:
        alternate = True
    raw = REPO / "rawData" / system
    edges = parse_edges(raw / "BranchData.dss")
    pv = parse_buses(raw / "PVSystem.dss", "PVsystem")
    bess = parse_buses(raw / "Storage.dss", "Storage")

    root = "1"
    order, depth, parent, children = build_tree(edges, root)
    feeders = children[root]
    rep, n_alike = pick_representative(feeders, children)

    rep_nodes = subtree_nodes(rep, children)
    rep_pv = [u for u in rep_nodes if u in pv]
    rep_bess = [u for u in rep_nodes if u in bess]
    rep_levels = max(depth[u] for u in rep_nodes) - depth[rep] + 1

    # ---------------- geometry ----------------
    # A bottom band carries the legend and panel labels so nothing overlaps the
    # drawings; both panels share one plot band and fill its full height.
    show_a = panel in ("both", "a")
    show_b = panel in ("both", "b")
    # Panel (a) alone carries no PV/BESS markers, so it needs no legend and only
    # a short caption -- a shallower band than the two-panel figure.
    # --bare drops in-figure captions and the legend so the LaTeX caption can
    # carry them; the bottom band collapses to a thin margin.
    band_h = 14.0 if bare else (140.0 if panel == "both" else 76.0)
    top_pad, plot_h = 26.0, 476.0
    if panel == "b" and feeder_style == "detailed":
        # Track the requested aspect so the panel is a strip when it will run the
        # full text width, and squarer when it sits beside the star.
        plot_h = min(560.0, (540.0 - 24.0) / max(feeder_aspect, 0.2) + 40.0)
    VB_H = top_pad + plot_h + band_h

    if panel == "both":
        VB_W = 1000.0
        pa_w, gutter = 430.0, 60.0
        pb_x0 = pa_w + gutter
        pb_w = VB_W - pb_x0 - 24.0
    elif panel == "a":
        pa_w = 460.0
        VB_W = pa_w
        pb_x0, pb_w = 0.0, 0.0
    else:
        pa_w, gutter = 0.0, 0.0
        VB_W = 540.0
        pb_x0, pb_w = 0.0, VB_W - 24.0

    # Radial labels sit outside the head-dot ring, so the star has to shrink by
    # roughly the widest label to keep the figure inside its box.
    #
    # Thinning with --label-every frees tangential room, so the type grows to use
    # it -- otherwise thinning would fix collisions while leaving the labels just
    # as unreadably small. At label_every=1 the budget is ~12 units, which is
    # ~3.3 pt at 88.9 mm; such a figure needs to be printed ~165 mm wide to reach
    # 6 pt. render() warns when the effective point size falls below that.
    if label_size:
        lab_fs = label_size
    else:
        lab_fs = min(26.0, 13.0 * math.sqrt(max(label_every, 1)))
    lab_room = (len(f"F{len(feeders)}") * lab_fs * 0.56 + 12.0) if labels else 0.0
    pa_r = (min(pa_w / 2 - 16.0 - lab_room, plot_h / 2 - 22.0 - lab_room)
            if show_a else 0.0)
    pa_cx, pa_cy = (VB_W / 2 if panel == "a" else 16.0 + pa_w / 2), top_pad + plot_h / 2

    lay = tidy_tree(rep, children, depth)
    slots = [s for s, _ in lay.values()]
    s_min, s_max = min(slots), max(slots)
    d_max = max(d for _, d in lay.values())
    s_span = max(s_max - s_min, 1e-9)

    pb_y0 = top_pad + 34.0   # headroom for the accent stub above the feeder head
    pb_h = plot_h - 56.0
    pb_iw = pb_w - 26.0

    def bx(s):
        return pb_x0 + 13.0 + (s - s_min) / s_span * pb_iw

    def by(d):
        return pb_y0 + d / max(d_max, 1) * pb_h

    p = [
        f'<svg {SVG_NS} version="1.1" '
        f'width="{size_mm}mm" height="{size_mm * VB_H / VB_W:.2f}mm" '
        f'viewBox="0 0 {VB_W:.2f} {VB_H:.2f}" font-family="{FONT}">',
        namedview(),
        f"<title>{system} two-level network diagram</title>",
        f"<desc>{len(order)} buses, {len(edges)} branches, {len(pv)} PV, "
        f"{len(bess)} BESS across {len(feeders)} feeders. Panel (a) macro view; "
        f"panel (b) feeder {rep} in full ({len(rep_nodes)} buses, {len(rep_pv)} PV, "
        f"{len(rep_bess)} BESS, {rep_levels} levels), the shape shared by "
        f"{n_alike} of {len(feeders)} feeders.</desc>",
    ]

    # ---------------- panel (a): macro ----------------
    n = len(feeders)
    if show_a:
        # Feeder index is positional: sorted by head-bus number, F1 outward from
        # seed_angle. It is a figure label, not a bus name -- the head buses are
        # 2, 1001, 2001, ... 100001.
        ordered = sorted(feeders, key=lambda b: int(b) if b.isdigit() else 0)
        tips, angles, index = {}, {}, {}
        for i, f in enumerate(ordered):
            a = math.radians(seed_angle) + 2 * math.pi * i / n
            tips[f] = (pa_cx + pa_r * math.cos(a), pa_cy + pa_r * math.sin(a))
            angles[f] = a
            index[f] = i + 1

        p.append('<g id="panel-a-feeders" fill="none" stroke-linecap="round">')
        for f, (x2, y2) in tips.items():
            if f == rep:
                continue  # drawn last, on top
            stroke = (STUB_DARK if index[f] % 2 else STUB_LIGHT) if alternate else LINE_STROKE
            p.append(
                f'<path d="M{pa_cx:.2f},{pa_cy:.2f}L{x2:.2f},{y2:.2f}" '
                f'stroke="{stroke}" stroke-width="1.6" stroke-opacity="0.85"/>'
            )
        p.append("</g>")

        # Feeder-head bus at every stub tip. Adjacent tips sit 2*pi*pa_r/n apart
        # (~12 units at n=102), so head_r must stay under ~5 or the ring fuses solid.
        p.append('<g id="panel-a-heads" fill="#000000" fill-opacity="0.83">')
        for f, (x2, y2) in tips.items():
            if f == rep:
                continue
            p.append(f'<circle cx="{x2:.2f}" cy="{y2:.2f}" r="{head_r:.2f}"/>')
        p.append("</g>")

        x2, y2 = tips[rep]
        p.append(
            f'<path d="M{pa_cx:.2f},{pa_cy:.2f}L{x2:.2f},{y2:.2f}" fill="none" '
            f'stroke="{ACCENT}" stroke-width="4.2" stroke-linecap="round"/>'
        )
        p.append(f'<circle cx="{x2:.2f}" cy="{y2:.2f}" r="7" fill="{ACCENT}"/>')

        if substation_style == "red":
            p.append(
                f'<circle cx="{pa_cx:.2f}" cy="{pa_cy:.2f}" r="13" fill="{SUBSTATION_FILL}" '
                f'stroke="{SUBSTATION_STROKE}" stroke-width="2.4"/>'
            )
        else:
            p.append(
                f'<circle cx="{pa_cx:.2f}" cy="{pa_cy:.2f}" r="{head_r:.2f}" '
                f'fill="#000000" fill-opacity="0.83"/>'
            )

        # Feeder labels, rotated to run along their own spoke. At n=102 the
        # tangential budget is 2*pi*pa_r/n ~ 13 units, so horizontal labels would
        # collide outright; running them radially spends that budget on glyph
        # height instead of width. Left-half labels are flipped 180 deg so none
        # of them read upside down.
        if labels:
            p.append(f'<g id="panel-a-labels" font-size="{lab_fs:.1f}">')
            lr = pa_r + head_r + 6.0
            for f in ordered:
                i = index[f]
                if i % label_every and f != rep:
                    continue
                a = angles[f]
                deg = math.degrees(a)
                flip = 90 < (deg % 360) < 270
                lx = pa_cx + lr * math.cos(a)
                ly = pa_cy + lr * math.sin(a)
                rot = deg + 180 if flip else deg
                anchor = "end" if flip else "start"
                fill = ACCENT if f == rep else "#171a18"
                weight = ' font-weight="bold"' if f == rep else ""
                p.append(
                    f'<text x="{lx:.2f}" y="{ly:.2f}" fill="{fill}"{weight} '
                    f'text-anchor="{anchor}" dominant-baseline="middle" '
                    f'transform="rotate({rot:.2f} {lx:.2f} {ly:.2f})">F{i}</text>'
                )
            p.append("</g>")

    # ---------------- panel (b): one feeder ----------------
    pv_r, bess_r = 11.0, 7.5
    combo_r = max(bess_r, pv_r * 0.78)

    if show_b and feeder_style == "detailed":
        _render_feeder_detailed(
            p, rep, rep_nodes, rep_pv, rep_bess, edges, parent, children, depth,
            pb_x0, pb_w, top_pad, plot_h, aspect=feeder_aspect,
            bus_labels=bus_labels,
        )
    elif show_b:
        p.append(
            f'<g id="panel-b-branches" fill="none" stroke="{LINE_STROKE}" '
            f'stroke-width="2.2" stroke-linecap="round">'
        )
        seg = []
        for u in rep_nodes:
            pu = parent[u]
            if pu is None or pu not in lay:
                continue
            s1, d1 = lay[pu]
            s2, d2 = lay[u]
            # elbow routing reads like a one-line diagram
            seg.append(
                f"M{bx(s1):.2f},{by(d1):.2f}L{bx(s1):.2f},{by(d2):.2f}"
                f"L{bx(s2):.2f},{by(d2):.2f}"
            )
        p.append(f'<path d="{"".join(seg)}"/>')
        p.append("</g>")

        # Accent stub + solid dot on the feeder head, mirroring the highlighted
        # spoke in panel (a). Solid accent fill keeps it distinct from a BESS
        # marker, which is pale-filled with the same teal stroke.
        s0, d0 = lay[rep]
        p.append(
            f'<path d="M{bx(s0):.2f},{by(d0) - 30:.2f}L{bx(s0):.2f},{by(d0):.2f}" '
            f'stroke="{ACCENT}" stroke-width="4.2" fill="none" stroke-linecap="round"/>'
        )
        p.append(f'<circle cx="{bx(s0):.2f}" cy="{by(d0):.2f}" r="7" fill="{ACCENT}"/>')

        # Every bus as a faint dot, so "101 buses" is visible and not just asserted.
        p.append(f'<g id="panel-b-buses" fill="{LINE_STROKE}" fill-opacity="0.55">')
        for u in rep_nodes:
            s, d = lay[u]
            p.append(f'<circle cx="{bx(s):.2f}" cy="{by(d):.2f}" r="2.1"/>')
        p.append("</g>")

        both = set(rep_pv) & set(rep_bess)
        p.append(
            f'<g id="panel-b-bess" fill="{BESS_FILL}" stroke="{BESS_STROKE}" '
            f'stroke-width="2.4">'
        )
        for u in rep_bess:
            s, d = lay[u]
            r = combo_r if u in both else bess_r
            p.append(f'<circle cx="{bx(s):.2f}" cy="{by(d):.2f}" r="{r:.2f}"/>')
        p.append("</g>")
        p.append(
            f'<g id="panel-b-pv" fill="{PV_FILL}" stroke="{PV_STROKE}" stroke-width="1.8">'
        )
        for u in rep_pv:
            s, d = lay[u]
            p.append(
                f'<rect x="{bx(s) - pv_r / 2:.2f}" y="{by(d) - pv_r / 2:.2f}" '
                f'width="{pv_r:.2f}" height="{pv_r:.2f}"/>'
            )
        p.append("</g>")

    # ---------------- labels + legend ----------------
    # Two short centred lines per panel; a single long line overflows panel (b),
    # which is only ~490 viewBox units wide.
    y1 = top_pad + plot_h + 28.0
    y2 = top_pad + plot_h + 56.0
    lab_y = top_pad + plot_h + 92.0
    pb_cx = pb_x0 + pb_w / 2
    captions = []
    if show_a and not bare:
        captions.append((pa_cx, f"{n} feeders", f"{len(order):,} buses"))
    if show_b and not bare:
        captions.append((pb_cx, "representative feeder",
                         f"{len(rep_nodes)} buses, {len(rep_pv)} PV, {len(rep_bess)} BESS"))
    for cx, l1, l2 in captions:
        p.append(f'<text x="{cx:.2f}" y="{y1:.2f}" text-anchor="middle" font-size="27">{l1}</text>')
        p.append(f'<text x="{cx:.2f}" y="{y2:.2f}" text-anchor="middle" font-size="27">{l2}</text>')
    if panel == "both":  # (a)/(b) tags are meaningless on a single-panel figure
        p.append(
            f'<text x="{pa_cx:.2f}" y="{lab_y:.2f}" text-anchor="middle" font-size="34">(a)</text>'
        )
        p.append(
            f'<text x="{pb_cx:.2f}" y="{lab_y:.2f}" text-anchor="middle" font-size="34">(b)</text>'
        )

    # Legend runs horizontally across the band so it never crowds either panel.
    # Panel (a) carries no PV/BESS markers, so on its own it needs no legend at
    # all; the substation row is dropped when the substation is a plain bus dot.
    fs = 28.0
    leg_y = top_pad + plot_h + 128.0
    items = []
    if show_b and not bare:
        if substation_style == "red":
            items.append(("sub", "Substation"))
        items += [("pv", "PV"), ("bess", "BESS")]
    widths = [len(lbl) * fs * 0.46 + 46.0 for _, lbl in items]
    cursor = (VB_W - sum(widths)) / 2
    p.append(f'<g id="legend" font-size="{fs}">')
    for (kind, label), wdt in zip(items, widths):
        cx = cursor + 13
        if kind == "sub":
            p.append(
                f'<circle cx="{cx:.2f}" cy="{leg_y - fs * 0.33:.2f}" r="9" '
                f'fill="{SUBSTATION_FILL}" stroke="{SUBSTATION_STROKE}" stroke-width="2.4"/>'
            )
        elif kind == "pv":
            p.append(
                f'<rect x="{cx - pv_r / 2:.2f}" y="{leg_y - fs * 0.33 - pv_r / 2:.2f}" '
                f'width="{pv_r:.2f}" height="{pv_r:.2f}" fill="{PV_FILL}" '
                f'stroke="{PV_STROKE}" stroke-width="1.8"/>'
            )
        else:
            p.append(
                f'<circle cx="{cx:.2f}" cy="{leg_y - fs * 0.33:.2f}" r="{bess_r:.2f}" '
                f'fill="{BESS_FILL}" stroke="{BESS_STROKE}" stroke-width="2.4"/>'
            )
        p.append(f'<text x="{cx + 20:.2f}" y="{leg_y:.2f}">{label}</text>')
        cursor += wdt
    p.append("</g>")
    p.append("</svg>")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(p))
    if show_b and feeder_style == "detailed" and _diag.get("bus_label_units"):
        pt = _diag["bus_label_units"] / VB_W * size_mm * 72.0 / 25.4
        print(f"  bus numbers at {_diag['bus_label_units']:.1f} units "
              f"= {pt:.1f} pt when printed {size_mm:.0f} mm wide")
        if pt < 5.0:
            need = 5.0 * 25.4 / 72.0 * VB_W / _diag["bus_label_units"]
            print(f"  WARNING: bus numbers below ~5 pt. Print this panel "
                  f">= {need:.0f} mm wide, or drop the numbers.")
    if labels and show_a:
        pt = lab_fs / VB_W * size_mm * 72.0 / 25.4
        n_shown = sum(1 for i in range(1, n + 1) if i % label_every == 0)
        print(f"  {n_shown} feeder labels at {lab_fs:.0f} units "
              f"= {pt:.1f} pt when printed {size_mm:.0f} mm wide")
        if pt < 6.0:
            need = 6.0 * 25.4 / 72.0 * VB_W / lab_fs
            print(f"  WARNING: below the ~6 pt floor for print. Either render this "
                  f"figure >= {need:.0f} mm wide, or thin the labels "
                  f"(--label-every 5 grows the type to fit).")
    print(
        f"  {system}: {len(order)} buses over {len(feeders)} feeders; "
        f"representative {rep} = {len(rep_nodes)} buses, {len(rep_pv)} PV, "
        f"{len(rep_bess)} BESS, {rep_levels} levels (shape shared by {n_alike}/{len(feeders)})"
    )
    try:
        shown = out_path.relative_to(REPO)
    except ValueError:
        shown = out_path
    print(f"  wrote {shown} ({out_path.stat().st_size / 1024:.0f} kB)")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("system", help="rawData subfolder, e.g. large10kC_1ph")
    ap.add_argument("--size-mm", type=float, default=88.9, help="IEEE column width")
    ap.add_argument("--seed-angle", type=float, default=-90.0,
                    help="bearing of the first feeder stub in panel (a)")
    ap.add_argument("--head-r", type=float, default=3.2,
                    help="radius of the black feeder-head bus dots in panel (a)")
    ap.add_argument("--substation", choices=("red", "black"), default="black",
                    help="substation glyph: plain black bus dot, or the red marker")
    ap.add_argument("--panel", choices=("both", "a", "b"), default="both",
                    help="'a' emits the macro star alone (the top-level figure)")
    ap.add_argument("--alternate", action=argparse.BooleanOptionalAction, default=None,
                    help="alternate consecutive feeder stubs between two greys (default on)")
    ap.add_argument("--labels", action=argparse.BooleanOptionalAction, default=None,
                    help="label feeder heads F1..Fn, rotated along their spokes "
                         "(default: on for --panel a, off for the two-panel figure)")
    ap.add_argument("--label-every", type=int, default=1, metavar="N",
                    help="label every Nth feeder (the highlighted one is always "
                         "labelled); type size grows to use the freed room")
    ap.add_argument("--label-size", type=float, default=None,
                    help="override label size in viewBox units")
    ap.add_argument("--feeder-style", choices=("detailed", "schematic"), default="detailed",
                    help="'detailed' draws every bus as a numbered circle in an "
                         "organic layout; 'schematic' is the old dendrogram")
    ap.add_argument("--bus-labels", action=argparse.BooleanOptionalAction, default=False,
                    help="render the per-bus numbers in the feeder panel; when off "
                         "they are still written to the SVG but display:none, so "
                         "they stay togglable in Inkscape")
    ap.add_argument("--feeder-aspect", type=float, default=2.5,
                    help="width:height of the feeder panel; 2.5 = wide strip for a "
                         "full-width figure, ~1.2 to sit beside the star")
    ap.add_argument("--bare", action="store_true",
                    help="network only: no in-figure captions or legend, for when "
                         "the LaTeX caption carries that text")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    suffix = {"both": "_two_level", "a": "_top_level", "b": "_low_level"}[args.panel]
    out = (Path(args.out) if args.out
           else REPO / "assets" / "networks" / f"{args.system}{suffix}.svg")
    render(args.system, out, args.size_mm, args.seed_angle,
           head_r=args.head_r, substation_style=args.substation, panel=args.panel,
           alternate=args.alternate, labels=args.labels, label_every=args.label_every,
           label_size=args.label_size, bare=args.bare, feeder_style=args.feeder_style,
           feeder_aspect=args.feeder_aspect, bus_labels=args.bus_labels)


if __name__ == "__main__":
    main()
