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
    SUBSTATION_FILL, SUBSTATION_STROKE, parse_edges, parse_buses, build_tree,
)

ACCENT = "#1f6f78"  # highlight for the representative feeder
FONT = "Times New Roman, Times, serif"


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


def render(system, out_path, size_mm, seed_angle, head_r=3.2, substation_style="black",
           panel="both"):
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
    band_h = 140.0 if panel == "both" else 76.0
    top_pad, plot_h = 26.0, 476.0
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

    pa_r = min(pa_w / 2 - 16.0, plot_h / 2 - 22.0) if show_a else 0.0
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
        f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" '
        f'width="{size_mm}mm" height="{size_mm * VB_H / VB_W:.2f}mm" '
        f'viewBox="0 0 {VB_W:.2f} {VB_H:.2f}" font-family="{FONT}">',
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
        tips = {}
        for i, f in enumerate(sorted(feeders, key=lambda b: int(b) if b.isdigit() else 0)):
            a = math.radians(seed_angle) + 2 * math.pi * i / n
            tips[f] = (pa_cx + pa_r * math.cos(a), pa_cy + pa_r * math.sin(a))

        p.append('<g id="panel-a-feeders" fill="none" stroke-linecap="round">')
        for f, (x2, y2) in tips.items():
            if f == rep:
                continue  # drawn last, on top
            p.append(
                f'<path d="M{pa_cx:.2f},{pa_cy:.2f}L{x2:.2f},{y2:.2f}" '
                f'stroke="{LINE_STROKE}" stroke-width="1.6" stroke-opacity="0.75"/>'
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

    # ---------------- panel (b): one feeder ----------------
    pv_r, bess_r = 11.0, 7.5
    combo_r = max(bess_r, pv_r * 0.78)
    if show_b:
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
    if show_a:
        captions.append((pa_cx, f"{n} feeders", f"{len(order):,} buses"))
    if show_b:
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
    if show_b:
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
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    suffix = {"both": "_two_level", "a": "_top_level", "b": "_feeder"}[args.panel]
    out = (Path(args.out) if args.out
           else REPO / "assets" / "networks" / f"{args.system}{suffix}.svg")
    render(args.system, out, args.size_mm, args.seed_angle,
           head_r=args.head_r, substation_style=args.substation, panel=args.panel)


if __name__ == "__main__":
    main()
