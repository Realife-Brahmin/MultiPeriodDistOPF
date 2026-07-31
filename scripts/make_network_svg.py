"""Render a single-phase OpenDSS radial feeder to SVG.

Large synthetic feeders (large10kC, ieee2552C) ship without bus coordinates, so
the layout is computed from the topology: a radial tree rooted at the substation,
with angular sector width allocated in proportion to subtree leaf count.

Marker language matches assets/networks/ieee123C_1ph.svg (the hand-drawn ieee123C
diagram) so the figures read as one family in the paper:
  PV   -> gold square   fill #d6a21a / stroke #6e4b1f
  BESS -> pale circle   fill #e8f0f1 / stroke #1f6f78
  line -> grey #7b7b7b

Usage:
    python scripts/make_network_svg.py large10kC_1ph
    python scripts/make_network_svg.py ieee2552C_1ph --arc 320 --rotate 100
"""

import argparse
import math
import re
import sys
from collections import defaultdict, deque
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

PV_FILL, PV_STROKE = "#d6a21a", "#6e4b1f"
BESS_FILL, BESS_STROKE = "#e8f0f1", "#1f6f78"
LINE_STROKE = "#7b7b7b"
SUBSTATION_FILL, SUBSTATION_STROKE = "#c51919", "#5c0000"


# The sodipodi URI must be sodipodi-0.dtd exactly; sodipodi-0.0.dtd is a
# different namespace and Inkscape then treats <sodipodi:namedview> as an
# unknown element ("unknown type: sodipodi0:namedview") and ignores it.
SVG_NS = (
    'xmlns="http://www.w3.org/2000/svg" '
    'xmlns:svg="http://www.w3.org/2000/svg" '
    'xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd" '
    'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"'
)


# Per-system ink for buses and branches: the paper's system colour (ieeeC/medC/
# lgC) mixed 45% into black. Applied to the network's own structure rather than
# a background wash, and kept dark on purpose -- at full saturation med2522's
# gold would be indistinguishable from the PV squares and ieee123's blue from
# the teal BESS rings. Dark enough to read as near-black, tinted enough to
# identify the system.
SYSTEM_BASE = {
    "ieee123": (31, 81, 134),    # ieeeC #1F5186  blue
    "ieee2552": (184, 134, 11),  # medC  #B8860B  gold
    "large10k": (34, 120, 34),   # lgC   #227822  green
}
INK_MIX = 0.45


def _hex(rgb):
    return "#%02x%02x%02x" % tuple(max(0, min(255, int(round(v)))) for v in rgb)


def ink_for(system, mix=INK_MIX):
    """System ink for buses/branches, or None if the system is unrecognised."""
    for key, rgb in SYSTEM_BASE.items():
        if system.lower().startswith(key):
            return _hex(tuple(c * mix for c in rgb))
    return None


def lighten(hex_colour, amount):
    """Blend a colour toward white; used for the alternating feeder stubs."""
    r, g, b = (int(hex_colour[i:i + 2], 16) for i in (1, 3, 5))
    return _hex(tuple(c + (255 - c) * amount for c in (r, g, b)))


def namedview(page="#ffffff", desk="#c8c8c8"):
    """Inkscape document view settings, so the editing canvas opens white.

    Without this element Inkscape supplies its own defaults, which follow a dark
    UI theme -- and any page colour set by hand is lost the next time the file is
    regenerated. pageopacity matters as much as pagecolor: at 0 the page is
    transparent and the (dark) desk shows straight through it, so a white
    pagecolor alone still renders dark.

    pageopacity=1 also makes exports land on opaque white rather than
    transparent, which is what these figures want -- they are going onto a white
    page, and a transparent background renders dark in some PDF viewers.
    """
    return (
        f'<sodipodi:namedview id="namedview" pagecolor="{page}" '
        f'bordercolor="#666666" borderopacity="1.0" '
        f'inkscape:pageopacity="1" inkscape:pagecheckerboard="0" '
        f'inkscape:showpageshadow="0" inkscape:deskcolor="{desk}" '
        f'inkscape:document-units="mm"/>'
    )


def parse_edges(path):
    """Return [(bus1, bus2)] from `New Line.` records."""
    edges = []
    for raw in path.read_text().splitlines():
        if not raw.strip().lower().startswith("new line"):
            continue
        b1 = re.search(r"Bus1=([^\s.]+)", raw, re.I)
        b2 = re.search(r"Bus2=([^\s.]+)", raw, re.I)
        if b1 and b2:
            edges.append((b1.group(1), b2.group(1)))
    return edges


def parse_buses(path, kind):
    """Return the set of buses hosting a `New <kind>.` element."""
    if not path.exists():
        return set()
    buses = set()
    pat = re.compile(rf"^\s*New\s+{kind}\.", re.I)
    for raw in path.read_text().splitlines():
        if not pat.match(raw):
            continue
        m = re.search(r"Bus1=([^\s.]+)", raw, re.I)
        if m:
            buses.add(m.group(1))
    return buses


def build_tree(edges, root):
    adj = defaultdict(list)
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    if root not in adj:
        sys.exit(f"root bus {root!r} not present in branch data")

    depth = {root: 0}
    parent = {root: None}
    children = defaultdict(list)
    order = []
    q = deque([root])
    while q:
        u = q.popleft()
        order.append(u)
        for v in adj[u]:
            if v not in depth:
                depth[v] = depth[u] + 1
                parent[v] = u
                children[u].append(v)
                q.append(v)

    stranded = set(adj) - set(depth)
    if stranded:
        print(f"  warning: {len(stranded)} buses unreachable from root, dropped")
    return order, depth, parent, children


def radial_layout(order, depth, children, root, arc_deg, rotate_deg, dr,
                  depth_scale="linear"):
    """Angular sector per subtree proportional to its leaf count.

    `depth_scale` maps hop count to radius. large10kC has a single 120-bus
    lateral reaching depth 101 while the other 101 feeders stop at depth 40;
    "sqrt" compresses that outlier so the dense core keeps most of the canvas.
    Both mappings are monotone in hop distance from the substation.
    """
    leaves = {}
    for u in reversed(order):  # BFS order reversed => children before parents
        kids = children[u]
        leaves[u] = 1 if not kids else sum(leaves[k] for k in kids)

    if depth_scale == "sqrt":
        radius_of = lambda d: math.sqrt(d) * dr
    elif depth_scale == "linear":
        radius_of = lambda d: d * dr
    else:
        raise ValueError(f"unknown depth scale {depth_scale!r}")

    arc = math.radians(arc_deg)
    rot = math.radians(rotate_deg)
    span = {root: (rot, rot + arc)}
    angle = {}
    radius = {}
    pos = {}

    for u in order:  # BFS: parent's span already assigned
        lo, hi = span[u]
        angle[u] = 0.5 * (lo + hi)
        r = radius_of(depth[u])
        radius[u] = r
        pos[u] = (r * math.cos(angle[u]), r * math.sin(angle[u]))

        kids = children[u]
        if not kids:
            continue
        total = sum(leaves[k] for k in kids)
        cursor = lo
        for k in kids:
            width = (hi - lo) * leaves[k] / total
            span[k] = (cursor, cursor + width)
            cursor += width
    return pos, angle, radius


def force_layout(order, edges, pos0, iterations=300, seed=0):
    """Fruchterman-Reingold, seeded from the radial layout.

    Deep feeders (ieee2552C: 275 levels, ~9 buses per level) look wrong under a
    radial tree -- the angular sweeps dominate and the interior goes empty. A
    spring layout instead spreads laterals into the organic sprawl that reads
    like a real distribution feeder. O(n^2) repulsion, fine at a few thousand
    buses; the radial seed keeps the iteration count low.
    """
    import numpy as np

    idx = {u: i for i, u in enumerate(order)}
    n = len(order)
    p = np.array([pos0[u] for u in order], dtype=np.float32)
    # normalise the seed into a unit square
    p -= p.min(axis=0)
    span = p.max()
    if span > 0:
        p /= span
    rng = np.random.default_rng(seed)
    p += rng.normal(0, 1e-4, p.shape).astype(np.float32)  # break exact ties

    e = np.array([(idx[a], idx[b]) for a, b in edges if a in idx and b in idx],
                 dtype=np.int32)
    src, dst = e[:, 0], e[:, 1]

    k = float(np.sqrt(1.0 / n))
    temp = 0.10
    cool = temp / (iterations + 1)

    for _ in range(iterations):
        diff = p[:, None, :] - p[None, :, :]          # n x n x 2
        dist = np.sqrt((diff ** 2).sum(-1))
        np.fill_diagonal(dist, np.inf)
        rep = (k * k) / dist
        disp = (diff / dist[:, :, None] * rep[:, :, None]).sum(axis=1)

        d_edge = p[src] - p[dst]
        len_edge = np.maximum(np.sqrt((d_edge ** 2).sum(-1)), 1e-9)
        att = (len_edge ** 2) / k
        contrib = d_edge / len_edge[:, None] * att[:, None]
        np.add.at(disp, src, -contrib)
        np.add.at(disp, dst, contrib)

        mag = np.maximum(np.sqrt((disp ** 2).sum(-1)), 1e-9)
        p += disp / mag[:, None] * np.minimum(mag, temp)[:, None]
        temp -= cool

    return {u: (float(p[i, 0]), float(p[i, 1])) for u, i in idx.items()}


def legend_svg(vb_w, vb_h, corner, pv_size, bess_r, line_w,
               n_pv_only, n_bess_only, n_both, combo_r, hidden=False):
    """Legend block, sized in viewBox units, in Times to match the paper body font.

    Rows with a zero count are dropped, so ieee2552C (every DER bus hosts both)
    shows a single "PV + BESS" row rather than three misleading ones.
    """
    scale = 2.6  # legend markers read larger than the in-plot ones
    pv_s = pv_size * scale
    be_r = bess_r * scale
    co_r = combo_r * scale
    slot = max(be_r, co_r, pv_s / 2)
    fs = 26.0
    row = max(fs * 1.55, slot * 2.3)
    pad = 0.03 * vb_w

    rows = [("substation", "Substation")]
    if n_pv_only:
        rows.append(("pv", f"PV ({n_pv_only:,})"))
    if n_bess_only:
        rows.append(("bess", f"BESS ({n_bess_only:,})"))
    if n_both:
        rows.append(("both", f"PV + BESS ({n_both:,})"))

    x = pad if corner.endswith("left") else vb_w - pad - 0.34 * vb_w
    y = pad + slot if corner.startswith("top") else vb_h - pad - len(rows) * row

    # Hidden rather than omitted: display:none keeps the group in the file and
    # togglable from Inkscape's Objects panel, so the legend can be brought back
    # without regenerating. The LaTeX caption carries this information instead.
    vis = 'style="display:none" ' if hidden else ""
    out = [f'<g id="legend" {vis}font-family="Times New Roman, Times, serif" '
           f'font-size="{fs}">']
    for i, (kind, label) in enumerate(rows):
        cx, cy = x + slot, y + i * row
        if kind == "substation":
            out.append(
                f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{be_r * 1.15:.2f}" '
                f'fill="{SUBSTATION_FILL}" stroke="{SUBSTATION_STROKE}" '
                f'stroke-width="{line_w * 2:.2f}"/>'
            )
        else:
            if kind in ("bess", "both"):
                r = co_r if kind == "both" else be_r
                out.append(
                    f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r:.2f}" '
                    f'fill="{BESS_FILL}" stroke="{BESS_STROKE}" '
                    f'stroke-width="{line_w * 2:.2f}"/>'
                )
            if kind in ("pv", "both"):
                out.append(
                    f'<rect x="{cx - pv_s / 2:.2f}" y="{cy - pv_s / 2:.2f}" '
                    f'width="{pv_s:.2f}" height="{pv_s:.2f}" fill="{PV_FILL}" '
                    f'stroke="{PV_STROKE}" stroke-width="{line_w * 1.6:.2f}"/>'
                )
        out.append(
            f'<text x="{cx + slot + fs * 0.55:.2f}" y="{cy + fs * 0.35:.2f}" '
            f'fill="#000000">{label}</text>'
        )
    out.append("</g>")
    return "\n".join(out)


def render(system, arc_deg, rotate_deg, size_mm, pv_size, bess_r, line_w, out_path,
           depth_scale="linear", legend="top-left", edge_style="dendrogram",
           layout="radial", force_iters=300, legend_visible=False, der_scale=1.0,
           tint=None):
    raw = REPO / "rawData" / system
    edges = parse_edges(raw / "BranchData.dss")
    if not edges:
        sys.exit(f"no Line records found in {raw / 'BranchData.dss'}")
    pv_buses = parse_buses(raw / "PVSystem.dss", "PVsystem")
    bess_buses = parse_buses(raw / "Storage.dss", "Storage")

    root = "1"
    order, depth, parent, children = build_tree(edges, root)
    pos, ang, rad = radial_layout(order, depth, children, root, arc_deg, rotate_deg,
                                  dr=1.0, depth_scale=depth_scale)
    if layout == "force":
        pos = force_layout(order, edges, pos, iterations=force_iters)
        edge_style = "straight"  # arcs are meaningless once positions are free

    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)
    w, h = maxx - minx, maxy - miny

    # Normalise so the longer axis spans 1000 viewBox units; marker and stroke
    # sizes below are expressed in those units. The viewBox tracks the true
    # bounding box aspect so a lopsided tree does not pad itself with whitespace.
    scale = 1000.0 / max(w, h)
    pad = 0.04 * 1000.0

    def sx(x):
        return (x - minx) * scale + pad

    def sy(y):
        return (y - miny) * scale + pad

    vb_w = w * scale + 2 * pad
    vb_h = h * scale + 2 * pad
    out_w_mm = size_mm
    out_h_mm = size_mm * vb_h / vb_w

    parts = [
        f'<svg {SVG_NS} version="1.1" '
        f'width="{out_w_mm:.2f}mm" height="{out_h_mm:.2f}mm" '
        f'viewBox="0 0 {vb_w:.2f} {vb_h:.2f}">',
        namedview(),
        f"<title>{system} radial network diagram</title>",
        f"<desc>{len(order)} buses, {len(edges)} branches, "
        f"{len(pv_buses)} PV, {len(bess_buses)} BESS. "
        f"Layout computed from topology (no bus coordinates in source data).</desc>",
    ]

    # Cosmetics follow ieee123C so all three system figures read as one family.
    # Ratios measured from that file against its bus-dot radius of 3.142:
    #   branch stroke 1.50 -> 0.48   (black, opacity 0.83)
    #   PV square     10.12 -> 3.22  (concentric, behind the dot)
    #   BESS ring      6.29 -> 2.00  (concentric, behind the dot)
    # The dot itself is sized from the layout's own median edge length, since
    # these feeders differ hugely in density (2.5k vs 10.3k buses).
    elens = sorted(
        math.hypot(sx(pos[u][0]) - sx(pos[parent[u]][0]),
                   sy(pos[u][1]) - sy(pos[parent[u]][1]))
        for u in order if parent[u] is not None
    )
    med_e = elens[len(elens) // 2] if elens else 10.0
    bus_r = max(0.55, min(5.0, 0.30 * med_e))

    # --- branches -----------------------------------------------------------
    ink = tint or "#000000"
    parts.append(
        f'<g id="branches" fill="none" stroke="{ink}" stroke-opacity="0.83" '
        f'stroke-width="{bus_r * 0.48:.2f}" stroke-linecap="round" '
        f'stroke-linejoin="round">'
    )
    seg = []
    for u in order:
        p = parent[u]
        if p is None:
            continue
        if edge_style == "straight":
            x1, y1 = pos[p]
            x2, y2 = pos[u]
            seg.append(f"M{sx(x1):.2f},{sy(y1):.2f}L{sx(x2):.2f},{sy(y2):.2f}")
            continue

        # Dendrogram routing: sweep along the parent's radius to the child's
        # bearing, then run radially outward. Guarantees no edge crossings,
        # which matters on deep feeders where parent and child can sit at very
        # different bearings and a straight chord would cut across the plot.
        rp, tp, tc = rad[p], ang[p], ang[u]
        ax, ay = rp * math.cos(tc), rp * math.sin(tc)
        x1, y1 = pos[p]
        x2, y2 = pos[u]
        d = f"M{sx(x1):.2f},{sy(y1):.2f}"
        if abs(tc - tp) > 1e-9 and rp > 1e-9:
            rr = rp * scale
            large = 1 if abs(tc - tp) > math.pi else 0
            sweep = 1 if tc > tp else 0
            d += f"A{rr:.2f},{rr:.2f} 0 {large} {sweep} {sx(ax):.2f},{sy(ay):.2f}"
        d += f"L{sx(x2):.2f},{sy(y2):.2f}"
        seg.append(d)
    parts.append(f'<path d="{"".join(seg)}"/>')
    parts.append("</g>")

    # A bus can host both. In ieee2552C every PV bus also hosts a BESS, so a
    # plain overlay would bury the teal circles entirely; co-located buses get a
    # circle wide enough to stay visible as a ring behind the gold square.
    both = pv_buses & bess_buses
    bess_only = bess_buses - both
    pv_only = pv_buses - both
    # ieee123's PV/BESS ratios hold at 128 buses. At 2.5k-10k buses the dot is
    # tiny and the DER markers inherit that, becoming invisible at column width,
    # so they scale independently. A deliberate departure from ieee123 -- the
    # alternative is markers no reader can resolve.
    pv_size = bus_r * 3.22 * der_scale
    bess_r = bus_r * 2.00 * der_scale
    combo_r = bess_r

    # PV first (largest), then BESS, then every bus as a black dot on top -- the
    # ieee123C stacking, so no marker buries the bus it belongs to.
    parts.append(
        f'<g id="pv-indicators" fill="{PV_FILL}" fill-opacity="0.802" '
        f'stroke="{PV_STROKE}" stroke-width="{bus_r * 0.48:.2f}">'
    )
    for b in sorted(pv_buses):
        if b not in pos:
            continue
        x, y = pos[b]
        parts.append(
            f'<rect x="{sx(x) - pv_size / 2:.2f}" y="{sy(y) - pv_size / 2:.2f}" '
            f'width="{pv_size:.2f}" height="{pv_size:.2f}"/>'
        )
    parts.append("</g>")

    parts.append(
        f'<g id="bess-indicators" fill="{BESS_FILL}" fill-opacity="0.55" '
        f'stroke="{BESS_STROKE}" stroke-width="{bus_r * 0.64:.2f}">'
    )
    for b in sorted(bess_only):
        if b in pos:
            x, y = pos[b]
            parts.append(f'<circle cx="{sx(x):.2f}" cy="{sy(y):.2f}" r="{bess_r:.2f}"/>')
    for b in sorted(both):
        if b in pos:
            x, y = pos[b]
            parts.append(f'<circle cx="{sx(x):.2f}" cy="{sy(y):.2f}" r="{combo_r:.2f}"/>')
    parts.append("</g>")

    parts.append(
        f'<g id="buses" fill="{ink}" fill-opacity="0.83" stroke="{ink}" '
        f'stroke-width="{bus_r * 0.1:.2f}">'
    )
    for b in order:
        x, y = pos[b]
        parts.append(f'<circle cx="{sx(x):.2f}" cy="{sy(y):.2f}" r="{bus_r:.2f}"/>')
    parts.append("</g>")

    # --- substation ---------------------------------------------------------
    rx, ry = pos[root]
    parts.append(
        f'<g id="substation"><circle cx="{sx(rx):.2f}" cy="{sy(ry):.2f}" '
        f'r="{bus_r * 2.6:.2f}" fill="{SUBSTATION_FILL}" stroke="{SUBSTATION_STROKE}" '
        f'stroke-width="{bus_r * 0.5:.2f}"/></g>'
    )

    if legend:
        parts.append(
            legend_svg(vb_w, vb_h, legend, pv_size, bess_r, bus_r * 0.48,
                       len(pv_only), len(bess_only), len(both), combo_r,
                       hidden=not legend_visible)
        )

    parts.append("</svg>")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(parts))
    print(
        f"  {system}: {len(order)} buses, {len(edges)} branches, "
        f"{len(pv_buses)} PV, {len(bess_buses)} BESS, max depth {max(depth.values())}"
    )
    try:
        shown = out_path.relative_to(REPO)
    except ValueError:
        shown = out_path
    print(f"  wrote {shown} ({out_path.stat().st_size / 1024:.0f} kB)")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("system", help="rawData subfolder, e.g. large10kC_1ph")
    ap.add_argument("--arc", type=float, default=360.0, help="angular span in degrees")
    ap.add_argument("--rotate", type=float, default=-90.0, help="rotation in degrees")
    ap.add_argument("--size-mm", type=float, default=88.9, help="IEEE column width")
    ap.add_argument("--pv-size", type=float, default=4.0, help="PV square side, viewBox units")
    ap.add_argument("--bess-r", type=float, default=2.4, help="BESS circle radius")
    ap.add_argument("--line-w", type=float, default=0.7, help="branch stroke width")
    ap.add_argument("--depth-scale", choices=("linear", "sqrt"), default="linear",
                    help="hop-count to radius mapping")
    ap.add_argument("--layout", choices=("radial", "force"), default="radial",
                    help="radial suits hub-and-spoke feeders; force suits deep sprawling ones")
    ap.add_argument("--force-iters", type=int, default=300, help="force layout iterations")
    ap.add_argument("--edge-style", choices=("dendrogram", "straight"), default="dendrogram",
                    help="dendrogram routing avoids edge crossings on deep feeders")
    ap.add_argument("--legend", default="top-left",
                    choices=("top-left", "top-right", "bottom-left", "bottom-right", "none"),
                    help="legend corner, or 'none' to leave it out of the file entirely")
    ap.add_argument("--tint", nargs="?", const="auto", default=None, metavar="HEX",
                    help="colour the buses and branches with the system's ink "
                         "(paper colour mixed into black); bare flag derives it "
                         "from the system name, or pass a hex")
    ap.add_argument("--der-scale", type=float, default=1.0,
                    help="enlarge PV/BESS markers relative to the ieee123 ratio; "
                         "needed on dense feeders where the bus dot is tiny")
    ap.add_argument("--legend-visible", action=argparse.BooleanOptionalAction, default=False,
                    help="render the legend; when off it is still written to the SVG "
                         "but display:none, so it stays togglable in Inkscape")
    ap.add_argument("--out", default=None, help="output path (default assets/networks/<system>.svg)")
    args = ap.parse_args()

    out = Path(args.out) if args.out else REPO / "assets" / "networks" / f"{args.system}.svg"
    render(
        args.system, args.arc, args.rotate, args.size_mm,
        args.pv_size, args.bess_r, args.line_w, out,
        depth_scale=args.depth_scale,
        legend=None if args.legend == "none" else args.legend,
        edge_style=args.edge_style,
        layout=args.layout,
        force_iters=args.force_iters,
        legend_visible=args.legend_visible,
        der_scale=args.der_scale,
        tint=(ink_for(args.system) if args.tint == "auto" else args.tint),
    )


if __name__ == "__main__":
    main()
