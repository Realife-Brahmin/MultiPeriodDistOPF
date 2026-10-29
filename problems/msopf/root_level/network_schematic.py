# Network schematics of the MSOPF systems as layered, Inkscape-editable SVG (2026-10-29), in the palette of the
# SPAR-DDP paper's network figures (TPWRS-2027-Sparse-Parallel-DDP-for-MPOPF/figures/networks/README.md):
# slate network, amber PV squares, lavender battery rings, and the substations as concentric rings in the
# PESGM2027 paper's substation colours (figures/colors.tex cS1-cS4; an added substation drawn dashed).
#   python problems/msopf/root_level/network_schematic.py <ieee123|med3980|sf6937> <out.svg>
#   inkscape <out.svg> --export-area-drawing --export-type=pdf --export-filename=<out.pdf>
# Topology, PV, batteries, ties and substations come from the deck in test_networks/. Bus positions:
#   ieee123  the hand layout of the SPAR-DDP ieee123 figure: its bus dots matched to their number labels
#            (bus,x,y,lx,ly in network_layout_ieee123.csv beside this script, made with inkscape --query-all);
#   med3980  test_networks/ieee9500_source/base/BusCoords.dss;
#   sf6937   every Buscoords.dss under test_networks/sources/smartds_sfo_p1u (longitude, latitude).
# Buses without coordinates (service-transformer secondaries, internal switch nodes, substation buses) take
# the position of a neighbour, so a customer bus sits on its primary bus.
import re, sys, math, csv
from collections import defaultdict, deque
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
NET = REPO / "test_networks"
SLATE, SLATE2 = "#51435F", "#61566E"
INK, INK2 = "#333333", "#555555"                      # neutral lines for the large feeders (no DDP-paper slate)
PV_FILL, PV_LINE, BAT_FILL, BAT_LINE = "#C99746", "#8C612F", "#EEE7F6", "#7956A6"
SUBC = ["#1F3B73", "#E3A21A", "#169C78", "#874FAF"]          # cS1-cS4 of the PESGM2027 figures/colors.tex
# Buses inside a substation (transformer and regulator nodes), not named as its neighbours.
INTERNAL = re.compile(r"hvmv|regxfmr|_lsb|_hsb|^sb\d+_p1u")
TIES_9500 = ["a8645_48332_sw", "tsw803273_sw", "v7173_48332_sw", "a333_48332_sw", "tsw568613_sw"]

CONFIG = {
    "ieee123": dict(deck="ieee123_5poi_1ph", variant=None, closed=[], subs=["subs1", "subs2", "subs3", "subs4"], added={},
                    width=1500, lw=3.2, pv=19, bat=13, sub=34, dot=6.5, labels=True),
    "med3980": dict(deck="ieee9500xf_1ph", variant=None, closed=TIES_9500, subs=["subs1", "subs2", "subs3", "subs4"],
                    added={"subs4": "l2804256"}, width=3000, lw=3.0, pv=11, bat=6.5, sub=24, dot=0, labels=False, lens=230, tf=110, nf=88, mk=3.0),
    "sf6937": dict(deck="smartds_p1u_xf_1ph", variant="Master_4sub_typeIII.dss", closed=[], subs=["subs1", "subs2", "subs3", "subs4"],
                   added={}, width=3000, lw=2.6, pv=10, bat=6.0, sub=24, dot=0, labels=False, lens=230, tf=110, nf=88, mk=3.0,
                   names={"subs1": "hs2", "subs2": "hs1", "subs3": "hs3", "subs4": "hs0"}),
}

bus = lambda s: s.split(".")[0].lower()

def deck(cfg):
    "Branches (name, a, b, tie?), PV and battery buses, substation buses."
    d = NET / cfg["deck"]; lines = []; enabled = set(t.lower() for t in cfg["closed"])
    if cfg["variant"]:
        for l in (d / cfg["variant"]).read_text().splitlines():
            m = re.match(r"(?i)line\.(\S+)\.enabled=true", l.strip())
            if m: enabled.add(m.group(1).lower())
    for f, tie in (("Master.dss", False), ("BranchData.dss", False), ("Ties.dss", True)):   # Master: ieee123's POI lines
        if not (d / f).exists():
            continue
        for l in (d / f).read_text().splitlines():
            m = re.match(r"(?i)new line\.(\S+).*?bus1=(\S+).*?bus2=(\S+)", l)
            if not m: continue
            name = m.group(1).lower()
            off = re.search(r"(?i)enabled=false", l) is not None
            if off and name not in enabled: continue
            lines.append((name, bus(m.group(2)), bus(m.group(3)), tie or off))
    grab = lambda f: sorted({bus(b) for b in re.findall(r"(?i)bus1=(\S+)", (d / f).read_text())})
    src = {m.group(1).lower(): bus(m.group(2)) for m in re.finditer(r"(?im)^new vsource\.(\S+).*?bus1=(\S+)", (d / "Master.dss").read_text())}
    src.update(cfg["added"])
    return lines, grab("PVSystem.dss"), grab("Storage.dss"), {s: src[s] for s in cfg["subs"]}

def coords(system):
    if system == "ieee123":
        with open(Path(__file__).with_name("network_layout_ieee123.csv")) as f:
            r = list(csv.DictReader(f))
        return {x["bus"]: (float(x["x"]), float(x["y"])) for x in r}, {x["bus"]: (float(x["lx"]), float(x["ly"])) for x in r}, False
    if system == "med3980":
        xy = {}
        for l in (NET / "ieee9500_source/base/BusCoords.dss").read_text().splitlines():
            p = [q.strip() for q in l.split(",")]
            if len(p) == 3: xy[p[0].lower()] = (float(p[1]), float(p[2]))
        return xy, {}, True
    xy = {}
    for f in (NET / "sources/smartds_sfo_p1u").rglob("*"):
        if f.name.lower() != "buscoords.dss": continue
        for l in f.read_text().splitlines():
            p = l.split()
            if len(p) == 3: xy[p[0].lower()] = (float(p[1]), float(p[2]))
    lat0 = sum(v[1] for v in xy.values()) / len(xy)
    # rotated 90 degrees clockwise (north to the right), so the tall sub-region lies across the page
    return {k: (v[1], -v[0] * math.cos(math.radians(lat0))) for k, v in xy.items()}, {}, True

def place(lines, xy):
    "Give every bus without coordinates its nearest placed neighbour's position (breadth first)."
    adj = defaultdict(list)
    for _, a, b, _ in lines: adj[a].append(b); adj[b].append(a)
    pos = {b: xy[b] for b in adj if b in xy}; q = deque(pos)
    while q:
        u = q.popleft()
        for v in adj[u]:
            if v not in pos: pos[v] = pos[u]; q.append(v)
    return pos, [b for b in adj if b not in pos]

def short(system, b):
    "A bus name short enough to print in a lens."
    if system == "sf6937":
        m = re.match(r"sb\d+_p1u(hs\d)_1247_(node|bus)_(\S+)", b)
        if m: return f"{m.group(1)} {m.group(2)} {m.group(3)}"
        m = re.match(r"p1u(dt\d+)", b)
        if m: return m.group(1)
        return b.replace("p1u", "").replace("_1247", "")
    return b.replace("_hsb", "")

def body(lines, pos, pvs, bats, T, lw, pv, bat, near=None, fade=1.0):
    "Lines, ties, batteries and PV under the map T; with near, only what lies within near = (centre, radius)."
    inside = (lambda q: math.dist(q, near[0]) <= near[1]) if near else (lambda q: True)
    # faded: lighter solid colours (per-element opacity would make every element a PDF transparency group)
    sl, sl2, bf, bl, pf, pl = ((INK, INK2, BAT_FILL, BAT_LINE, PV_FILL, PV_LINE) if fade >= 1 else
                               ("#9C9C9C", "#9C9C9C", "#F1ECF7", "#A892C4", "#D2AF72", "#A98A5C"))
    o = []
    for name, a, b, tie in lines:
        if not (inside(pos[a]) or inside(pos[b])): continue
        (ax, ay), (bx, by) = T(pos[a]), T(pos[b])
        if abs(ax - bx) + abs(ay - by) < 1e-6: continue
        st = (f'stroke="{sl2}" stroke-width="{2.2 * lw:.2f}" stroke-dasharray="{6 * lw:.1f},{4 * lw:.1f}"' if tie
              else f'stroke="{sl}" stroke-width="{lw:.2f}"')
        o.append(f'<line x1="{ax:.1f}" y1="{ay:.1f}" x2="{bx:.1f}" y2="{by:.1f}" {st} stroke-linecap="round"/>')
    for b in bats:
        if not inside(pos[b]): continue
        x, y = T(pos[b])
        o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{bat * 1.35:.1f}" fill="{bf}" stroke="{bl}" stroke-width="{bat * 0.28:.2f}"/>')
    for b in pvs:
        if not inside(pos[b]): continue
        x, y = T(pos[b])
        o.append(f'<rect x="{x - pv / 2:.1f}" y="{y - pv / 2:.1f}" width="{pv}" height="{pv}" fill="{pf}" stroke="{pl}" stroke-width="{pv * 0.1:.2f}"/>')
    return o

def ring(x, y, R, colr, dashed):
    dash = ' stroke-dasharray="6,4"' if dashed else ""
    return [f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{R}" fill="white" fill-opacity="0.75" stroke="{colr}" stroke-width="{R * 0.16:.1f}"{dash}/>',
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{R * 0.62:.1f}" fill="none" stroke="{colr}" stroke-width="{R * 0.12:.1f}"{dash}/>']

def lens_main(system, out, cfg, lines, pvs, bats, subs, xy, pos):
    """Large networks: the whole feeder drawn light, and each substation in a lens outside it that magnifies its
    neighbourhood, with the substation's bus and a few nearby buses named; a leader joins the lens to the spot."""
    adj = defaultdict(list)
    for _, a, b, _ in lines: adj[a].append(b); adj[b].append(a)
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys); k = cfg["width"] / (x1 - x0)
    T = lambda p: ((p[0] - x0) * k, (y1 - p[1]) * k)
    W, H = (x1 - x0) * k, (y1 - y0) * k
    # Small lenses that only magnify the geometry; the names go in a short caption under each lens, so the
    # lenses need not be large enough to hold text and the feeder keeps most of the panel.
    RL, gap = cfg["lens"], cfg["lens"] * 0.3
    tf, nf, F = cfg["tf"], cfg["nf"], cfg["mk"]                       # title and name font sizes, mark scale
    above, below = 1.25 * tf, 2.5 * nf                                # text heights over and under a lens
    placed = []
    for i, (s, b) in enumerate(subs.items()):
        px, py = T(pos[b])
        side = min((px, "left"), (W - px, "right"), (py, "top"), (H - py, "bottom"))[1] if cfg.get("sides") is None else cfg["sides"][i]
        c = {"left": (-gap - RL, py), "right": (W + gap + RL, py), "top": (px, -gap - RL - below),
             "bottom": (px, H + gap + RL + above)}[side]
        placed.append([i, s, b, side, list(c)])
    for side in ("left", "right", "top", "bottom"):                 # lenses on one side must not overlap
        L = sorted((p for p in placed if p[3] == side), key=lambda p: p[4][1] if side in ("left", "right") else p[4][0])
        ax = 1 if side in ("left", "right") else 0
        need = 2 * RL + above + below + 0.5 * gap if ax == 1 else max(2 * RL, 12 * nf) + gap
        for a_, b_ in zip(L, L[1:]):
            if b_[4][ax] - a_[4][ax] < need: b_[4][ax] = a_[4][ax] + need
    o = ['<g inkscape:groupmode="layer" inkscape:label="network" id="layer-network">']
    o += body(lines, pos, pvs, bats, T, cfg["lw"], cfg["pv"], cfg["bat"], fade=0.55)   # light, but solid
    o.append("</g>")
    o.append('<g inkscape:groupmode="layer" inkscape:label="lenses" id="layer-lenses"><defs>')
    for i, s, b, side, (cx, cy) in placed:
        o.append(f'<clipPath id="lens{i + 1}"><circle cx="{cx:.1f}" cy="{cy:.1f}" r="{RL * 0.97:.1f}"/></clipPath>')
    o.append("</defs>")
    allx, ally = [0, W], [0, H]
    txt = lambda x, y, size, colr, t, bold=False: (f'<text x="{x:.1f}" y="{y:.1f}" font-family="Times New Roman" font-size="{size:.0f}" '
                                                    f'{"font-weight=" + chr(34) + "bold" + chr(34) if bold else ""} fill="{colr}" text-anchor="middle">{t}</text>')
    for i, s, b, side, (cx, cy) in placed:
        colr, added = SUBC[i], s in cfg["added"]
        p = pos[b]; px, py = T(p)
        # the neighbourhood: radius set by the 10th nearest distinct placed bus
        d = sorted({round(math.dist(p, xy[q]), 9) for q in xy if q in pos and math.dist(p, xy[q]) > 0})
        rn = 1.15 * d[min(9, len(d) - 1)]; z = 0.9 * RL / (rn * k)
        Tl = lambda q, z=z, px=px, py=py, cx=cx, cy=cy: (cx + (T(q)[0] - px) * z, cy + (T(q)[1] - py) * z)
        dx, dy = cx - px, cy - py; L = math.hypot(dx, dy)
        o.append(f'<line x1="{px:.1f}" y1="{py:.1f}" x2="{cx - dx / L * RL:.1f}" y2="{cy - dy / L * RL:.1f}" stroke="{colr}" stroke-width="{1.6 * F:.1f}"/>')
        o.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{9 * F:.1f}" fill="white" stroke="{colr}" stroke-width="{3 * F:.1f}"/>')
        o.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{RL:.1f}" fill="white" stroke="{colr}" stroke-width="{2 * F:.1f}"/>')
        o.append(f'<g clip-path="url(#lens{i + 1})">')
        o += body(lines, pos, pvs, bats, Tl, 2.2 * F, 12 * F, 7 * F, near=(p, 1.6 * rn))
        for q in xy:
            if q in pos and math.dist(pos[q], p) <= 1.6 * rn:
                x, y = Tl(pos[q]); o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{4 * F:.1f}" fill="{INK}"/>')
        o += ring(cx, cy, cfg["sub"] * F, colr, added)
        o.append("</g>")
        # the substation's bus and its two nearest named buses (breadth first), under the lens
        seen, order, q = {b}, [], deque([b])
        while q and len(order) < 2:
            u = q.popleft()
            for v in adj[u]:
                if v in seen: continue
                seen.add(v); q.append(v)
                if v in xy and not INTERNAL.search(v) and short(system, v) not in [short(system, w) for w in order]:
                    order.append(v)                                   # the first feeder buses, not station nodes
        title = f"S{i + 1}" + (f" ({cfg['names'][s]})" if cfg.get("names") else "") + (" (added)" if added else "")
        o.append(txt(cx, cy - RL - 0.35 * tf, tf, colr, title, True))
        o.append(txt(cx, cy + RL + 1.05 * nf, nf, colr, short(system, b), True))
        o.append(txt(cx, cy + RL + 2.15 * nf, 0.85 * nf, INK, "near " + ", ".join(short(system, v) for v in order[:2])))
        allx += [cx - max(RL, 4.5 * nf), cx + max(RL, 4.5 * nf)]; ally += [cy - RL - above, cy + RL + below]
    o.append("</g>")
    m = 20; vx, vy = min(allx) - m, min(ally) - m; vw, vh = max(allx) - vx + m, max(ally) - vy + m
    head = [f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
            f'width="{vw:.0f}" height="{vh:.0f}" viewBox="{vx:.1f} {vy:.1f} {vw:.1f} {vh:.1f}">',
            f'<!-- {system}: {cfg["deck"]}; written by problems/msopf/root_level/network_schematic.py -->']
    Path(out).write_text("\n".join(head + o + ["</svg>"]), encoding="utf-8")
    print(f"{system}: lenses on {[p[3] for p in placed]} -> {out}")

def main(system, out):
    cfg = CONFIG[system]
    if cfg.get("lens"):
        lines, pvs, bats, subs = deck(cfg)
        xy, _, _ = coords(system)
        pos, _ = place(lines, xy)
        return lens_main(system, out, cfg, lines, pvs, bats, subs, xy, pos)
    lines, pvs, bats, subs = deck(cfg)
    xy, lab, flip = coords(system)
    pos, lost = place(lines, xy)
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys); k = cfg["width"] / (x1 - x0); pad = 2 * cfg["sub"]
    T = lambda p: (pad + (p[0] - x0) * k, pad + ((y1 - p[1]) if flip else (p[1] - y0)) * k)
    if system == "ieee123": k = 1.0; T = lambda p: (pad + p[0] - x0, pad + p[1] - y0)
    W, H = T((x1, y0 if flip else y1))[0] + pad, T((x0, y0 if flip else y1))[1] + pad
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
         f'width="{W:.0f}" height="{H:.0f}" viewBox="0 0 {W:.1f} {H:.1f}">',
         f'<!-- {system}: {cfg["deck"]}; written by problems/msopf/root_level/network_schematic.py -->']
    layer = lambda name: o.append(f'<g inkscape:groupmode="layer" inkscape:label="{name}" id="layer-{name}">')
    layer("network")
    for name, a, b, tie in lines:
        if tie: continue
        (ax, ay), (bx, by) = T(pos[a]), T(pos[b])
        if abs(ax - bx) + abs(ay - by) < 1e-6: continue
        o.append(f'<line x1="{ax:.1f}" y1="{ay:.1f}" x2="{bx:.1f}" y2="{by:.1f}" stroke="{SLATE}" stroke-width="{cfg["lw"]}" stroke-linecap="round"/>')
    o.append("</g>"); layer("ties")
    for name, a, b, tie in lines:
        if not tie: continue
        (ax, ay), (bx, by) = T(pos[a]), T(pos[b])
        o.append(f'<line x1="{ax:.1f}" y1="{ay:.1f}" x2="{bx:.1f}" y2="{by:.1f}" stroke="{SLATE2}" stroke-width="{2.2 * cfg["lw"]}" '
                 f'stroke-dasharray="{6 * cfg["lw"]},{4 * cfg["lw"]}" stroke-linecap="round"/>')
    o.append("</g>")
    if cfg["dot"]:
        layer("buses")
        for b in xy:
            if b in pos:
                x, y = T(pos[b]); o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{cfg["dot"]}" fill="{SLATE}"/>')
        o.append("</g>")
    layer("batteries")
    for b in bats:
        x, y = T(pos[b]); r = cfg["bat"]
        o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r * 1.35:.1f}" fill="{BAT_FILL}" fill-opacity="0.85" stroke="{BAT_LINE}" stroke-width="{r * 0.28:.2f}"/>')
    o.append("</g>"); layer("pv")
    for b in pvs:
        x, y = T(pos[b]); s = cfg["pv"]
        o.append(f'<rect x="{x - s / 2:.1f}" y="{y - s / 2:.1f}" width="{s}" height="{s}" fill="{PV_FILL}" fill-opacity="0.8" stroke="{PV_LINE}" stroke-width="{s * 0.1:.2f}"/>')
    o.append("</g>")
    layer("substations")
    for i, (s, b) in enumerate(subs.items()):
        x, y = T(pos[b]); R = cfg["sub"]; dash = ' stroke-dasharray="6,4"' if s in cfg["added"] else ""
        o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{R}" fill="white" fill-opacity="0.6" stroke="{SUBC[i]}" stroke-width="{R * 0.16:.1f}"{dash}/>')
        o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{R * 0.62:.1f}" fill="none" stroke="{SUBC[i]}" stroke-width="{R * 0.12:.1f}"{dash}/>')
        o.append(f'<text x="{x + 1.15 * R:.1f}" y="{y - 1.0 * R:.1f}" font-family="Times New Roman" font-weight="bold" font-size="{max(26, R * 1.2):.0f}" fill="{SUBC[i]}">S{i + 1}</text>')
    o.append("</g>")
    if cfg["labels"]:
        layer("labels")
        for b, (lx, ly) in lab.items():
            x, y = T((lx, ly)); o.append(f'<text x="{x:.1f}" y="{y + 9:.1f}" font-family="Times New Roman" font-size="26" fill="#3B3346" text-anchor="middle">{b}</text>')
        o.append("</g>")
    o.append("</svg>")
    Path(out).write_text("\n".join(o), encoding="utf-8")
    nt = sum(t for *_, t in lines)
    print(f"{system}: {len(lines) - nt} lines, {nt} ties, {len(pvs)} PV, {len(bats)} batteries, substations {subs}; "
          f"{len(xy)} coordinates, {len(lost)} buses unplaced -> {out}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
