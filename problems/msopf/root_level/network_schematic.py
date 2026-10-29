# Network schematics of the MSOPF systems as layered, Inkscape-editable SVG (2026-10-29), in the palette of the
# SPAR-DDP paper's network figures (TPWRS-2027-Sparse-Parallel-DDP-for-MPOPF/figures/networks/README.md):
# slate network, amber PV squares, lavender battery rings, and the substations as concentric rings in the
# PESGM2027 paper's substation colours (figures/colors.tex cS1-cS4). The large feeders (med3980, sf6937) are drawn
# in grey with the ieee123 figure's markers around dark bus dots, and carry a table of their substation buses.
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
# The large feeders' markers: the ieee123 figure's colours at full strength (paler ones vanish at this scale), a
# slightly deeper lavender fill, and a light grey for the lines.
PV_SOFT, PV_LINE_SOFT, BAT_SOFT, BAT_FILL_BIG, GREY = PV_FILL, PV_LINE, BAT_LINE, "#DDD0EE", "#9A9A9A"
SUBC = ["#1F3B73", "#E3A21A", "#169C78", "#874FAF"]          # cS1-cS4 of the PESGM2027 figures/colors.tex
# Buses inside a substation (transformer and regulator nodes), not named as its neighbours.
INTERNAL = re.compile(r"hvmv|regxfmr|_lsb|_hsb|_int$|^sb\d+_p1u")
TIES_9500 = ["a8645_48332_sw", "tsw803273_sw", "v7173_48332_sw", "a333_48332_sw", "tsw568613_sw"]

CONFIG = {
    "ieee123": dict(deck="ieee123_5poi_1ph", variant=None, closed=[], subs=["subs1", "subs2", "subs3", "subs4"], added={},
                    width=1500, lw=3.2, pv=19, bat=13, sub=34, dot=6.5, labels=True),
    "med3980": dict(deck="ieee9500xf_1ph", variant=None, closed=TIES_9500, subs=["subs1", "subs2", "subs3", "subs4"],
                    added={"subs4": "l2804256"}, width=3000, lw=4.6, pv=18, bat=10.8, sub=24, dot=0, labels=False, height_mm=60, mm=dict(mark=2.0, label=3.6, table=2.5)),
    "sf6937": dict(deck="smartds_p1u_xf_1ph", variant="Master_4sub_typeIII.dss", closed=[], subs=["subs1", "subs2", "subs3", "subs4"],
                   added={}, width=3000, lw=2.6, pv=10, bat=6.0, sub=24, dot=0, labels=False, height_mm=60, mm=dict(mark=2.0, label=3.6, table=2.5),
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
    "A bus name short enough for the substation table."
    if system == "sf6937":
        m = re.match(r"sb\d+_p1u(hs\d)_1247_(node|bus)_(\S+)", b)
        if m: return f"{m.group(1)} {m.group(2)} {m.group(3)}"
        m = re.match(r"p1u(dt\d+)", b)
        if m: return m.group(1)
        return b.replace("p1u", "").replace("_1247", "")
    return b.replace("_hsb", "")

def map_main(system, out, cfg, lines, pvs, bats, subs, pos):
    """Large networks: the whole feeder, its substations marked (half-transparent white inside each ring, so the
    buses and lines beneath stay visible), and a small table of the substation buses in the emptiest corner of the drawing. Marks, table
    and text are sized in printed millimetres: the feeder is printed cfg["height_mm"] tall."""
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys); k = cfg["width"] / (x1 - x0)
    T = lambda p: ((p[0] - x0) * k, (y1 - p[1]) * k)
    W, H = (x1 - x0) * k, (y1 - y0) * k
    u = H / cfg["height_mm"]                                          # drawing units per printed mm
    S = cfg["mm"]
    R, fs, tf = S["mark"] * u, S["label"] * u, S["table"] * u
    txt = lambda x, y, size, colr, t, bold=False, anchor="start", halo=True: (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Times New Roman" font-size="{size:.1f}"'
        f'{" font-weight=" + chr(34) + "bold" + chr(34) if bold else ""} fill="{colr}" text-anchor="{anchor}"'
        + (f' stroke="white" stroke-opacity="0.6" stroke-width="{0.12 * size:.1f}" paint-order="stroke"' if halo else "") + f'>{t}</text>')
    icon = lambda x, y, r, colr, fop=1.0: (
        # a wider white stroke under the ring first, so it stays apart from markers of a similar colour
        [f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" fill="none" stroke="white" stroke-width="{0.48 * r:.1f}"/>',
         f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" fill="white" fill-opacity="{fop}" stroke="{colr}" stroke-width="{0.2 * r:.1f}"/>',
         f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{0.6 * r:.1f}" fill="none" stroke="{colr}" stroke-width="{0.16 * r:.1f}"/>',
         f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{0.2 * r:.1f}" fill="{colr}"/>'])

    o = ['<g inkscape:groupmode="layer" inkscape:label="network" id="layer-network">']
    for name, a, b, tie in lines:
        (ax, ay), (bx, by) = T(pos[a]), T(pos[b])
        if abs(ax - bx) + abs(ay - by) < 1e-6: continue
        lw = cfg["lw"]
        st = (f'stroke="{GREY}" stroke-width="{2.2 * lw:.2f}" stroke-dasharray="{6 * lw:.1f},{4 * lw:.1f}"' if tie
              else f'stroke="{GREY}" stroke-width="{lw:.2f}"')
        o.append(f'<line x1="{ax:.1f}" y1="{ay:.1f}" x2="{bx:.1f}" y2="{by:.1f}" {st} stroke-linecap="round"/>')
    o.append("</g>")
    # PV and batteries as in the ieee123 figure: a pale marker around a dark bus dot
    o.append('<g inkscape:groupmode="layer" inkscape:label="der" id="layer-der">')
    rb, sp = cfg["bat"] * 1.35, cfg["pv"]
    for b in bats:
        x, y = T(pos[b])
        o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{rb:.1f}" fill="{BAT_FILL_BIG}" stroke="{BAT_SOFT}" stroke-width="{0.2 * rb:.2f}"/>')
    for b in pvs:
        x, y = T(pos[b])
        o.append(f'<rect x="{x - sp / 2:.1f}" y="{y - sp / 2:.1f}" width="{sp:.1f}" height="{sp:.1f}" fill="{PV_SOFT}" '
                 f'stroke="{PV_LINE_SOFT}" stroke-width="{0.1 * sp:.2f}"/>')
    for b in sorted(set(pvs) | set(bats)):
        x, y = T(pos[b])
        o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{0.3 * min(rb, sp / 2):.1f}" fill="{INK}"/>')
    o.append("</g>")

    # The substations: the feeder shows through each ring; the name beside it, outward from the centre.
    o.append('<g inkscape:groupmode="layer" inkscape:label="substations" id="layer-substations">')
    cx_, cy_ = W / 2, H / 2
    for i, (s, b) in enumerate(subs.items()):
        x, y = T(pos[b]); colr = SUBC[i]
        o += icon(x, y, R, colr, fop=0.45)
        right = x >= cx_ or x < 0.15 * W
        o.append(txt(x + (1.5 if right else -1.5) * R, y - 0.9 * R, fs, colr, f"S{i + 1}", True, "start" if right else "end"))

    # The table of substation buses, in the emptiest corner: points sampled along every line and at every marker
    # must stay out of its box (with a margin); if no corner is clear it goes beside the drawing, bottom right.
    rows = [(f"S{i + 1}", short(system, b)) for i, (s, b) in enumerate(subs.items())]
    # Times New Roman advance widths (em) of the characters used, for the column widths
    adv = lambda t, bold=False: sum(0.25 if c == " " else 0.78 if c in "mM" else 0.56 if c.isupper() else
                                    0.5 if c.isdigit() or c == "_" else 0.45 for c in t) * (1.05 if bold else 1.0)
    c1 = 0.25 * tf + 1.2 * R + 0.35 * tf + adv("S4", True) * tf + 0.8 * tf    # icon, name, gap
    tw = c1 + max(adv(r[1]) for r in rows) * tf + 0.3 * tf
    rh = 1.3 * tf; th = rh * (len(rows) + 1) + 0.55 * tf
    pts = []
    for _, a, b, _ in lines:
        (ax, ay), (bx, by) = T(pos[a]), T(pos[b]); n = max(1, int(math.hypot(bx - ax, by - ay) / (0.4 * u)))
        pts += [(ax + (bx - ax) * j / n, ay + (by - ay) * j / n) for j in range(n + 1)]
    pts += [T(pos[b]) for b in subs.values()]
    pad, over = 0.8 * u, 4.0 * u                                      # clearance; overhang allowed past the right edge
    clear = lambda X, Y: not any(X - pad <= px <= X + tw + pad and Y - pad <= py <= Y + th + pad for px, py in pts)
    best = None
    for X in [W + over - tw - j * 0.5 * u for j in range(int((W + over - tw) / (0.5 * u)) + 1)]:   # from the right
        for Y in [H - th - j * 0.5 * u for j in range(int((H - th) / (0.5 * u)) + 1)]:            # from the bottom
            if clear(X, Y):
                score = (W + over - tw - X) + (H - th - Y)            # nearest the bottom-right corner
                if best is None or score < best[0]: best = (score, X, Y)
    TX, TY = (best[1], best[2]) if best else (W + 2 * u, H - th)
    print(f"  table {tw / u:.1f} x {th / u:.1f} mm at ({TX / u:.1f}, {TY / u:.1f}) mm"
          + ("" if best else " (no clear corner: beside the drawing)"))
    rule = lambda y, w: f'<line x1="{TX:.1f}" y1="{y:.1f}" x2="{TX + tw:.1f}" y2="{y:.1f}" stroke="{INK}" stroke-width="{w * u:.2f}"/>'
    o.append('<g id="substation-table">')
    o.append(rule(TY, 0.18))
    y = TY + 1.0 * tf
    o.append(txt(TX + 0.25 * tf, y, tf, INK, "Sub.", False, "start", False))
    o.append(txt(TX + c1, y, tf, INK, "Bus", False, "start", False))
    y += 0.35 * tf; o.append(rule(y, 0.1))
    for i, (lab, bb) in enumerate(rows):
        y += rh; colr = SUBC[i]
        o += icon(TX + 0.25 * tf + 0.6 * R, y - 0.33 * tf, 0.6 * R, colr)
        o.append(txt(TX + 0.25 * tf + 1.2 * R + 0.35 * tf, y, tf, colr, lab, True, "start", False))
        o.append(txt(TX + c1, y, tf, INK, bb, False, "start", False))
    o.append(rule(y + 0.45 * tf, 0.18))
    o.append("</g></g>")
    right = max(W, TX + tw)
    m = 0.5 * u; vx, vy, vw, vh = -m, -m, right + 2 * m, max(H, TY + th) + 2 * m
    print(f"  {u:.1f} drawing units per printed mm; panel {vw / u:.0f} x {vh / u:.0f} mm")
    head = [f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
            f'width="{vw:.0f}" height="{vh:.0f}" viewBox="{vx:.1f} {vy:.1f} {vw:.1f} {vh:.1f}">',
            f'<!-- {system}: {cfg["deck"]}; written by problems/msopf/root_level/network_schematic.py -->']
    Path(out).write_text("\n".join(head + o + ["</svg>"]), encoding="utf-8")
    print(f"{system}: -> {out}")

def main(system, out):
    cfg = CONFIG[system]
    if cfg.get("mm"):                                                 # the large feeders: marks and a table
        lines, pvs, bats, subs = deck(cfg)
        xy, _, _ = coords(system)
        pos, _ = place(lines, xy)
        return map_main(system, out, cfg, lines, pvs, bats, subs, pos)
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
