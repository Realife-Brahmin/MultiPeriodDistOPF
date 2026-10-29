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
PV_FILL, PV_LINE, BAT_FILL, BAT_LINE = "#C99746", "#8C612F", "#EEE7F6", "#7956A6"
SUBC = ["#1F3B73", "#E3A21A", "#169C78", "#874FAF"]          # cS1-cS4 of the PESGM2027 figures/colors.tex
TIES_9500 = ["a8645_48332_sw", "tsw803273_sw", "v7173_48332_sw", "a333_48332_sw", "tsw568613_sw"]

CONFIG = {
    "ieee123": dict(deck="ieee123_5poi_1ph", variant=None, closed=[], subs=["subs1", "subs2", "subs3", "subs4"], added={},
                    width=1500, lw=3.2, pv=19, bat=13, sub=34, dot=6.5, labels=True),
    "med3980": dict(deck="ieee9500xf_1ph", variant=None, closed=TIES_9500, subs=["subs1", "subs2", "subs3", "subs4"],
                    added={"subs4": "l2804256"}, width=1500, lw=1.3, pv=7, bat=4.6, sub=26, dot=0, labels=False),
    "sf6937": dict(deck="smartds_p1u_xf_1ph", variant="Master_4sub_typeIII.dss", closed=[], subs=["subs1", "subs2", "subs3", "subs4"],
                   added={}, width=1500, lw=1.1, pv=6, bat=4.0, sub=26, dot=0, labels=False),
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
    return {k: (v[0] * math.cos(math.radians(lat0)), v[1]) for k, v in xy.items()}, {}, True

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

def main(system, out):
    cfg = CONFIG[system]
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
