# The MSOPF ieee123 schematic from the SPAR-DDP Inkscape figure (2026-10-29): same layout, fonts, line widths
# and bus dots, with this deck's PV and batteries (test_networks/ieee123_5poi_1ph) and its four substations.
#   python problems/msopf/root_level/network_ieee123_template.py <template.svg> <query.csv> <out.svg>
# template: TPWRS-2027-Sparse-Parallel-DDP-for-MPOPF/figures/networks/ieee123C_1ph.svg
# query:    inkscape --query-all <template.svg> (object bounding boxes in px), used to map bus positions
#           (network_layout_ieee123.csv) into each marker group's own coordinates.
# Edits: the children of "pv-indicator-squares" and "battery-bus-circles" are redrawn at this deck's buses in the
# template's own marker styles; "Substation1_if_alone" is hidden and the four substations (buses 1, 33, 120, 116) get the
# icon of the large feeders' figures at the same printed size; and the template's labels 33 and 35 are swapped,
# since the deck's branches are 28-33, 33-34 and 29-35.
import csv, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DECK = HERE.parents[2] / "test_networks" / "ieee123_5poi_1ph"
SUBC = ["#1F3B73", "#E3A21A", "#169C78", "#874FAF"]
SUBBUS = ["1", "33", "120", "116"]
# Substation icons as in network_schematic.py: radius and name size in printed mm, at the figure's print scale.
SCALE = 88.9 / 400.75                               # one IEEE column over the drawing's width (mm per mm)
R_PRINT, NAME_PRINT = 1.5, 3.1
NAMEAT = [(-1.4, 2.4, "end"), (-1.2, -1.6, "end"), (1.4, -1.1, "start"), (0.0, 2.9, "middle")]   # name offsets, in radii
PX = 96 / 25.4                                      # px per mm (the template's user unit)

def block(s, gid):
    "Start and end of the <g> element with this id."
    st = s.rfind("<g", 0, s.index(f'id="{gid}"')); d = 0
    for m in re.finditer(r"<g\b|</g>", s[st:]):
        d += 1 if m.group(0) == "<g" else -1
        if d == 0:
            return st, st + m.end()

def fit(elems, query, center):
    "Least-squares px = a * local + b per axis, from the group's elements without their own transform."
    import numpy as np
    L, P = [], []
    for e in elems:
        if "transform=" in e: continue
        i = re.search(r'\sid="([^"]+)"', e).group(1)
        if i not in query: continue
        x, y, w, h = query[i]; L.append(center(e)); P.append((x + w / 2, y + h / 2))
    L, P = np.array(L), np.array(P)
    return [np.polyfit(L[:, k], P[:, k], 1) for k in (0, 1)]

def main(template, qfile, out):
    s = Path(template).read_text(encoding="utf-8")
    query = {r[0]: tuple(map(float, r[1:])) for r in csv.reader(open(qfile))}
    dots = {r["bus"]: (float(r["x"]), float(r["y"])) for r in csv.DictReader(open(HERE / "network_layout_ieee123.csv"))}
    grab = lambda f: sorted({re.match(r"[^.\s]+", b).group(0) for b in re.findall(r"(?i)bus1=(\S+)", (DECK / f).read_text())}, key=int)
    pv, bat = grab("PVSystem.dss"), grab("Storage.dss")
    num = lambda e, k: float(re.search(rf'\s{k}="([-\d.e]+)"', e).group(1))

    for gid, tag, buses in (("g1", "rect", pv), ("g109435", "ellipse", bat)):
        a, b = block(s, gid); g = s[a:b]
        elems = re.findall(rf"<{tag}\b[^>]*/>", g, flags=re.S)
        if tag == "rect":
            center = lambda e: (num(e, "x") + num(e, "width") / 2, num(e, "y") + num(e, "height") / 2)
        else:
            center = lambda e: (num(e, "cx"), num(e, "cy"))
        (ax, bx), (ay, by) = fit(elems, query, center)
        proto = next(e for e in elems if "transform=" not in e)
        new = []
        for k, bus in enumerate(buses):
            X, Y = (dots[bus][0] - bx) / ax, (dots[bus][1] - by) / ay
            e = re.sub(r'\sid="[^"]+"', f' id="{tag}-msopf-{bus}"', proto)
            e = re.sub(r'inkscape:label="[^"]*"', f'inkscape:label="{bus}"', e)
            if tag == "rect":
                e = re.sub(r'\sx="[^"]+"', f' x="{X - num(proto, "width") / 2:.4f}"', e)
                e = re.sub(r'\sy="[^"]+"', f' y="{Y - num(proto, "height") / 2:.4f}"', e)
            else:
                e = re.sub(r'\scx="[^"]+"', f' cx="{X:.4f}"', e); e = re.sub(r'\scy="[^"]+"', f' cy="{Y:.4f}"', e)
            new.append(e)
        head = g[:g.index(">") + 1]
        s = s[:a] + head + "".join(new) + "</g>" + s[b:]
        print(f"{gid}: {len(elems)} template markers replaced by {len(new)} ({tag}, scale {ax:.3f} px per unit)")

    # Substations: hide the template's ring at bus 1; draw the same icon as the large feeders' figures (concentric
    # rings in the substation's colour, bold name), at the same printed size: the 400 mm drawing prints one column
    # (88.9 mm) wide, so R_PRINT mm on paper is R_PRINT / SCALE mm here. A top-level layer in document mm.
    a, b = block(s, "g32466"); g = s[a:b]
    s = s[:a] + g.replace('style="display:inline"', 'style="display:none"', 1) + s[b:]
    R, fs = R_PRINT / SCALE, NAME_PRINT / SCALE
    icons = ['<g inkscape:groupmode="layer" inkscape:label="substations" id="layer-msopf-substations">']
    for k, bus in enumerate(SUBBUS):
        x, y = dots[bus][0] / PX, dots[bus][1] / PX; c = SUBC[k]
        icons.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{R:.3f}" style="fill:none;stroke:{c};stroke-width:{0.16 * R:.3f}"/>')
        icons.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{0.62 * R:.3f}" style="fill:none;stroke:{c};stroke-width:{0.12 * R:.3f}"/>')
        ox, oy, anchor = NAMEAT[k]
        icons.append(f'<text x="{x + ox * R:.3f}" y="{y + oy * R + 0.35 * fs:.3f}" style="font-family:\'Times New Roman\';font-weight:bold;'
                     f'font-size:{fs:.3f}px;text-anchor:{anchor};fill:{c};stroke:#ffffff;stroke-width:{0.14 * fs:.3f};paint-order:stroke">S{k + 1}</text>')
    icons.append("</g>")
    s = s[:s.rindex("</svg>")] + "\n".join(icons) + "\n</svg>\n"

    # Labels 33 and 35: swap their texts.
    s = re.sub(r"(<tspan\b[^>]*>)(33|35)(</tspan>)", lambda m: m.group(1) + ("35" if m.group(2) == "33" else "33") + m.group(3), s)
    s = s.replace('sodipodi:docname="ieee123C_1ph.svg"', 'sodipodi:docname="ieee123.svg"')
    Path(out).write_text(s, encoding="utf-8")
    print(f"PV {pv}\nbatteries {bat}\nsubstations at {SUBBUS} -> {out}")

if __name__ == "__main__":
    main(*sys.argv[1:4])
