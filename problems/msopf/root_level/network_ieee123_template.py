# The MSOPF ieee123 schematic from the SPAR-DDP Inkscape figure (2026-10-29): same layout, fonts, line widths
# and bus dots, with this deck's PV and batteries (test_networks/ieee123_5poi_1ph) and its four substations.
#   python problems/msopf/root_level/network_ieee123_template.py <template.svg> <query.csv> <out.svg>
# template: TPWRS-2027-Sparse-Parallel-DDP-for-MPOPF/figures/networks/ieee123C_1ph.svg
# query:    inkscape --query-all <template.svg> (object bounding boxes in px), used to map bus positions
#           (network_layout_ieee123.csv) into each marker group's own coordinates.
# Edits: the children of "pv-indicator-squares" and "battery-bus-circles" are redrawn at this deck's buses in the
# template's own marker styles; "Substation1_if_alone" is hidden and copied to buses 1, 33, 120 and 116 in the
# S1-S4 colours of the PESGM2027 figures, each with its name; and the template's labels 33 and 35 are swapped,
# since the deck's branches are 28-33, 33-34 and 29-35.
import csv, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DECK = HERE.parents[2] / "test_networks" / "ieee123_5poi_1ph"
SUBC = ["#1F3B73", "#E3A21A", "#169C78", "#874FAF"]
SUBBUS = ["1", "33", "120", "116"]
NAMEAT = [(-10, 15), (-14, -14), (12, -10), (11, -12)]   # where each name sits, from its ring's centre (local units)
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

    # Substations: hide the template's ring at bus 1, copy it to the four substation buses.
    a, b = block(s, "g32466"); g = s[a:b]
    x, y, w, h = query["g32466"]; cx, cy = x + w / 2, y + h / 2
    m = re.search(r'transform="matrix\(([^)]+)\)"', g).group(1)
    copies = []
    for k, bus in enumerate(SUBBUS):
        dx, dy = (dots[bus][0] - cx) / PX, (dots[bus][1] - cy) / PX
        c = g.replace('transform="matrix(', f'transform="translate({dx:.4f},{dy:.4f}) matrix(', 1)
        c = c.replace('id="g32466"', f'id="subs{k + 1}-msopf"').replace('inkscape:label="Substation1_if_alone"', f'inkscape:label="S{k + 1}"')
        c = re.sub(r'id="(path\d+)"', rf'id="\1-s{k + 1}"', c)
        c = c.replace("stroke:#2b2b2b", f"stroke:{SUBC[k]}").replace("stroke:#8f8f8f", f"stroke:{SUBC[k]}")
        c = c.replace("stroke-opacity:0.494118", "stroke-opacity:0.95").replace("opacity:0.83", "opacity:1")
        # its name, in the template's label font, in a free spot beside the ring (local units of the copy; ring centre
        # at (-48.07, 206.06), outer radius 9.5)
        ox, oy = NAMEAT[k]
        c = c.replace("</g>", f'<text x="{-48.07 + ox:.2f}" y="{206.06 + oy + 2.3:.2f}" style="font-family:\'Times New Roman\';'
                              f'font-weight:bold;font-size:6.5px;text-anchor:middle;fill:{SUBC[k]}" id="subs{k + 1}-name">S{k + 1}</text></g>', 1)
        copies.append(c)
    hidden = g.replace('style="display:inline"', 'style="display:none"', 1)
    s = s[:a] + hidden + "".join(copies) + s[b:]

    # Labels 33 and 35: swap their texts.
    s = re.sub(r"(<tspan\b[^>]*>)(33|35)(</tspan>)", lambda m: m.group(1) + ("35" if m.group(2) == "33" else "33") + m.group(3), s)
    s = s.replace('sodipodi:docname="ieee123C_1ph.svg"', 'sodipodi:docname="ieee123.svg"')
    Path(out).write_text(s, encoding="utf-8")
    print(f"PV {pv}\nbatteries {bat}\nsubstations at {SUBBUS} -> {out}")

if __name__ == "__main__":
    main(*sys.argv[1:4])
