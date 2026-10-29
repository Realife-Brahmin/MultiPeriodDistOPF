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
R_PRINT, NAME_PRINT = 2.0, 3.6
GHOST, TABLE_PRINT, INK = 2.6, 2.5, "#333333"         # fade radius (in ring radii); table font (printed mm); table ink
NAMEAT = [(-1.4, 2.4, "end"), (-2.15, -2.55, "end"), (0.0, -2.15, "middle"), (1.2, -2.15, "end")]   # name offsets, in radii, clear of the bus labels
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
    R, fs, tf = R_PRINT / SCALE, NAME_PRINT / SCALE, TABLE_PRINT / SCALE
    ring = lambda x, y, r, c: [
        f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{r:.3f}" style="fill:none;stroke:#ffffff;stroke-width:{0.48 * r:.3f}"/>',
        f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{r:.3f}" style="fill:none;stroke:{c};stroke-width:{0.2 * r:.3f}"/>',
        f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{0.6 * r:.3f}" style="fill:none;stroke:{c};stroke-width:{0.16 * r:.3f}"/>']
    text = lambda x, y, size, c, t, anchor="start", bold=False: (
        f'<text x="{x:.3f}" y="{y:.3f}" style="font-family:\'Times New Roman\';{"font-weight:bold;" if bold else ""}'
        f'font-size:{size:.3f}px;text-anchor:{anchor};fill:{c}">{t}</text>')
    # Each substation sits in a soft white disk, so the buses and lines right around it fade and the ring stands out.
    icons = ['<defs><radialGradient id="msopf-ghost"><stop offset="0" stop-color="#ffffff" stop-opacity="0.75"/>'
             '<stop offset="0.55" stop-color="#ffffff" stop-opacity="0.6"/><stop offset="1" stop-color="#ffffff" stop-opacity="0"/>'
             '</radialGradient></defs>',
             '<g inkscape:groupmode="layer" inkscape:label="substations" id="layer-msopf-substations">']
    for k, bus in enumerate(SUBBUS):
        x, y = dots[bus][0] / PX, dots[bus][1] / PX; c = SUBC[k]
        icons.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{GHOST * R:.3f}" style="fill:url(#msopf-ghost);stroke:none"/>')
        icons += ring(x, y, R, c)
        ox, oy, anchor = NAMEAT[k]
        icons.append(text(x + ox * R, y + oy * R + 0.35 * fs, fs, c, f"S{k + 1}", anchor, True)
                     .replace('">S', ';stroke:#ffffff;stroke-width:%.3f;paint-order:stroke">S' % (0.14 * fs)))
    # The substation buses in a one-row table under the drawing (its box: the query's first row, the document root,
    # in px), in the style of the large feeders' tables.
    x0, y0, w, h = [v / PX for v in next(iter(query.values()))]
    c0, cw = 8.0 / SCALE, 8.5 / SCALE                                   # first column, substation columns
    W = c0 + len(SUBBUS) * cw; tx, ty = x0 + (w - W) / 2, y0 + h + 2.5 / SCALE
    rule = lambda y, t: f'<line x1="{tx:.3f}" y1="{y:.3f}" x2="{tx + W:.3f}" y2="{y:.3f}" style="stroke:{INK};stroke-width:{t / SCALE:.3f}"/>'
    icons.append(rule(ty, 0.18)); y = ty + 1.05 * tf
    icons.append(text(tx + 0.25 * tf, y, tf, INK, "Sub."))
    for k in range(len(SUBBUS)):
        cx = tx + c0 + (k + 0.5) * cw
        icons += ring(cx - 0.55 * tf, y - 0.33 * tf, 0.6 * R, SUBC[k])
        icons.append(text(cx - 0.55 * tf + 0.6 * R + 0.3 * tf, y, tf, SUBC[k], f"S{k + 1}", "start", True))
    y += 0.4 * tf; icons.append(rule(y, 0.1)); y += 1.25 * tf
    icons.append(text(tx + 0.25 * tf, y, tf, INK, "Bus"))
    for k, bus in enumerate(SUBBUS):
        icons.append(text(tx + c0 + (k + 0.5) * cw, y, tf, INK, bus, "middle"))
    icons.append(rule(y + 0.45 * tf, 0.18))
    icons.append("</g>")
    s = s[:s.rindex("</svg>")] + "\n".join(icons) + "\n</svg>\n"

    # Labels 33 and 35: swap their texts.
    s = re.sub(r"(<tspan\b[^>]*>)(33|35)(</tspan>)", lambda m: m.group(1) + ("35" if m.group(2) == "33" else "33") + m.group(3), s)
    # The rings cover the substation buses' labels: blank them, the table under the drawing gives the buses.
    s = re.sub(r"(<tspan\b[^>]*>)(" + "|".join(SUBBUS) + r")(</tspan>)", r"\1\3", s)
    s = s.replace('sodipodi:docname="ieee123C_1ph.svg"', 'sodipodi:docname="ieee123.svg"')
    Path(out).write_text(s, encoding="utf-8")
    print(f"PV {pv}\nbatteries {bat}\nsubstations at {SUBBUS} -> {out}")

if __name__ == "__main__":
    main(*sys.argv[1:4])
