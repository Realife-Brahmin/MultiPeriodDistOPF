"""Recolour the structural elements of a hand-edited SVG to a system ink.

ieee123C_1ph.svg and large10kC_1ph_top_level.svg carry hand drawing that the
generators cannot reproduce, so they must never be regenerated -- they are
edited in place instead.

Only buses and branches are touched. DER markers (gold PV, teal BESS), load
indicators and the substation keep their colours, because those encode meaning
rather than system identity. Substitutions are scoped to named groups so a
colour used elsewhere for something else is left alone.

Usage:
    python scripts/tint_svg_elements.py assets/networks/ieee123C_1ph.svg "#0e243c"
"""

import re
import sys
from pathlib import Path

# group id / inkscape:label -> list of (regex, replacement-template) applied
# only within that group. {ink} is the system ink; {d}/{l} its stub steps.
SCOPES = {
    # ieee123C_1ph.svg (hand-drawn)
    "branches": [(r"stroke:#000000", "stroke:{ink}")],
    "buses": [(r"fill:#000000", "fill:{ink}"), (r"stroke:#000000", "stroke:{ink}")],
    # large10kC_1ph_top_level.svg (generated star + hand-added rings/label)
    "panel-a-feeders": [(r'stroke="#666666"', 'stroke="{d}"'),
                        (r'stroke="#9a9a9a"', 'stroke="{l}"')],
    "panel-a-heads": [(r'fill="#000000"', 'fill="{ink}"')],
}


def _hex(rgb):
    return "#%02x%02x%02x" % tuple(max(0, min(255, int(round(v)))) for v in rgb)


def lighten(hex_colour, amount):
    r, g, b = (int(hex_colour[i:i + 2], 16) for i in (1, 3, 5))
    return _hex(tuple(c + (255 - c) * amount for c in (r, g, b)))


def find_group(svg, name):
    """Return (start, end) of the <g> whose id or inkscape:label is `name`."""
    for attr in ('id="%s"' % name, 'inkscape:label="%s"' % name):
        i = svg.find(attr)
        if i < 0:
            continue
        start = svg.rfind("<g", 0, i)
        depth, j = 0, start
        while j < len(svg):
            if svg.startswith("<g", j):
                depth += 1
            elif svg.startswith("</g>", j):
                depth -= 1
                if depth == 0:
                    return start, j + 4
            j += 1
    return None


path, ink = Path(sys.argv[1]), sys.argv[2]
svg = path.read_text(encoding="utf-8")
subs = {"ink": ink, "d": lighten(ink, 0.36), "l": lighten(ink, 0.68)}

total = 0
for name, rules in SCOPES.items():
    span = find_group(svg, name)
    if not span:
        continue
    start, end = span
    body = svg[start:end]
    for pattern, repl in rules:
        body, n = re.subn(pattern, repl.format(**subs), body)
        if n:
            print(f"  {name}: {n} x {pattern} -> {repl.format(**subs)}")
            total += n
    svg = svg[:start] + body + svg[end:]

if not total:
    sys.exit(f"no structural groups matched in {path.name}; nothing written")

path.write_text(svg, encoding="utf-8")
print(f"  wrote {path.name} ({total} substitutions)")
