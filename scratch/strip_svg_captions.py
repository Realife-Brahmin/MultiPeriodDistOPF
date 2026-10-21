"""Drop the generator's in-figure caption lines from a hand-edited SVG.

Targets the two <text> elements by id so every hand-added element -- the
substation rings, the '1' bus label, the F1..F102 feeder labels -- survives
untouched. Textual removal rather than an XML round-trip, because ElementTree
would rewrite the sodipodi/inkscape namespace prefixes on the way out.
"""
import re
import sys

src, dst = sys.argv[1], sys.argv[2]
ids = sys.argv[3].split(",")

svg = open(src, encoding="utf-8").read()
before = len(svg)

for i in ids:
    pat = re.compile(r'[ \t]*<text\b[^>]*?\bid="' + re.escape(i) + r'"[^>]*?>.*?</text>\s*',
                     re.DOTALL)
    svg, n = pat.subn("", svg)
    print(f"  {i}: removed {n} element(s)")
    if n != 1:
        sys.exit(f"expected exactly 1 match for id={i}, got {n}")

open(dst, "w", encoding="utf-8").write(svg)
print(f"  {before} -> {len(svg)} bytes")

# Anything still referencing the stripped captions would mean a missed node.
for probe in ("102 feeders", "10,321 buses"):
    if probe in svg:
        print(f"  WARNING: '{probe}' still present in output")
