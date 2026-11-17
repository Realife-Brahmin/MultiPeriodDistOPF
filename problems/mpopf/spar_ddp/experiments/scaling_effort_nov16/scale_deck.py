"""Scale a single-substation deck's loads and DER fleets by lambda.

    python scale_deck.py <variant deck dir> <source deck dir> <lambda>

Writes <variant>/Loads.dss, Capacitor.dss, PVSystem.dss and Storage.dss from the source deck's
files with every load's kW and kVAr, every capacitor's kvar, every PV unit's kVA and Pmpp, and every
battery's kVA, kWrated and kWhrated multiplied by lambda (so each fleet stays at
the same share of the load, the DER placement rule), and points the variant's
Master.dss at those copies. Lines, ties and the source are untouched.
"""
import re, sys, pathlib

var, src, lam = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), float(sys.argv[3])
KEYS = {"Loads.dss": ("kW", "kVAr"), "Capacitor.dss": ("kvar",), "PVSystem.dss": ("kVA", "Pmpp"), "Storage.dss": ("kVA", "kWrated", "kWhrated")}
for name, keys in KEYS.items():
    out = [f"! {var.name}: {src.name}/{name} with {', '.join(keys)} times {lam:g} (scale_deck.py)\n"]
    n = 0
    for line in (src / name).read_text().splitlines(keepends=True):
        if line.lstrip().lower().startswith("new "):
            for k in keys:
                pat = re.compile(r"(?i)(\s%s=)([-+0-9.eE]+)" % re.escape(k))
                line, c = pat.subn(lambda m: m.group(1) + repr(float(m.group(2)) * lam), line)
                n += c
        out.append(line)
    (var / name).write_text("".join(out), newline="\n")
    print(f"{name}: {n} values scaled by {lam:g}")
m = (var / "Master.dss").read_text()
for name in KEYS:
    m = re.sub(r"(?im)^Redirect \.\./%s/%s\s*$" % (re.escape(src.name), re.escape(name)), f"Redirect {name}", m)
m = re.sub(r"(?im)^! Load scale: .*\n", "", m)
m = m.replace("\nNew Circuit", f"\n! Load scale: loads, PV and batteries are the source deck's times {lam:g} (scale_deck.py).\nNew Circuit", 1)
(var / "Master.dss").write_text(m, newline="\n")
print("Master.dss redirects:", [l for l in m.splitlines() if l.startswith("Redirect")])
