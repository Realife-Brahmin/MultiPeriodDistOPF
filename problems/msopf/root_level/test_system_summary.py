# The MSOPF test systems at a glance (2026-10-29), for the paper's test-system table: buses, branches, service
# transformers, closed ties, substations, loads, capacitor banks, PV and batteries, read from the decks the studies use.
#   python problems/msopf/root_level/test_system_summary.py
# Topology as in network_schematic.py (the same decks, closed ties and substations); ratings from Loads.dss,
# PVSystem.dss (Pmpp), Storage.dss (kWrated, kWhrated) and Capacitor.dss (kvar; counted only where Master.dss redirects it). Writes results/multiperiod/test_system_summary.csv.
import csv, re
from pathlib import Path
from network_schematic import CONFIG, NET, deck

OUT = Path(__file__).resolve().parents[1] / "results" / "multiperiod" / "test_system_summary.csv"
num = lambda l, k: float(re.search(rf"(?i)\b{k}=([-\d.e]+)", l).group(1))
news = lambda d, f, kind: [l for l in (d / f).read_text().splitlines() if l.lower().startswith(f"new {kind}")]

rows = []
for system, cfg in CONFIG.items():
    d = NET / cfg["deck"]
    lines, pvs, bats, subs = deck(cfg)
    # A closed switch inside a substation (sf6937's sb9_p1uhs0_1247_267, joining hs0's two transformer buses) is a
    # line, not a tie between substations.
    lines = [(n, a, b, t and not n.startswith("sb")) for n, a, b, t in lines]
    buses = {b for _, a, b2, _ in lines for b in (a, b2)}
    pv, st, ld = news(d, "PVSystem.dss", "pvsystem"), news(d, "Storage.dss", "storage"), news(d, "Loads.dss", "load")
    redirected = any(l.strip().lower() == "redirect capacitor.dss" for l in (d / "Master.dss").read_text().splitlines())
    cp = news(d, "Capacitor.dss", "capacitor") if redirected else []
    rows.append(dict(system=system, buses=len(buses), branches=sum(not t for *_, t in lines),
                     service_transformers=sum(n.startswith("svc_") for n, *_ in lines), ties_closed=sum(t for *_, t in lines),
                     substations=len(subs), substations_added=len(cfg["added"]) + (3 if system == "ieee123" else 0),
                     loads=len(ld), load_kW=round(sum(num(l, "kW") for l in ld), 1),
                     capacitors=len(cp), capacitor_kvar=round(sum(num(l, "kvar") for l in cp), 1),
                     pv=len(pv), pv_kW=round(sum(num(l, "Pmpp") for l in pv), 1),
                     batteries=len(st), battery_kW=round(sum(num(l, "kWrated") for l in st), 1),
                     battery_kWh=round(sum(num(l, "kWhrated") for l in st), 1)))
    print(rows[-1])
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("wrote", OUT)
