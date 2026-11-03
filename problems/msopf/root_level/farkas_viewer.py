# Interactive page of the four-source Farkas bracketing results (2026-10-27): for every direction, the validated
# dispatch at 98% of the hour's edge, the measured edge and the certified edge, drawn as an inner (proven feasible)
# and an outer (proven infeasible beyond) shell over the predicted polytope of the sensitivity screen.
#   python problems/msopf/root_level/farkas_viewer.py
# Writes results/multiperiod/viewer/farkas_maps.html (self-contained apart from plotly). Same look as four_regions.html,
# kept as a separate page so as not to touch that viewer's uncommitted edits.
# Inputs: results/kamiak/farkas_oct26/<system>/fb_A_*.csv, fb_C_*.csv; results/multiperiod/multi_windows_fib78.csv;
# the T = 24 hourly sensitivities multi_sens_<stage>_T24.csv of the regions campaign (results/kamiak/regions_oct25).
from pathlib import Path
from itertools import combinations
import csv, json, math, os
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "results"
FB = ROOT / "kamiak" / "farkas_oct26"
OUT = ROOT / "multiperiod" / "viewer" / "farkas_maps.html"
LIMIT = 20.0                                         # deg, the search limit of the region campaign
SYSTEMS = [  # (tab id, title, bracketing folder, multi_windows system, sensitivity stage, note)
    ("ieee123x4", "ieee123, four substations", "ieee123_four", "ieee123_5poi_1ph", "four_ieee123a",
     "IEEE 123-node feeder, substations 1-4 (the paper's four-source case); band at every non-source bus."),
    ("sf6937", "sf6937, four substations", "sfp1u_four", "smartds_p1u_xf_typeIII", "four_smartds_p1u_typeiii",
     "Synthetic San Francisco feeder (NREL SMART-DS, urban sub-region P1U), service transformers, band at customer buses."),
    ("med3980", "med3980, four substations", "med3980_four", "ieee9500xf_1ph", "four_ieee9500xf",
     "IEEE 9500-node network with its service transformers (3980 buses), three substations plus one added; band at customer buses."),
]

def rows(p):
    with open(p, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def sens(stage):
    f = next(ROOT.joinpath("kamiak").rglob(f"multi_sens_{stage}_T24.csv"))
    out = []
    for r in rows(f):
        ks = [float(v) for k, v in r.items() if k.startswith("K_")]
        out.append((r["substation"], int(r["hour"]), float(r["import_kW"]), ks))
    return out

def polytope(hs):
    "Faces of {d: P + K d >= 0 for every substation and hour}, as {substation: [triangles of 3-D points]}."
    A = np.array([[-k for k in h[3]] for h in hs]); b = np.array([h[2] for h in hs])     # A d <= b
    tri = np.array(list(combinations(range(len(hs)), 3)))
    M = A[tri]; det = np.linalg.det(M); ok = np.abs(det) > 1e-9 * np.abs(A).max() ** 3
    V = np.linalg.solve(M[ok], b[tri[ok]][..., None])[..., 0]
    V = V[(V @ A.T <= b + 1e-6 * np.abs(b).max()).all(axis=1)]
    V = V[np.unique(np.round(V, 7), axis=0, return_index=True)[1]]   # dedupe, but keep the unrounded vertices for the on-face test
    faces = {}
    for i, h in enumerate(hs):
        on = V[np.abs(V @ A[i] - b[i]) <= 1e-6 * max(1.0, abs(b[i]))]
        if len(on) < 3:
            continue
        c = on.mean(axis=0); n = A[i] / np.linalg.norm(A[i])
        e1 = on[0] - c; e1 -= n * (e1 @ n); e1 /= np.linalg.norm(e1); e2 = np.cross(n, e1)
        on = on[np.argsort(np.arctan2((on - c) @ e2, (on - c) @ e1))]
        faces.setdefault(h[0], []).extend([[c.tolist(), on[j].tolist(), on[(j + 1) % len(on)].tolist()] for j in range(len(on))])
    return faces

def sphere_triangles(U):
    "Triangulation of unit vectors U (n x 3) on the sphere: the facets of their convex hull."
    out = []
    for t in combinations(range(len(U)), 3):
        p = U[list(t)]; n = np.cross(p[1] - p[0], p[2] - p[0])
        if np.linalg.norm(n) < 1e-12:
            continue
        s = (U - p[0]) @ n
        if (s <= 1e-9).all() or (s >= -1e-9).all():
            out.append(list(t))
    return out

def three_source():
    "ieee123 with substations 1, 3, 4 (2-D): the 16 bracketed directions and the predicted simplex."
    C = sorted((r for f in (FB / "ieee123_three").glob("fb_C_*.csv") for r in rows(f)), key=lambda r: int(r["id"]))
    P = [r for r in rows(ROOT / "multiperiod" / "simplex_check.csv") if (r["study"], r["pv"], r["case"]) == ("three", "qfree", "nobatt")]
    dirs = []
    for c, p in zip(C, P):
        th = math.radians(float(p["direction"])); r = float(c["hour_boundary_deg"]); rc = float(c["r_cert_deg"])
        dirs.append(dict(id=int(c["id"]), u=[math.cos(th), math.sin(th)], theta=float(p["direction"]), status="certified",
                         inner=0.98 * r, edge=r, cert=rc, pred=float(p["predicted_deg"]), hour=int(c["tstar"]), day=float(p["measured_deg"]),
                         band=100 * (rc / r - 1), farkas=c["farkas_at_cert"] not in ("", "NaN"), relax=int(c["relaxations"]), secs=float(c["seconds"])))
    S = [r for r in rows(ROOT / "multiperiod" / "simplex_figure_data.csv") if (r["pv"], r["case"], r["kind"]) == ("qfree", "nobatt", "predicted")]
    simplex = {"x": [float(r["d3_minus_d1_deg"]) for r in S], "y": [float(r["d4_minus_d1_deg"]) for r in S]}
    return dict(id="ieee123", sys="ieee123", n=3, combo="1+3+4", title="ieee123, substations 1+3+4", dim=2, dirs=dirs, simplex=simplex,
                axes=["δ₃ − δ₁ (°)", "δ₄ − δ₁ (°)"], sources=["subs1", "subs3", "subs4"],
                note="IEEE 123-node feeder, substations 1, 3 and 4; band at every non-source bus. Hour 12 limits every direction.")

# Two and three sources (2026-10-29, results/kamiak/farkas_oct29): no measured campaign, so the day's reach comes from
# phase A itself, and there is no predicted simplex.
# Every pair and triple of each four-source case (user, 2026-10-29): folder <prefix>_s<digits>, or the first batch's names.
DESC = {"ieee123": ("ieee123", "IEEE 123-node feeder; the other substations out; band at every non-source bus."),
        "med3980": ("med3980", "IEEE 9500-node network with its service transformers (3980 buses), ties closed; substation 4 is the added one at l2804256; dropped substations' sections fed through the ties; band at customer buses."),
        "sf6937": ("sfp1u", "Synthetic San Francisco feeder (NREL SMART-DS P1U), Type III ties; substations 1-4 are hs2, hs1, hs3, hs0; dropped substations' sections fed through the ties; band at customer buses.")}
ALIAS = {"ieee123_s34": "ieee123_two", "med3980_s12": "med3980_two", "med3980_s123": "med3980_three",
         "sfp1u_s14": "sfp1u_two", "sfp1u_s124": "sfp1u_three"}
SUBSETS = [(f"{sysname}_s{''.join(c)}", sysname, len(c), ALIAS.get(f"{prefix}_s{''.join(c)}", f"{prefix}_s{''.join(c)}"),
            [f"subs{k}" for k in c], note)
           for sysname, (prefix, note) in DESC.items() for n in (2, 3) for c in combinations("1234", n)
           if not (sysname == "ieee123" and c == ("1", "3", "4"))]          # that one is three_source(), from farkas_oct26
NAMES = {2: "two", 3: "three", 4: "four"}

def subset(tab, sysname, n, folder, sources, note):
    d = ROOT / "kamiak" / "farkas_oct29" / folder
    A = {int(r["id"]): r for f in d.glob("fb_A_*.csv") for r in rows(f)}
    if not A:
        return None
    C = {int(r["id"]): r for f in d.glob("fb_C_*.csv") for r in rows(f)}
    dirs = []
    for i in sorted(A):
        a = A[i]; u = [float(a[f"u{k}"]) for k in range(1, n)]; day = float(a["day_reach_deg"])
        x = dict(id=i, u=u, pred=None, day=day)
        if n == 2:
            x["label"] = ("+" if u[0] > 0 else "−") + f"({sources[1]} − {sources[0]})"
        else:
            x["theta"] = round(math.degrees(math.atan2(u[1], u[0])) % 360, 1)
        c = C.get(i)
        if a["at_limit"] == "true":
            x.update(status="limit", inner=LIMIT, edge=LIMIT, cert=LIMIT, hour=None)
        elif c is not None and c["r_cert_deg"] not in ("", "NaN"):
            r = float(c["hour_boundary_deg"]); rc = float(c["r_cert_deg"])
            x.update(status="certified", inner=0.98 * r, edge=r, cert=rc, hour=int(c["tstar"]), band=100 * (rc / r - 1),
                     farkas=c["farkas_at_cert"] not in ("", "NaN"), relax=int(c["relaxations"]), secs=float(c["seconds"]))
        else:
            x.update(status="pending", inner=0.98 * day, edge=day, cert=day, hour=None)
        dirs.append(x)
    ax = [f"δ{s.replace('subs', '')} − δ{sources[0].replace('subs', '')} (°)" for s in sources[1:]]
    combo = "+".join(s.replace("subs", "") for s in sources)
    return dict(id=tab, sys=sysname, n=n, combo=combo, title=f"{sysname}, substations {combo}", dim=n - 1, dirs=dirs, simplex=None,
                axes=ax, sources=sources, note=note)

SHOW = os.environ.get("FV_SYSTEMS", "ieee123,ieee123x4,sf6937,med3980,subsets").split(",")
data = [three_source()] if "ieee123" in SHOW else []
print(f"ieee123: {len(data[0]['dirs']) if data else 0} certified")
for case in SUBSETS:
    s = subset(*case) if "subsets" in SHOW else None
    if s is not None:
        data.append(s)
        print(f"{case[0]}: {sum(x['status'] == 'certified' for x in s['dirs'])} certified, {sum(x['status'] == 'limit' for x in s['dirs'])} at the limit, "
              f"{sum(x['status'] == 'pending' for x in s['dirs'])} pending")
for tab, title, folder, system, stage, note in SYSTEMS:
    d = FB / folder
    if tab not in SHOW or not d.exists():
        continue
    A = {int(r["id"]): r for f in d.glob("fb_A_*.csv") for r in rows(f)}
    C = {int(r["id"]): r for f in d.glob("fb_C_*.csv") for r in rows(f)}
    W = {int(r["direction"]): r for r in rows(ROOT / "multiperiod" / "multi_windows_fib78.csv") if r["system"] == system and r["sources"] == "subs1:subs2:subs3:subs4" and r["case"] == "nobatt"}
    dirs = []
    for i in sorted(W):
        w = W[i]; u = [float(w["u1"]), float(w["u2"]), float(w["u3"])]
        day, pred = float(w["measured_deg"]), float(w["predicted_deg"])
        if w["at_limit"] == "true":
            dirs.append(dict(id=i, u=u, status="limit", inner=LIMIT, edge=LIMIT, cert=LIMIT, pred=pred, hour=None, day=day)); continue
        c = C.get(i)
        if c is not None and c["r_cert_deg"] not in ("", "NaN"):
            r = float(c["hour_boundary_deg"]); rc = float(c["r_cert_deg"])
            dirs.append(dict(id=i, u=u, status="certified", inner=0.98 * r, edge=r, cert=rc, pred=pred, hour=int(c["tstar"]), day=day,
                             band=100 * (rc / r - 1), farkas=c["farkas_at_cert"] not in ("", "NaN"), relax=int(c["relaxations"]), secs=float(c["seconds"])))
        else:
            dirs.append(dict(id=i, u=u, status="pending", inner=0.98 * day, edge=day, cert=day, pred=pred, hour=None, day=day))
    U = np.array([x["u"] for x in dirs])
    data.append(dict(id=tab, sys=tab.replace("x4", ""), n=4, combo="1+2+3+4", title=title, note=note, dim=3, dirs=dirs, tri=sphere_triangles(U),
                     faces=polytope(sens(stage)), sources=["subs1", "subs2", "subs3", "subs4"]))
    print(f"{tab}: {sum(x['status'] == 'certified' for x in dirs)} certified, {sum(x['status'] == 'limit' for x in dirs)} at the limit, "
          f"{sum(x['status'] == 'pending' for x in dirs)} pending; predicted faces {sorted(data[-1]['faces'])}")

data.sort(key=lambda s: (["ieee123", "med3980", "sf6937"].index(s["sys"]), s["n"], s["combo"]))
html = (Path(__file__).with_name("farkas_viewer_template.html").read_text(encoding="utf-8")
        .replace("__DATA__", json.dumps(data, separators=(",", ":"))).replace("__LIMIT__", str(LIMIT)))
OUT.write_text(html, encoding="utf-8")
print("wrote", OUT)
