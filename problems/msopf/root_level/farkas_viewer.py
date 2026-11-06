# Interactive page of the Farkas bracketing results (2026-10-27; rebuilt 2026-11-02 on the decks regenerated 2026-10-30):
# for every case and direction, the validated dispatch at 98% of the limiting hour's edge and the certified edge with
# batteries idle, and, where the whole-day campaign has run, the reach with batteries used over the day and its
# certified edge. Colours and layers follow the paper (PESGM2027 figures/colors.tex, scripts/mapstyle.py): green fill
# and crimson mesh for batteries idle, light blue fill and blue mesh for batteries used. The predicted polytope of the
# sensitivity screen is no longer drawn (user, 2026-11-02).
#   python problems/msopf/root_level/farkas_viewer.py
# Writes results/multiperiod/viewer/farkas_maps.html (self-contained apart from plotly); docs/farkas_maps.html is a copy.
# Inputs, all under results/kamiak/:
#   med3980, sf6937: farkas_oct30b/<prefix>_s<digits>|<prefix>_four/fb_A_*.csv, fb_C_*.csv (the regenerated decks)
#   ieee123 (deck unchanged since): farkas_oct30b/ieee123_s14, farkas_oct29/ieee123_two (S3+S4), farkas_oct26/ieee123_three,
#     farkas_oct26/ieee123_four
#   batteries used: farkas_disp_nov02/{ieee123_s14, ieee123_two, med3980_s123}, farkas_disp_oct30b/sfp1u_four (fb_AD_*, fb_CD_*)
from pathlib import Path
from itertools import combinations
import csv, json, math, os

ROOT = Path(__file__).resolve().parents[1] / "results"
K = ROOT / "kamiak"
OUT = ROOT / "multiperiod" / "viewer" / "farkas_maps.html"
LIMIT = 20.0                                         # deg, the search limit of the region campaigns

def rows(p):
    with open(p, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def num(x):
    return float(x) if x not in (None, "", "NaN") else None

def sphere_triangles(U):
    "Triangulation of unit vectors U (n x 3) on the sphere: the facets of their convex hull."
    out = []
    for t in combinations(range(len(U)), 3):
        p = [U[k] for k in t]
        a = [p[1][j] - p[0][j] for j in range(3)]; b = [p[2][j] - p[0][j] for j in range(3)]
        n = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
        if sum(v * v for v in n) < 1e-24:
            continue
        s = [sum(n[j] * (U[k][j] - p[0][j]) for j in range(3)) for k in range(len(U))]
        if all(v <= 1e-9 for v in s) or all(v >= -1e-9 for v in s):
            out.append(list(t))
    return out

NOTES = {"ieee123": "IEEE 123-node feeder at its own 4.16 kV, each substation behind the feeder's 5 MVA 1+j8% transformer; ieee123C's PV and batteries (44% of load); lines at their conductors' ampacity (336 ACSR 530 A, 1/0 ACSR 230 A, 1/0 AA cable 195 A, x k/3); the other substations out; band at every non-source bus.",
         "med3980": "IEEE 9500-node network with its service transformers (3980 buses), ties closed; substation 4 is the added one at l2804256; dropped substations' sections fed through the ties; band at customer buses.",
         "sf6937": "Synthetic San Francisco feeder (NREL SMART-DS P1U), Type III ties; substations 1-4 are hs2, hs1, hs3, hs0; its own SMART-DS customer loads on the dataset's peak day (28 May) and PV from that day's sunrise to sunset; band at customer buses."}
PREFIX = {"med3980": "med3980", "sf6937": "sfp1u"}
DISP = {("ieee123", "14"): K / "farkas_disp_nov02" / "ieee123_s14", ("ieee123", "34"): K / "farkas_disp_nov02" / "ieee123_two",
        ("med3980", "123"): K / "farkas_disp_nov02" / "med3980_s123", ("sf6937", "1234"): K / "farkas_disp_oct30b" / "sfp1u_four"}

# The paper's three cases with rated line currents (2026-11-04): ieee123 S1+S3 (4.16 kV), med3980 S1-S3, sf6937 all four
# (day 148); every line at its rating, the idle certificates with derived current bounds and lsec, all exact-checked.
# Per case: (folder of fb_A and fb_AD, glob of the idle certificates, folder of fb_CD). FV_ALL=1 draws the earlier
# unrated pages instead (every pair and triple; the folders in the header).
RATED = {("ieee123", "13"): (ROOT / "ieee123_416_nov03" / "fleetC_tiers" / "ieee123_s13",
                             ROOT / "ieee123_416_nov03" / "derived_bounds" / "ieee123_s13", "fb_C_*.csv",
                             ROOT / "ieee123_416_nov03" / "derived_bounds" / "ieee123_s13"),
         ("med3980", "123"): (K / "farkas_rated_nov02" / "med3980_s123", K / "farkas_rated_nov02" / "med3980_s123_derived", "fb_C_*.csv",
                              K / "farkas_rated_nov02" / "med3980_s123"),
         ("sf6937", "1234"): (K / "farkas_day148b_nov03" / "sfp1u_four", K / "farkas_day148b_nov03" / "sfp1u_four", "fb_Cderived_*.csv",
                              K / "farkas_day148b_nov03" / "sfp1u_four")}

def cases():
    yield "ieee123", "14", K / "farkas_oct30b" / "ieee123_s14"
    yield "ieee123", "34", K / "farkas_oct29" / "ieee123_two"
    yield "ieee123", "134", K / "farkas_oct26" / "ieee123_three"
    yield "ieee123", "1234", K / "farkas_oct26" / "ieee123_four"
    for sysname, prefix in PREFIX.items():
        for n in (2, 3):
            for c in combinations("1234", n):
                yield sysname, "".join(c), K / "farkas_oct30b" / f"{prefix}_s{''.join(c)}"
        yield sysname, "1234", K / "farkas_oct30b" / f"{prefix}_four"

ALL = os.environ.get("FV_ALL") == "1"

def case(sysname, digits, d, cdir = None, cglob = "fb_C_*.csv", dd = None):
    A = {int(r["id"]): r for f in d.glob("fb_A_*.csv") for r in rows(f)}
    if not A:
        return None
    if all(num(a["day_reach_deg"]) is None for a in A.values()):
        print(f"{sysname}_s{digits}: no dispatch even at equal angles (the remaining substations cannot serve the day); left out")
        return None
    C = {int(r["id"]): r for f in (cdir or d).glob(cglob) for r in rows(f)}
    AD = {int(r["id"]): r for f in dd.glob("fb_AD_*.csv") for r in rows(f)} if dd else {}
    CD = {int(r["id"]): r for f in dd.glob("fb_CD_*.csv") for r in rows(f)} if dd else {}
    n = len(digits); sources = [f"subs{k}" for k in digits]
    dirs = []
    for i in sorted(A):
        a = A[i]; u = [float(a[f"u{k}"]) for k in range(1, n)]; day = float(a["day_reach_deg"])
        x = dict(id=i, u=u, day=day)
        if n == 2:
            x["label"] = ("+" if u[0] > 0 else "−") + f"(δ{digits[1]} − δ{digits[0]})"
        elif n == 3:
            x["theta"] = round(math.degrees(math.atan2(u[1], u[0])) % 360, 1)
        c = C.get(i)
        if a["at_limit"] == "true":
            x.update(status="limit", inner=LIMIT, edge=LIMIT, cert=LIMIT, hour=None)
        elif c is not None and num(c["r_cert_deg"]) is not None:
            r = float(c["hour_boundary_deg"]); rc = float(c["r_cert_deg"])
            x.update(status="certified", inner=0.98 * r, edge=r, cert=rc, hour=int(c["tstar"]), band=100 * (rc / r - 1),
                     farkas=num(c["farkas_at_cert"]) is not None, relax=int(c["relaxations"]), secs=float(c["seconds"]))
        else:
            x.update(status="pending", inner=0.98 * day, edge=day, cert=day, hour=None)
        ad = AD.get(i)
        if ad is not None:
            lim = ad["at_limit"] == "true"
            x.update(dreach=LIMIT if lim else float(ad["reach_deg"]), dlimit=lim)
            cd = CD.get(i)
            if cd is not None and num(cd["r_cert_deg"]) is not None:
                x.update(dcert=float(cd["r_cert_deg"]), dunres=float(cd["r_unresolved_deg"]),
                         dfarkas=num(cd["farkas_at_cert"]) is not None, dsecs=float(cd["seconds"]))
        dirs.append(x)
    combo = "+".join(digits)
    s = dict(id=f"{sysname}_s{digits}", sys=sysname, n=n, combo=combo, title=f"{sysname}, substations {', '.join(digits)}", dim=n - 1,
             dirs=dirs, sources=sources, note=NOTES[sysname], disp=bool(AD),
             axes=[f"δ{k} − δ{digits[0]} (°)" for k in digits[1:]])
    if n == 4:
        s["tri"] = sphere_triangles([x["u"] for x in dirs])
    return s

data = []
todo = ([(s_, g, d, None, "fb_C_*.csv", DISP.get((s_, g))) for s_, g, d in cases()] if ALL else
        [(s_, g, a, c, cg, cd) for (s_, g), (a, c, cg, cd) in RATED.items()])
for sysname, digits, d, cdir, cglob, dd in todo:
    s = case(sysname, digits, d, cdir, cglob, dd)
    if s is None:
        continue
    data.append(s); ds = s["dirs"]
    print(f"{s['id']}: idle {sum(x['status'] == 'certified' for x in ds)} certified, {sum(x['status'] == 'limit' for x in ds)} at the limit, "
          f"{sum(x['status'] == 'pending' for x in ds)} pending" +
          (f"; batteries used {sum('dcert' in x for x in ds)} of {len(ds)} certified, {sum(x.get('dlimit', False) for x in ds)} at the limit" if s["disp"] else ""))

data.sort(key=lambda s: (["ieee123", "med3980", "sf6937"].index(s["sys"]), s["n"], s["combo"]))
html = (Path(__file__).with_name("farkas_viewer_template.html").read_text(encoding="utf-8")
        .replace("__DATA__", json.dumps(data, separators=(",", ":"))).replace("__LIMIT__", str(LIMIT)))
OUT.write_text(html, encoding="utf-8")
print("wrote", OUT)
