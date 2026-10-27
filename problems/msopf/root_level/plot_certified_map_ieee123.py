# The certified three-source map of ieee123 from Farkas bracketing (2026-10-27): substations 1, 3, 4, batteries idle,
# T = 24, band at every non-source bus, current limits = the 1.0x envelope of the validated dispatches.
#   python problems/msopf/root_level/plot_certified_map_ieee123.py
# (a) the angle plane: the predicted simplex (sensitivities, no current limits), validated feasible points at 98% of
#     each direction's hour boundary, the hour boundary itself, and the certified edge beyond which the relaxation
#     proves the model with the limits infeasible.
# (b) the same along each direction, relative to the hour boundary.
# Inputs: simplex_figure_data.csv (prediction), results/kamiak/farkas_oct26/ieee123_three/fb_C_*.csv (bracketing).
from pathlib import Path
import math
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1] / "results"
RES = ROOT / "multiperiod"
FB = ROOT / "kamiak" / "farkas_oct26" / "ieee123_three"
pred = pd.read_csv(RES / "simplex_figure_data.csv")
pred = pred[(pred.pv == "qfree") & (pred.case == "nobatt") & (pred.kind == "predicted")]
win = pd.concat([pd.read_csv(f) for f in FB.glob("fb_C_*.csv")]).sort_values("id")
win["theta_deg"] = win.label.str.replace(" deg", "").astype(float)
tstar = sorted(win.tstar.unique())

th = [math.radians(t) for t in win.theta_deg]
def ring(r):
    xs = [ri * math.cos(t) for ri, t in zip(r, th)]; ys = [ri * math.sin(t) for ri, t in zip(r, th)]
    return xs + xs[:1], ys + ys[:1]

fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 4.6), gridspec_kw={"width_ratios": [1.05, 1]})
# (a) the plane
a.fill(*ring([1e3] * len(th)), color="#f6d5d1", zorder=0)          # everything far out, then cleared inside the certified edge
a.fill(*ring(win.r_cert_deg), color="#fff4cc", zorder=1)           # unresolved band
a.fill(*ring(0.98 * win.hour_boundary_deg), color="#d9efe0", zorder=2)
a.plot(pred.d3_minus_d1_deg, pred.d4_minus_d1_deg, color="#555555", lw=1.2, ls="--", zorder=3, label="predicted simplex")
x, y = ring(0.98 * win.hour_boundary_deg); a.plot(x, y, "o-", color="#2e8b57", ms=3.5, lw=1, zorder=4, label="validated feasible (98% of edge)")
x, y = ring(win.hour_boundary_deg); a.plot(x[:-1], y[:-1], "o", mfc="white", mec="black", ms=4, zorder=5, label="measured edge (local solves)")
x, y = ring(win.r_cert_deg); a.plot(x, y, "s-", color="#c0392b", ms=3, lw=1, zorder=6, label="certified infeasible beyond")
lim = 0.36
a.set_xlim(-lim, lim); a.set_ylim(-lim, lim); a.set_aspect("equal")
a.set_xlabel("δ₃ − δ₁ (deg)"); a.set_ylabel("δ₄ − δ₁ (deg)")
a.set_title("(a) ieee123, substations 1, 3, 4, idle batteries", fontsize=10)
a.legend(fontsize=7.5, loc="lower left", framealpha=0.9)
a.grid(alpha=0.25)
# (b) along each direction, relative to the measured edge
t = win.theta_deg; gap = 100 * (win.r_cert_deg / win.hour_boundary_deg - 1)
b.axhspan(0, 3, color="#f6d5d1", zorder=0)
b.fill_between(t, -2.0, gap, color="#fff4cc", zorder=1)
b.plot(t, [-2.0] * len(t), "o-", color="#2e8b57", ms=3.5, lw=1, label="validated feasible (98%)")
b.plot(t, gap, "s-", color="#c0392b", ms=3, lw=1, label="certified infeasible beyond")
b.axhline(0, color="black", lw=0.8, label="measured edge")
b.set_xlabel("direction θ (deg)"); b.set_ylabel("distance from the measured edge (%)")
b.set_xticks(range(0, 360, 45)); b.set_ylim(-2.6, 1.6)
b.set_title("(b) the unresolved band, per direction", fontsize=10)
b.legend(fontsize=7.5, loc="center right"); b.grid(alpha=0.25)
fig.text(0.5, 0.005, f"Hour {', '.join(map(str, tstar))} (the limiting hour in all {len(win)} directions); current limits = largest currents "
         "of the validated dispatches (1.0x). Certificates hold for the model with these limits.", ha="center", fontsize=7.5, color="#555555")
fig.tight_layout(rect=(0, 0.03, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(RES / f"fig_certified_map_ieee123.{ext}", dpi=200)
print("wrote", RES / "fig_certified_map_ieee123.png")
