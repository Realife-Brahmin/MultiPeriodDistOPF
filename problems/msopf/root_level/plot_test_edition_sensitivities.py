"""Plot sensitivity predictions, with measured boundaries explicitly marked pending.

Usage: python problems/msopf/root_level/plot_test_edition_sensitivities.py <run-dir> [...]
Writes PNG, SVG and predicted daily-mean vertices alongside each run's evidence.
"""
import csv
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def read(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def plot(root):
    validation = json.loads((root / "validation.json").read_text())
    assert validation["validated"]
    sf, = root.glob("multi_sens_*.csv")
    pf, = root.glob("multi_predicted_*.csv")
    data, rays = read(sf), read(pf)
    names = sorted({r["substation"] for r in data})
    keys = [k for k in data[0] if k.startswith("K_")]
    K = np.array([[np.mean([float(r[k]) for r in data if r["substation"] == n]) for k in keys] for n in names])
    P = np.array([np.mean([float(r["import_kW"]) for r in data if r["substation"] == n]) for n in names])
    assert np.linalg.matrix_rank(K) == 3
    # Enumerate all intersections of three daily-mean no-backflow planes. Keep only
    # vertices satisfying all four planes, so a tetrahedron is verified, not assumed.
    verts, active = [], []
    for inds in itertools.combinations(range(4), 3):
        if np.linalg.cond(K[list(inds)]) > 1e10:
            continue
        v = np.linalg.solve(K[list(inds)], -P[list(inds)])
        if np.min(P + K @ v) >= -1e-5:
            verts.append(v)
            active.append([names[i] for i in inds])
    assert len(verts) == 4, f"Expected bounded daily-mean tetrahedron; got {len(verts)} vertices"
    verts = np.array(verts)
    with (root / "predicted_daily_vertices.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["vertex", *[k.removeprefix("K_").removesuffix("_kW_per_deg") + "_minus_subs1_deg" for k in keys], "zero_mean_import_planes"])
        for i, (v, faces) in enumerate(zip(verts, active), 1):
            w.writerow([i, *v, ";".join(faces)])
    idle = np.array([[float(r[k])*float(r["predicted_deg"]) for k in ("u1", "u2", "u3")] for r in rays if r["case"] == "nobatt"])
    blue, orange, ink = "#2a78d6", "#eb6834", "#303030"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(9, 7), facecolor="white")
    ax = fig.add_subplot(111, projection="3d")
    faces = [[verts[i] for i in c] for c in itertools.combinations(range(4), 3)]
    ax.add_collection3d(Poly3DCollection(faces, facecolor=blue, edgecolor=blue, alpha=0.08, linewidth=0.5))
    for i, j in itertools.combinations(range(4), 2):
        ax.plot(*verts[[i, j]].T, color=blue, linewidth=1.6)
    ax.scatter(*idle.T, facecolors="none", edgecolors=orange, s=18, linewidth=0.9, depthshade=False)
    ax.scatter([0], [0], [0], color=ink, marker="+", s=45)
    for i, v in enumerate(verts, 1):
        ax.text(*v, f"  V{i}", fontsize=9, color=ink)
    for axis, key in zip((ax.xaxis, ax.yaxis, ax.zaxis), keys):
        name = key.removeprefix("K_").removesuffix("_kW_per_deg")
        axis.set_label_text(f"δ{name.removeprefix('subs')} − δ1 (deg)")
        axis.labelpad = 6
        axis.set_major_locator(MaxNLocator(4))
        axis.set_tick_params(labelsize=9)
        axis.pane.fill = False
        axis._axinfo["grid"]["color"] = (0.87, 0.87, 0.87, 1)
    allpoints = np.vstack([verts, idle])
    lo, hi = allpoints.min(axis=0), allpoints.max(axis=0)
    for setter, l, h in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), lo, hi):
        setter(l-0.1*(h-l), h+0.1*(h-l))
    ax.set_box_aspect(np.maximum(hi-lo, 1e-3))
    title = "IEEE9500 with service transformers" if "9500" in root.name else ("SMART-DS P1U, Type II" if root.name.endswith("_typeii") else "SMART-DS P1U, Type III")
    fig.suptitle(title + " — four substations", fontsize=15, y=0.97)
    ax.set_title("Sensitivity predictions at T = 24; measured boundaries pending", fontsize=11, pad=16)
    fig.legend(handles=[Line2D([0], [0], color=blue, lw=1.6, label="Daily-mean import planes (cyclic-battery estimate)"),
                        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="none", markeredgecolor=orange,
                               label="Worst-hour import estimate (idle batteries), 78 rays")],
               loc="lower center", bbox_to_anchor=(0.5, 0.055), frameon=False, ncol=1, fontsize=9)
    fig.text(0.5, 0.022, "PV reactive power dispatched; customer-only voltage band. Estimates omit voltage-boundary cuts.", ha="center", fontsize=9, color=ink)
    fig.subplots_adjust(left=0.0, right=0.94, bottom=0.24, top=0.88)
    for suffix in ("png", "svg"):
        fig.savefig(root / f"predicted_region.{suffix}", dpi=180)
    plt.close(fig)
    print(f"Wrote predicted_region.png/.svg and four verified daily-mean vertices: {root.name}")


if __name__ == "__main__":
    for directory in sys.argv[1:]:
        plot(Path(directory))
