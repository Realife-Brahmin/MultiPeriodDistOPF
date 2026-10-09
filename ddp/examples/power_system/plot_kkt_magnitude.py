"""Sparsity plots of the stage KKT matrix, full against thresholded by entry
magnitude (kkt_magnitude_analysis.jl pattern dumps), and a magnitude histogram.

For each system: panel 1 is every stored entry, panels 2-3 keep only entries
with |K_ij| >= 1e-4 and >= 1e-2 (diagonal always kept, as in the analysis),
panel 4 shows where the entries dropped at 1e-2 sit. Pixels count entries per
bin (log colour), so dense regions stay visible at large10k scale.

    python ddp/examples/power_system/plot_kkt_magnitude.py <pattern_dir> <figure_dir> <iter_by_system>
    e.g. ... captures/magnitude/patterns ddp/results/kkt_magnitude/figures ieee123C_1ph:67,ieee2522C_1ph:68,large10kC_1ph:100
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

patdir, figdir, spec = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(figdir, exist_ok=True)
LABEL = {"ieee123C_1ph": "ieee123", "ieee2522C_1ph": "med2522", "large10kC_1ph": "large10k"}

for item in spec.split(","):
    system, it = item.split(":")
    path = os.path.join(patdir, f"pattern_{system}_iter{it}.csv")
    if not os.path.exists(path):
        continue
    data = np.genfromtxt(path, delimiter=",", skip_header=1, usecols=(0, 1, 2, 3))
    i, j, lr, la = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    n = int(max(i.max(), j.max()))
    diag = i == j
    bins = min(400, n)
    panels = [("all entries", np.ones_like(la, bool)),
              ("|K_ij| >= 1e-4", diag | (la >= -4)),
              ("|K_ij| >= 1e-2", diag | (la >= -2)),
              ("dropped at 1e-2", ~diag & (la < -2))]
    fig, axes = plt.subplots(1, 5, figsize=(17, 3.6),
                             gridspec_kw={"width_ratios": [1, 1, 1, 1, 1.25]})
    for ax, (title, m) in zip(axes[:4], panels):
        h, _, _ = np.histogram2d(i[m], j[m], bins=bins, range=[[0.5, n + 0.5], [0.5, n + 0.5]])
        h = np.ma.masked_equal(h, 0)
        ax.imshow(h, origin="upper", cmap="viridis", norm=LogNorm(vmin=1, vmax=max(h.max(), 2)),
                  extent=[0.5, n + 0.5, n + 0.5, 0.5], interpolation="nearest")
        ax.set_title(f"{title}\n{int(m.sum()):,} entries ({100 * m.mean():.1f}%)", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
    ax = axes[4]
    ax.hist(la[~diag], bins=80, range=(-16, 12), color="C0", alpha=0.7, label="off-diagonal |K_ij|")
    ax.hist(la[diag], bins=80, range=(-16, 12), color="C3", alpha=0.6, label="diagonal |K_ii|")
    ax.hist(lr[~diag], bins=80, range=(-16, 12), histtype="step", color="k", label="off-diag, scaled r_ij")
    ax.set_yscale("log"); ax.set_xlabel("log10 magnitude", fontsize=8); ax.tick_params(labelsize=7)
    ax.legend(fontsize=6, loc="upper left")
    fig.suptitle(f"{LABEL.get(system, system)} stage-1 KKT, T=6, iteration {it} "
                 f"(n = {n:,}, {len(la):,} stored entries)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig(os.path.join(figdir, f"kkt_magnitude_{LABEL.get(system, system)}_iter{it}.png"), dpi=150)
    plt.close(fig)
    print("wrote", LABEL.get(system, system), it)
