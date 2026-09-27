"""Render binned KKT sparsity patterns from kkt_ordering_patterns.jl.

One row per configuration: the matrix after its ordering, its factor, and a zoom
on the factor's trailing rows and columns, where every ordering places the dense
battery-curvature block and where the fill differs. The captured K is shown once
at the top left. Colour is log10 of the nonzero count per bin; empty bins are
white, and any nonzero bin is at least light blue so isolated entries stay
visible at large10k scale.

    python ddp/examples/power_system/plot_kkt_ordering_patterns.py \
        <bins_dir> <tag> <out.png> "<title>" <config>[,<config>...]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

bins_dir, tag, out_png, title, configs = sys.argv[1:6]
configs = configs.split(",")
cmap = LinearSegmentedColormap.from_list("navy", ["#aebfe3", "#3a5a9c", "#0b1f4d"])
cmap.set_bad("#ffffff")


def load(name):
    path = Path(bins_dir) / f"{tag}__{name}.csv"
    header = path.open().readline().lstrip("# ").split()
    meta = dict(kv.split("=") for kv in header)
    counts = np.loadtxt(path, delimiter=",", skiprows=1)
    return counts, meta


def draw(ax, name, label):
    counts, meta = load(name)
    img = np.ma.masked_where(counts == 0, np.log10(np.maximum(counts, 1)))
    ax.imshow(img, cmap=cmap, interpolation="nearest", aspect="equal",
              vmin=0, vmax=max(1.0, float(img.max())))
    if "tail" in meta:
        sub = f"trailing {int(meta['tail']):,} rows/cols: nnz {int(meta['nnz']):,}"
    else:
        sub = f"nnz {int(meta['nnz']):,}   bandwidth {int(meta['bandwidth']):,}"
    ax.set_title(f"{label}\n{sub}", fontsize=8.5)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#8a8f99"); s.set_linewidth(0.6)


rows = len(configs)
fig, axes = plt.subplots(rows, 4, figsize=(13.5, 3.5 * rows), squeeze=False)
for r, c in enumerate(configs):
    if r == 0:
        draw(axes[r][0], "K", "captured K")
    elif r == 1:
        draw(axes[r][0], "K__tail", "captured K (zoom)")
    else:
        axes[r][0].axis("off")
    draw(axes[r][1], f"{c}__reordered", f"{c}: reordered")
    draw(axes[r][2], f"{c}__factor", f"{c}: factor")
    draw(axes[r][3], f"{c}__factor__tail", f"{c}: factor (zoom)")
fig.suptitle(title, fontsize=11)
fig.tight_layout()
fig.savefig(out_png, dpi=140)
print("wrote", out_png)
