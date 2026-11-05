import csv
import statistics
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(r"C:/Users/Aryan Ritwajeet Jha/Documents/documents_general/MultiPeriodDistOPF/.codex/worktrees/pr171-tree-ablation/problems/mpopf/spar_ddp")
out = Path(__file__).parent
samples = defaultdict(list)
with (root / "results/proper_logging_campaign_oct16/near_opt_completed.csv").open() as f:
    for row in csv.DictReader(f):
        if (row["system"], row["solver"], row["schedule"], row["gap_percent"], row["metric"]) == ("large10kC_1ph", "ddp", "balanced", "0.1", "original"):
            samples[int(row["T"])].append(float(row["full_seconds"]))
assert set(samples) == {24, 96, 144, 288} and all(len(v)==2 for v in samples.values())
series = {1020: {t: statistics.median(v) for t,v in samples.items()}}
with (root / "experiments/battery_count/results_oct21/medians.csv").open() as f:
    for row in csv.DictReader(f):
        if row["method"] == "ddp":
            batteries = {"10k_Y": 100, "10k_X": 250}[row["variant"]]
            series.setdefault(batteries, {})[int(row["T"])] = float(row["full_optimizer_s"])
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), layout="constrained")
for b in (100, 250, 1020):
    horizons = sorted(series[b])
    axes[0].plot(horizons, [series[b][t] for t in horizons], "o-", label=f"{b:,} batteries")
for t in sorted(samples):
    batteries = sorted(series)
    axes[1].plot(batteries, [series[b][t] for b in batteries], "o-", label=f"T = {t}")
for ax, label in zip(axes, ("Horizon (periods)", "Battery count")):
    ax.set_xlabel(label)
    ax.set_ylabel("Full-convergence optimizer time (s)")
    ax.set_yscale("log")
    ax.grid(alpha=.22)
    ax.legend(frameon=False)
fresh = {12: 18.518, 48: 28.1705, 72: 42.027, 288: 199.038}
axes[0].scatter(list(fresh), list(fresh.values()), marker="D", color="black", zorder=4, label="1,020 batteries, fresh runs")
axes[0].legend(frameon=False)
axes[0].set_xticks([12, 48, 96, 144, 288])
axes[1].set_xticks([100,250,1020])
fig.suptitle("SPAR-DDP on large10k")
fig.text(.5, .005, "Two-repeat medians. 1,020 batteries: Oct 16. diamonds: fresh Nov runs. The Oct 21 variants have different installed storage.", ha="center", fontsize=8)
fig.get_layout_engine().set(rect=(0,.035,1,.965))
fig.savefig(out/"large10k-scaling.png", dpi=180)
with (out/"large10k-scaling.csv").open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["batteries","horizon","full_convergence_s","campaign"])
    for b in sorted(series):
        for t,y in sorted(series[b].items()):
            writer.writerow([b,t,y,"Oct16" if b==1020 else "Oct21"])

    for t,y in fresh.items():
        writer.writerow([1020,t,y,'Nov04 fresh'])
