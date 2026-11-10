from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

base=Path(__file__).resolve().parent
fig,axes=plt.subplots(1,2,figsize=(11,4.6),layout="constrained")
for ax,horizon,labels in ((axes[0],96,["Balanced\ncurrent 63x1","Grouped\ncurrent 63x1","Grouped\nright 63x1","Grouped\nright 31x2"]),(axes[1],288,["Balanced\ncurrent 63x1","Grouped\nright 31x2"])):
    with (base/f"results_T{horizon}/medians.csv").open() as handle:
        data=list(csv.DictReader(handle))
    bottom=[0.]*len(data)
    for field,label,color in (("A_s","Phase A","#3274a1"),("B_s","Phase B","#e1812c"),("C_s","Phase C","#3a923a"),("other","Other optimizer work","#999999")):
        values=[float(r[field]) if field!="other" else float(r["full_s"])-sum(float(r[k]) for k in ("A_s","B_s","C_s")) for r in data]
        ax.bar(range(len(data)),values,bottom=bottom,label=label,color=color)
        bottom=[a+b for a,b in zip(bottom,values)]
    ax.set_xticks(range(len(data)),labels,fontsize=8)
    ax.set_title(f"T = {horizon}")
    ax.set_ylabel("Full optimizer time after warm-up (s)")
    for i,r in enumerate(data):
        ax.text(i,float(r["full_s"])+1,f'{float(r["full_s"]):.2f}',ha="center",fontsize=9)
    ax.set_ylim(0,max(float(r["full_s"]) for r in data)*1.12)
axes[1].legend(frameon=False,fontsize=8,loc="upper right")
fig.suptitle("Original large10k, 1,020 batteries — two-repeat medians")
fig.savefig(base/"results_T288/comparison.png",dpi=180)
