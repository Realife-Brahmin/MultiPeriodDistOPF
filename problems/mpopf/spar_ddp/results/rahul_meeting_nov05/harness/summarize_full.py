from pathlib import Path
import csv,json,statistics,re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
root=Path(__file__).parent
rows=[]
for campaign in ("phaseb_full_profile_nov04","rahul_ablation_nov05","rahul_fresh_288_nov05"):
    for folder in sorted((root/campaign/"results").iterdir()):
        if not (folder/"COMPLETE").exists():
            raise RuntimeError(f"Incomplete {folder}")
        m=re.fullmatch(r"(10k_Y|10k_1020)_T(\d+)_ddp_r(.+)",folder.name)
        assert m
        variant,T,rep=m.groups()
        summary=json.loads((folder/"timing_summary.json").read_text())
        with (folder/"passes.csv").open() as f:
            passes=list(csv.DictReader(f))
        with (folder/"iteration_timing.csv").open() as f:
            iterations=list(csv.DictReader(f))
        assert len(passes)==len(iterations)>0
        physical=(folder/"physical.log").read_text()
        primal=float(re.search(r"original_max=(\S+)",physical)[1])
        assert primal<=1e-5
        dss=(folder/"opendss.log").read_text()
        assert re.search(r"below limit: 0",dss) and re.search(r"above limit: 0",dss)
        row=dict(campaign=campaign,variant=variant,T=int(T),repeat=rep,
                 profiled=(campaign=="phaseb_full_profile_nov04" and rep=="3"),
                 passes=len(passes),final_iteration=passes[-1]["iteration"],
                 full_s=summary["full_optimizer_s"],A_s=summary["phase_A_s"],
                 B_s=summary["scan_wall_s"],C_s=summary["phase_C_s"],
                 original_primal=primal,objective=float(re.search(r"objective=(\S+)",physical)[1]),
                 B_per_pass_s=summary["scan_wall_s"]/len(passes),
                 gc_worker_sum_s=sum(float(p["gc_sum_s"]) for p in passes))
        row["other_s"]=row["full_s"]-sum(row[k] for k in ("A_s","B_s","C_s"))
        row["B_fraction"]=row["B_s"]/row["full_s"]
        rows.append(row)
with (root/"full_runs.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
fig,axes=plt.subplots(1,2,figsize=(10,4.8),layout="constrained")
labels=["100 batteries\nT=96","1,020 batteries\nT=96","100 batteries\nT=144","1,020 batteries\nT=144"]
groups=[("10k_Y",96),("10k_1020",96),("10k_Y",144),("10k_1020",144)]
bottom=[0]*4
for key,label,color in (("A_s","Phase A","#3274a1"),("B_s","Phase B","#e1812c"),("C_s","Phase C","#3a923a"),("other_s","Other solve work","#9a9a9a")):
    values=[statistics.median(r[key] for r in rows if r["campaign"]=="phaseb_full_profile_nov04" and not r["profiled"] and (r["variant"],r["T"])==g) for g in groups]
    axes[0].bar(range(4),values,bottom=bottom,label=label,color=color)
    bottom=[a+b for a,b in zip(bottom,values)]
axes[0].set_xticks(range(4),labels,fontsize=8)
axes[0].set_ylabel("Full-convergence optimizer time (s)")
axes[0].legend(frameon=False)
points={}
for r in rows:
    if r["variant"]=="10k_1020" and r["repeat"] in ("1","2") and not r["profiled"]:
        points.setdefault(r["T"],[]).append(r["full_s"])
xs=sorted(points)
axes[1].plot(xs,[statistics.median(points[t]) for t in xs],"o-")
axes[1].set_xticks(xs)
axes[1].set_xlabel("Horizon (periods)")
axes[1].set_ylabel("Full-convergence optimizer time (s)")
axes[1].grid(alpha=.2)
fig.suptitle("Fresh large10k solves — two-repeat medians")
fig.savefig(root/"fresh_full_solve_profile.png",dpi=180)
print("Validated",len(rows),"full solves, including four separately sampled profiles.")
