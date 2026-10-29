"""Check the matched primary/secondary controls and preserve the scenario-site selection."""
import collections
import hashlib
import heapq
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
NETWORKS = ROOT / "test_networks"


def definitions(deck, file):
    return [x.strip() for x in (NETWORKS/deck/file).read_text().splitlines() if x.lower().startswith("new ")]


def fields(line):
    return {k.lower(): v for k,v in re.findall(r"([%\w]+)=([^\s]+)",line)}


def main():
    primary="ieee8500_1ph"; secondary="ieee8500xf_1ph"
    a=definitions(primary,"BranchData.dss")
    b=definitions(secondary,"BranchData.dss")
    assert a == [l for l in b if not l.lower().startswith("new line.svc_")]
    for file in ("PVSystem.dss","Storage.dss"):
        assert definitions(primary,file)==definitions(secondary,file)
    services={}
    for line in b:
        if line.lower().startswith("new line.svc_"):
            d=fields(line)
            services[d["bus2"].split(".")[0].lower()]=d["bus1"].split(".")[0].lower()
    def loads(deck):
        out=collections.defaultdict(complex)
        for line in definitions(deck,"Loads.dss"):
            d=fields(line)
            out[d["bus"].split(".")[0].lower()]+=complex(float(d["kw"]),float(d["kvar"]))
        return out
    pl=loads(primary); sl=collections.defaultdict(complex)
    for bus,s in loads(secondary).items():sl[services[bus]]+=s
    discrepancy=max(abs(pl[bus]-sl[bus]) for bus in pl.keys()|sl.keys())
    assert discrepancy<1e-3
    adjacency=collections.defaultdict(list)
    for line in a:
        d=fields(line);u=d["bus1"].split(".")[0].lower();v=d["bus2"].split(".")[0].lower()
        z=complex(float(d["r1"]),float(d["x1"]))
        adjacency[u].append((v,abs(z)));adjacency[v].append((u,abs(z)))
    assert len(a)==len(adjacency)-1
    def distance(source):
        out={source:0};q=[(0,source)]
        while q:
            dist,u=heapq.heappop(q)
            if dist>out[u]:continue
            for v,w in adjacency[u]:
                if dist+w<out.get(v,float("inf")):
                    out[v]=dist+w;heapq.heappush(q,(dist+w,v))
        assert len(out)==len(adjacency)
        return out
    ds=[distance("hvmv_sub_hsb")];picked=[]
    candidates=[b for b in adjacency if len(adjacency[b])>=3 and not b.endswith("_cap")]
    for _ in range(3):
        bus=max(candidates,key=lambda b:(min(d[b] for d in ds),b))
        picked.append({"bus":bus,"min_path_distance_proxy_ohm":min(d[bus] for d in ds)})
        ds.append(distance(bus));candidates.remove(bus)
    assert [r["bus"] for r in picked]==["m1047401","m1108258","m1142880"]
    result={"validated":True,"primary_lines":len(a),"primary_buses":len(adjacency),
        "secondary_delivery_points":len(services),"primary_delivery_proxies":len(set(services.values())),
        "max_load_aggregation_discrepancy_kVA":discrepancy,"identical_DER_definitions":True,
        "distance_rule":"greedy max-min sum(abs(series_Z)) along the radial path, among primary junctions of degree >=3",
        "sites":picked,"new_source_series_equivalent_ohm":{"R":0.0759998,"X":0.877025},
        "voltage_policy":"0.95-1.05 pu at secondary delivery points or their primary attachment proxies; no transit-bus band",
        "controls":"PV q dispatched, fixed PV p, existing batteries, hard cyclic SOC; constant source EMF angles for 24 hours"}
    out=ROOT/"problems/msopf/results/secondary_8500_oct25/input_check.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+"\n")
    print("PASS: identical primary impedances and DER fleets; load aggregation differs by at most",discrepancy,"kVA")
    print("Scenario source sites:",[r["bus"] for r in picked])


if __name__=="__main__":main()
