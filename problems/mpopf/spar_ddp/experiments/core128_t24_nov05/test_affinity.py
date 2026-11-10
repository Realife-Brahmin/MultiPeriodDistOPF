"""Resource-map regression: 64/128 cores, disjoint teams, no SMT/socket crossing."""
import os,runpy,sys,tempfile
from pathlib import Path
from unittest.mock import patch
allocator=Path(sys.argv[1]).resolve()
for budget,w,b in [(64,24,1),(64,24,2),(128,63,1),(128,96,1),(128,48,2)]:
    with tempfile.TemporaryDirectory() as d:
        old=Path.cwd()
        os.chdir(d)
        topology="# CPU,CORE,SOCKET\n"+"".join(f"{i},{i},{i//64}\n" for i in range(128))
        try:
            with patch.dict(os.environ,{"STUDY_CORE_BUDGET":str(budget)}),patch.object(sys,"argv",[str(allocator),str(w),str(b)]),patch("os.sched_getaffinity",return_value=set(range(128))),patch("subprocess.check_output",return_value=topology):
                runpy.run_path(str(allocator),run_name="__main__")
            groups=[[int(x) for x in line.split(",")] for line in Path("cpu_groups.csv").read_text().splitlines()]
            assert len(groups)==w+1 and len(groups[0])==1
            assert all(len(g)==b and len({x//64 for x in g})==1 for g in groups[1:])
            flattened=[x for g in groups for x in g]
            assert len(set(flattened))==len(flattened)==w*b+1<=budget
            print(f"MAP_PASS budget={budget} workers={w} threads={b}")
        finally:
            os.chdir(old)

