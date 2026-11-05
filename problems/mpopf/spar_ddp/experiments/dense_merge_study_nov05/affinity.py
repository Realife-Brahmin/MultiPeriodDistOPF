import os,subprocess,sys
from collections import defaultdict
from pathlib import Path
workers=int(sys.argv[1])
threads=[int(b) for b in sys.argv[2:]]
allowed=os.sched_getaffinity(0)
sockets=defaultdict(list)
used=set()
text=subprocess.check_output(["lscpu","-p=CPU,CORE,SOCKET"],text=True)
for line in text.splitlines():
    if line.startswith("#"):continue
    cpu,core,socket=map(int,line.split(","))
    if cpu in allowed and (socket,core) not in used:
        used.add((socket,core));sockets[socket].append(cpu)
assert sum(map(len,sockets.values()))>=64
# Keep precisely 64 physical CPUs, taking an even share of each socket.
pool=defaultdict(list)
while sum(map(len,pool.values()))<64:
    for socket in sorted(sockets):
        if sockets[socket] and sum(map(len,pool.values()))<64:
            pool[socket].append(sockets[socket].pop(0))
coord_socket=min(pool)
coord=pool[coord_socket][0]
for b in threads:
    remaining={k:list(v) for k,v in pool.items()}
    remaining[coord_socket].remove(coord)
    groups=[[coord]]
    for _ in range(workers):
        eligible=[k for k,v in remaining.items() if len(v)>=b]
        assert eligible,(workers,b,remaining)
        socket=max(eligible,key=lambda k:len(remaining[k]))
        groups.append(remaining[socket][:b])
        remaining[socket]=remaining[socket][b:]
    assert len(set(x for g in groups for x in g))==1+workers*b<=64
    Path(f"cpu_groups_{b}.csv").write_text("".join(",".join(map(str,g))+"\n" for g in groups))
Path("cpu_groups.csv").write_text(Path(f"cpu_groups_{threads[0]}.csv").read_text())
Path("cpu_topology.csv").write_text(text)
Path("allocation_affinity.txt").write_text(str(sorted(allowed)))
