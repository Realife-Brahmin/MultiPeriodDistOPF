import re,sys
from pathlib import Path
root=Path(sys.argv[1])
def medians(folder):
    text=(root/folder/"replay.log").read_text()
    assert "BALANCED_COMPLETE" in text
    return {m:float(t) for m,t in re.findall(r"BALANCED_RESULT scheduling=(\S+) workers=\d+ median_s=(\S+)",text)}
base=medians("kernel63")
options=[(base["current"],63,1)]
for label,workers in (("threads31",31),("threads15",15)):
    for mode,time in medians(label).items():
        options.append((time,workers,int(mode.split("_")[1])))
best=min(options)
if base["right"] < base["current"] and best[0]<base["current"] and best[1]!=63:
    print(best[1],best[2])
else:
    print(0,0)
