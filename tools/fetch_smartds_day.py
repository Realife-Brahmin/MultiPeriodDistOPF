"""One day of the SMART-DS P1U customer load profiles (user, 2026-11-03).

Every profile the P1U deck references (LoadShapes.dss: kW and kvar multipliers, 15-minute, one year from
1 January) is read from the public OEDI bucket only up to the end of the requested day (HTTP range), and the
day's 96 values are written to test_networks/sources/smartds_sfo_p1u/profiles_day<DAY>.csv, one row per file.
Default day 148 (28 May 2018), the dataset's own peak day (analysis/Summary_data.csv).

    python -I tools/fetch_smartds_day.py [day]
"""
import concurrent.futures as cf
import os
import re
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "test_networks", "sources", "smartds_sfo_p1u")
URL = "https://oedi-data-lake.s3.amazonaws.com/SMART-DS/v1.0/2018/SFO/P1U/profiles/"
DAY = int(sys.argv[1]) if len(sys.argv) > 1 else 148
FIRST, LAST = (DAY - 1) * 96, DAY * 96          # rows [FIRST, LAST) of the 35040


def day_values(name):
    nbytes = 25 * LAST                            # lines are under 25 bytes; read more if not
    while True:
        req = urllib.request.Request(URL + name, headers={"Range": f"bytes=0-{nbytes - 1}"})
        with urllib.request.urlopen(req, timeout=120) as r:
            text = r.read().decode()
        lines = text.split("\n")
        if len(lines) > LAST:                     # the last line may be cut; rows below LAST are whole
            return [float(x) for x in lines[FIRST:LAST]]
        nbytes *= 2


def main():
    shapes = open(os.path.join(SRC, "LoadShapes.dss")).read()
    names = sorted(set(re.findall(r"file=(?:\.\./)+profiles/([^)\s]+)", shapes)))
    out = os.path.join(SRC, f"profiles_day{DAY}.csv")
    rows = {}
    with cf.ThreadPoolExecutor(16) as ex:
        for k, (n, v) in enumerate(zip(names, ex.map(day_values, names)), 1):
            rows[n] = v
            if k % 500 == 0:
                print(f"  {k} of {len(names)}", flush=True)
    with open(out, "w") as f:
        f.write("profile," + ",".join(f"q{i}" for i in range(96)) + "\n")
        for n in names:
            f.write(n + "," + ",".join(repr(x) for x in rows[n]) + "\n")
    print(f"wrote {out}: {len(names)} profiles, day {DAY}")


if __name__ == "__main__":
    main()
