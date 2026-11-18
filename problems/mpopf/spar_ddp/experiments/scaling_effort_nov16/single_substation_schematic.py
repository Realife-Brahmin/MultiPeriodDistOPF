"""Schematics of the single-substation med3980 and sf6937 (SPAR-DDP), drawn by the MSOPF generator
(problems/msopf/root_level/network_schematic.py, unchanged) with configs for these decks:
substation 1 only, the ties closed in test_networks/<deck>_sub1_1ph/Master.dss (dashed), and the
other substations' transformers (disabled there) left out. PV squares and battery rings as in the
PESGM maps; bus positions from the source models' coordinates.

    python single_substation_schematic.py <med3980|sf6937> <out.svg>
    inkscape <out.svg> --export-area-drawing --export-type=pdf --export-filename=<out.pdf>
"""
import re, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "problems" / "msopf" / "root_level"))
import network_schematic as ns

VARIANT = {"med3980": ("ieee9500xf_1ph", "ieee9500xf_sub1_1ph"), "sf6937": ("smartds_p1u_xf_1ph", "smartds_p1u_xf_sub1_1ph")}
# The SPAR-DDP paper's markers (TPWRS-2027-Sparse-Parallel-DDP-for-MPOPF/figures/networks/README.md):
# amber PV squares and lavender battery rings, not the MSOPF paper's green batteries.
ns.BAT_SOFT, ns.BAT_FILL_BIG = "#7956A6", "#EEE7F6"
ns.PV_SOFT, ns.PV_LINE_SOFT = "#C99746", "#8C612F"

def main(system, out):
    src, var = VARIANT[system]
    master = (REPO / "test_networks" / var / "Master.dss").read_text()
    closed = re.findall(r"(?im)^line\.(\S+)\.enabled=true", master)
    dropped = {n.lower() for n in re.findall(r"(?im)^line\.(\S+)\.enabled=false", master)}
    cfg = dict(ns.CONFIG[system], deck=src, variant=None, closed=closed, subs=["subs1"], added={})
    lines, pvs, bats, subs = ns.deck(cfg)
    lines = [l for l in lines if l[0] not in dropped]
    xy, _, _ = ns.coords(system)
    pos, _ = ns.place(lines, xy)
    print(f"{system}: {len(lines)} branches ({sum(t for *_, t in lines)} closed ties), {len(pvs)} PV, {len(bats)} batteries, "
          f"substation {subs}; {len(dropped)} transformers left out")
    ns.map_main(system, out, cfg, lines, pvs, bats, subs, pos)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
