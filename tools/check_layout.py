"""Verify problem-first layout and preserved evidence; no simulation or file writes."""
import json
import argparse
from pathlib import Path
import re
import subprocess
import tomllib

ROOT = Path(__file__).resolve().parents[1]
BASE = "5742a5d10bf26149bf8f4800251da451efe39b82"
MOVES = (
    ("envs/spar_ddp/Project.toml", "problems/mpopf/spar_ddp/Project.toml"),
    ("envs/spar_ddp/README.md", "problems/mpopf/spar_ddp/ENVIRONMENT.md"),
    ("envs/msopf/Project.toml", "problems/msopf/Project.toml"),
    ("envs/msopf/README.md", "problems/msopf/ENVIRONMENT.md"),
    ("envs/multi_poi/Project.toml", "problems/msopf/archive/legacy-environment/Project.toml"),
    ("envs/msopf_sdp/", "problems/msopf/environments/sdp/"),
    ("envs/multi_poi/", "problems/msopf/"),
    ("envs/tadmm/parse_opendss.jl", "test_networks/parse_opendss.jl"),
    ("envs/tadmm/", "problems/mpopf/tadmm/"),
    ("spar_ddp/", "problems/mpopf/spar_ddp/"),
    ("rawData/", "test_networks/"),
    ("test_systems/", "test_networks/sources/"),
    ("results/", "problems/mpopf/tadmm/results/"),
)


def destination(name):
    for old, new in MOVES:
        if name == old or old.endswith("/") and name.startswith(old):
            return new + name[len(old):]
    return name


def git(*args, **kwargs):
    return subprocess.check_output(["git", *args], cwd=ROOT, **kwargs)


def main(strict_migration=False):
    try:
        tree = git("ls-tree", "-r", BASE, text=True)
    except subprocess.CalledProcessError:
        raise SystemExit("Evidence audit needs migration baseline history; fetch full Git history first.")
    expected = []
    for line in tree.splitlines():
        header, name = line.split("\t", 1)
        if not (name.startswith(("rawData/", "test_systems/", "results/", "figures/"))
                or any(x in name for x in ("/results/", "/logs/", "/patches/"))
                or name.startswith("spar_ddp/DDP4OPF.jl/src/")):
            continue
        if name in ("rawData/README.md", "figures/networks/README.md"):
            continue  # visitor navigation, not electrical or scientific evidence
        if name.startswith("spar_ddp/DDP4OPF.jl/src/") and not strict_migration:
            continue  # normal development may change maintained solver source
        target = destination(name)
        assert (ROOT / target).is_file(), f"Missing preserved artifact: {target}"
        expected.append((target, header.split()[2]))
    # Batch hashing uses the checkout's Git text filters, not platform-specific
    # CRLF bytes; binary payloads and evidence Git blobs must remain identical.
    payload = "".join(json.dumps(name, ensure_ascii=False)+"\n" for name, _ in expected)
    hashes = git("hash-object", "--stdin-paths", input=payload, text=True).splitlines()
    assert len(hashes) == len(expected)
    errors = [name for (name, old), actual in zip(expected, hashes) if old != actual]
    assert not errors, f"Changed preserved evidence/data/solver blobs: {errors[:20]}"
    ddp = ROOT / "problems/mpopf/spar_ddp"
    project = tomllib.loads((ddp / "Project.toml").read_text())
    assert (ddp / project["sources"]["DDP4OPF"]["path"] / "Project.toml").is_file()
    assert {"OpenDSSDirect", "Parameters"} <= project["deps"].keys()
    assert (ROOT / "problems/msopf/environments/sdp/Project.toml").is_file()
    # Check literal repository roots and includes without evaluating any driver.
    roots = includes = 0
    root_files = set()
    for base in (ROOT / "problems", ROOT / "test_networks"):
        for path in base.rglob("*.jl"):
            if not path.is_file():
                continue
            if any(x in path.parts for x in ("results", "experiments", "archive", "archived", "sources")):
                continue
            text = path.read_text(encoding="utf-8-sig")
            root_match = re.search(r'const (?:REPO|REPO_ROOT) = normpath\(joinpath\(@__DIR__, ((?:"\.\."(?:, )?)+)\)\)', text)
            if root_match:
                resolved = path.parent.joinpath(*re.findall(r'"([^"]+)"', root_match[1])).resolve()
                assert resolved == ROOT, f"Wrong repository root in {path.relative_to(ROOT)}: {resolved}"
                roots += 1
                root_files.add(path.relative_to(ROOT).as_posix())
            for match in re.finditer(r'include\(joinpath\((REPO|REPO_ROOT|@__DIR__), ((?:"[^"]+"(?:, )?)+)\)\)', text):
                parent = path.parent if match[1] == "@__DIR__" else ROOT
                target = parent.joinpath(*re.findall(r'"([^"]+)"', match[2])).resolve()
                assert target.is_file(), f"Missing literal include in {path.relative_to(ROOT)}: {target}"
                includes += 1
    assert {
        "problems/msopf/full_angle_pf.jl",
        "problems/msopf/root_level/convert_to_1ph.jl",
        "problems/mpopf/spar_ddp/examples/power_system/ieee123c_filterddp.jl",
        "problems/mpopf/spar_ddp/examples/power_system/centralized_ipopt_matched.jl",
        "problems/mpopf/spar_ddp/examples/power_system/opendss_check.jl",
    } <= root_files
    assert includes >= 20
    if strict_migration:
        old = git("show", f"{BASE}:envs/tadmm/parse_opendss.jl", text=True).partition("using OpenDSSDirect")[2]
        old = old.replace('        wd = @__DIR__\n        rawDataFolderPath = joinpath(wd, "..", "..", "rawData")',
                          '        rawDataFolderPath = @__DIR__')
        new = (ROOT / "test_networks/parse_opendss.jl").read_text().partition("using OpenDSSDirect")[2]
        assert old == new, "Shared parser numerical body changed"
    print(f"PASS: {len(expected)} preserved evidence/data/solver Git blobs; "
          f"{roots} repository roots; {includes} literal includes; colocated environments.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict-migration", action="store_true",
                        help="also prove unchanged active solver/parser source for the layout PR")
    main(parser.parse_args().strict_migration)
