"""Export one adapted OpenDSS deck with shared profiles; never change input files."""
import argparse
from pathlib import Path
import re
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PROFILES = ("LoadShapeDefault.dss", "LoadShapePVDefault.dss", "LoadShapePSubsCostDefault.dss")


def package_network(system, output, root=ROOT):
    root = Path(root).resolve()
    data = root / "rawData"
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", system):
        raise ValueError("Use a network folder name, not a path")
    deck = data / system
    if deck.is_symlink() or not (deck / "Master.dss").is_file():
        raise ValueError(f"No adapted network Master.dss for {system}")
    output = Path(output).resolve()
    if output.is_relative_to(data.resolve()):
        raise ValueError("Write the ZIP outside rawData to keep the input tree unchanged")
    entries = {}
    for path in sorted(deck.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Refusing symbolic link in network: {path}")
        if path.is_file():
            entries[(Path("test_networks") / system / path.relative_to(deck)).as_posix()] = path
    for profile in PROFILES:
        entries[f"test_networks/{profile}"] = data / profile
    for filename in ("LICENSE", "CITATION.cff"):
        entries[filename] = root / filename
    # Read first: missing dependencies must not leave a partial export.
    payloads = {name: path.read_bytes() for name, path in entries.items()}
    payloads["README.txt"] = (
        f"OpenDSS test network: {system}\n\n"
        f"Compile test_networks/{system}/Master.dss, then Solve in OpenDSS.\n"
        "Keep the folder structure: the network redirects shared profiles from its parent.\n"
        "This is a power-flow deck, not an optimized battery schedule.\n\n"
        "Please cite Aryan Ritwajeet Jha and contributors, MultiPeriodDistOPF,\n"
        "https://github.com/Realife-Brahmin/MultiPeriodDistOPF.\n"
        "See CITATION.cff and the upstream notices retained in the deck.\n"
        "Network provenance: https://github.com/Realife-Brahmin/MultiPeriodDistOPF/blob/master/test_networks/README.md\n"
        "The citation request does not add restrictions to the MIT license.\n"
    ).encode("utf-8")
    # Exclusive creation and fixed metadata: no overwrite, reproducible bytes.
    with zipfile.ZipFile(output, "x", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, content in sorted(payloads.items()):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, content)
    return len(payloads)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("system", help="raw network folder, e.g. ieee123C_1ph")
    parser.add_argument("--output", required=True, type=Path, help="new output ZIP path")
    args = parser.parse_args()
    try:
        count = package_network(args.system, args.output)
    except (ValueError, OSError) as exc:
        parser.exit(1, f"Cannot package network: {exc}\n")
    print(f"Created {args.output} ({count} files); source inputs unchanged")


if __name__ == "__main__":
    main()
