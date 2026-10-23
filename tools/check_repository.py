"""Dependency-free navigation/syntax checks, not numerical validation or benchmarks."""
import ast
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
# Deliberately curated: frozen evidence and old notes retain historical paths.
NAVIGATION = (
    "README.md", "spar_ddp/README.md", "spar_ddp/experiments/README.md",
    "spar_ddp/experiments/value_coupling_audit/README.md",
    "spar_ddp/examples/power_system/README.md", "spar_ddp/results/README.md",
    "spar_ddp/notes/README.md", "scratch/README.md",
    "scratch/archive/workstation-queues/README.md", "tools/README.md",
    "docs/repository-maintenance.md",
)


def check(root=ROOT):
    errors = []
    count = 0
    for name in NAVIGATION:
        page = root / name
        if not page.is_file():
            errors.append(f"missing navigation page: {name}")
            continue
        for raw in re.findall(r"\[[^\]\n]*\]\(([^)]+)\)", page.read_text(encoding="utf-8")):
            target = urlsplit(raw.strip("<>"))
            if target.scheme or target.netloc or not target.path:
                continue
            path = (page.parent / unquote(target.path)).resolve()
            if not path.is_relative_to(root.resolve()) or not path.exists():
                errors.append(f"{name}: missing/outside-repository link {raw}")
    for folder in ("tools", "spar_ddp"):
        for path in sorted((root / folder).rglob("*.py")):
            if "__pycache__" in path.parts:
                continue
            try:
                ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
                count += 1
            except (SyntaxError, UnicodeError) as exc:
                errors.append(f"{path.relative_to(root)}: {exc}")
    for name in ("spar_ddp/DDP4OPF.jl/Project.toml", "envs/spar_ddp/Project.toml",
                 "envs/msopf/Project.toml", "envs/tadmm/Project.toml"):
        if not (root / name).is_file():
            errors.append(f"missing workstream environment/package: {name}")
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"PASS: {len(NAVIGATION)} navigation pages; {count} Python files parsed. "
          "No numerical tests or simulations executed.")


if __name__ == "__main__":
    check()
