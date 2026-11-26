"""Dependency-free navigation/syntax checks, not numerical validation or benchmarks."""
import ast
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
# Deliberately curated: frozen evidence and old notes retain historical paths.
NAVIGATION = (
    "README.md", "problems/mpopf/spar_ddp/README.md", "problems/mpopf/spar_ddp/experiments/README.md",
    "problems/mpopf/spar_ddp/experiments/value_coupling_audit/README.md",
    "problems/mpopf/spar_ddp/examples/power_system/README.md", "problems/mpopf/spar_ddp/results/README.md",
    "problems/mpopf/spar_ddp/notes/README.md", "scratch/README.md",
    "scratch/archive/workstation-queues/README.md", "tools/README.md",
    "docs/repository-maintenance.md", "test_networks/README.md", "problems/README.md",
    "AGENTS.md", ".claude/CLAUDE.md", "problems/mpopf/spar_ddp/CLAUDE.md",
    "docs/problem-layout-migration.md",
    "problems/msopf/README.md", "problems/msopf/ENVIRONMENT.md",
    "problems/mpopf/tadmm/README.md", "problems/mpopf/spar_ddp/ENVIRONMENT.md",
    "docs/first-time-installation-guide.md",
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
    for folder in ("tools", "problems"):
        for path in sorted((root / folder).rglob("*.py")):
            if "__pycache__" in path.parts:
                continue
            try:
                ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
                count += 1
            except (SyntaxError, UnicodeError) as exc:
                errors.append(f"{path.relative_to(root)}: {exc}")
    for name in ("problems/mpopf/spar_ddp/SPARDDP.jl/Project.toml", "problems/mpopf/spar_ddp/Project.toml",
                 "problems/msopf/Project.toml", "problems/mpopf/tadmm/Project.toml"):
        if not (root / name).is_file():
            errors.append(f"missing workstream environment/package: {name}")
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"PASS: {len(NAVIGATION)} navigation pages; {count} Python files parsed. "
          "No numerical tests or simulations executed.")


if __name__ == "__main__":
    check()
