# Claude Code repository guidance

This repository holds separate projects that share little beyond the test decks and the
OpenDSS validation. Each keeps its own Claude notes, loaded when you work in its folder:

- **MSOPF** (multi-source multi-period OPF, `problems/msopf/`): [envs/multi_poi/CLAUDE.md](../problems/msopf/CLAUDE.md).
- **SPAR-DDP** (`problems/mpopf/spar_ddp/`): [project notes](../problems/mpopf/spar_ddp/CLAUDE.md), which points to
  the current DDP instructions in [.agents/AGENTS.md](../.agents/AGENTS.md).
- **tADMM** (`problems/mpopf/tadmm/`): kept for reference; see [envs/tadmm/README.md](../problems/mpopf/tadmm/README.md).

Shared instructions are imported here so this hidden entry point loads them:

@../.agents/AGENTS.md

Historical Claude context is preserved
in [docs/history/CLAUDE-pre-cleanup-oct18.md](../docs/history/CLAUDE-pre-cleanup-oct18.md); old
launch requests there are not current authorization.

## Shared by all projects

- **Test decks** are the balanced single-phase OpenDSS decks in `test_networks/` (raw public feeders in
  `test_networks/sources/`, fetched by `bash tools/fetch_test_systems.sh`). Rules differ per project (for
  example terminal SOC is soft in SPAR-DDP and tADMM, hard in MSOPF); never carry one project's
  modelling rules into another.
- **OpenDSS validation:** check a dispatch against OpenDSS before trusting an instance or a result,
  and run each solve's OpenDSS check as part of the run.
- **Branches and commits:** each project works on its own dated branch and merges through a PR.
  Make focused commits that tell the story; check after every push that it landed. No
  `Co-Authored-By: Claude` trailer or other Claude attribution in commit messages (user, 2026-10-21). Logs over
  100 MB cannot be pushed: commit them gzipped.
- **Julia** through Bash, not PowerShell (PowerShell quoting has corrupted Julia commands). Pkg
  strips `Project.toml` comments, so put rationale in a README. Timed launchers unset `NO_COLOR`
  and `FORCE_COLOR`. Stop Julia by PID.
- **Kamiak** (partition `vcea`, account `dubey`) is shared between projects: check the queue
  before submitting, and do not take resources another project's jobs need.
