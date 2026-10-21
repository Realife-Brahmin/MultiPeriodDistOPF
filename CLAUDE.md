# Claude Code repository guidance

This repository holds separate projects that share little beyond the test decks and the
OpenDSS validation. Each keeps its own Claude notes, loaded when you work in its folder:

- **MSOPF** (multi-source multi-period OPF, `envs/multi_poi/`): [envs/multi_poi/CLAUDE.md](envs/multi_poi/CLAUDE.md).
- **SPAR-DDP** (`spar_ddp/`, `envs/spar_ddp/`): [spar_ddp/CLAUDE.md](spar_ddp/CLAUDE.md), which points to
  the current DDP instructions in [AGENTS.md](AGENTS.md).
- **tADMM** (`envs/tadmm/`): kept for reference; see [envs/tadmm/README.md](envs/tadmm/README.md).

Read [AGENTS.md](AGENTS.md) for the shared instructions. Historical Claude context is preserved
in [docs/history/CLAUDE-pre-cleanup-oct18.md](docs/history/CLAUDE-pre-cleanup-oct18.md); old
launch requests there are not current authorization.

## Shared by all projects

- **Test decks** are the balanced single-phase OpenDSS decks in `rawData/` (raw public feeders in
  `test_systems/`, fetched by `bash tools/fetch_test_systems.sh`). Rules differ per project (for
  example terminal SOC is soft in SPAR-DDP and tADMM, hard in MSOPF); never carry one project's
  modelling rules into another.
- **OpenDSS validation:** check a dispatch against OpenDSS before trusting an instance or a result,
  and run each solve's OpenDSS check as part of the run.
- **Branches and commits:** each project works on its own dated branch and merges through a PR.
  Make focused commits that tell the story; check after every push that it landed. Logs over
  100 MB cannot be pushed: commit them gzipped.
- **Julia** through Bash, not PowerShell (PowerShell quoting has corrupted Julia commands). Pkg
  strips `Project.toml` comments, so put rationale in a README. Timed launchers unset `NO_COLOR`
  and `FORCE_COLOR`. Stop Julia by PID.
- **Kamiak** (partition `vcea`, account `dubey`) is shared between projects: check the queue
  before submitting, and do not take resources another project's jobs need.
