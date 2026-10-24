# SPAR-DDP: context for Claude

Claude Code loads this file for work under `problems/mpopf/spar_ddp/`. The current SPAR-DDP instructions and
evidence are in the "Current SPAR-DDP evidence" section of [.agents/AGENTS.md](../../../.agents/AGENTS.md);
read it first and keep it as the single place they are maintained. Repository-wide rules are in
[.claude/CLAUDE.md](../../../.claude/CLAUDE.md). The dated standing rules from before the October 18 cleanup
are in [docs/history/CLAUDE-pre-cleanup-oct18.md](../../../docs/history/CLAUDE-pre-cleanup-oct18.md):
consult them when relevant, as history, not as current authorization.

SPAR-DDP shares only the test decks and the OpenDSS checks with MSOPF (`problems/msopf/`); do not
apply MSOPF's modelling choices (hard terminal SOC, the full-angle model) here.
