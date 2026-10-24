# Research scratch work

Exploratory helpers, old experiment queues, diagnostic probes and ad hoc SVG
editing tools live here. These are not the supported starting point for users.
For SPAR-DDP use ../spar_ddp/README.md; for tADMM use ../envs/tadmm/README.md.

- `probes/`: one-off diagnostics and research experiments.
- `archive/workstation-queues/`: superseded machine-specific unattended queues
  and their monitor, preserved as history and not runnable entry points.
- Remaining `run_*.ps1`, `run_*.sh`: individual experimental/regression helpers;
  these still require explicit review before use.
- `extract_*.ps1`, `analyze_*.py`, `_export_traj_csv.py`: analysis/conversion.
- `make_network_svg*.py`, `strip_svg_captions.py`, `tint_svg_elements.py`:
  exploratory network illustration tools; committed Inkscape SVGs are the
  editable figure sources, not promised exact outputs of these generators.

Inspect inputs, output paths and side effects before running anything. Some
queues launch expensive jobs or commit/push results. Nothing here is an
automatic instruction to resume old experiments.

Review this directory after each research milestone. Promote maintained tools
to their owning workstream or tools/, archive scripts needed to explain saved
evidence, and remove obsolete scratch only after checking references. The
resource fetcher is maintained separately at ../tools/fetch_resources.sh.
