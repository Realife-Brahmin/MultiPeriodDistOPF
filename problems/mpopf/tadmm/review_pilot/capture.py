"""Freeze source before launching a case; later edits cannot alter provenance."""
from pathlib import Path
import hashlib
import shutil
import subprocess
import sys

pilot = Path(__file__).resolve().parent
root = pilot.parents[3]
out = Path(sys.argv[1]).resolve()
system = sys.argv[2]
files = [p for p in pilot.iterdir() if p.is_file() and p.suffix in ('.jl', '.py', '.ps1', '.sh', '.toml')]
files += [p for p in (pilot/'slurm').glob('*.toml') if p.is_file()]
files += [p for p in (pilot/'slurm/bin').glob('*') if p.is_file()]
files += sorted((root / 'test_networks' / system).glob('*.dss'))
files += [root / 'test_networks' / 'parse_opendss.jl']
shared = root / 'problems/mpopf/spar_ddp/examples/power_system'
files += [shared / f for f in ('export_ieee123c_data.jl', 'loadflow_start.jl', 'terminal_soc_penalty.jl', 'control_layout.jl', 'opendss_check.jl', 'voltage_screening.jl')]
rows = []
for p in files:
    rel = p.relative_to(root)
    dest = out / 'source' / rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(p, dest)
    rows.append(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {rel.as_posix()}\n')
(out / 'source.sha256').write_text(''.join(rows))
revision = (root / 'source.commit').read_text() if (root / 'source.commit').is_file() else subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True)
(out / 'base_commit.txt').write_text(revision)
