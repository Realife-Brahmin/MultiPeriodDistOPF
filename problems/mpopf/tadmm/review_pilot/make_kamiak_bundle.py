"""Freeze committed numerical code and only the inputs needed by this pilot."""
import gzip
import io
from pathlib import Path
import subprocess
import sys
import tarfile

root = Path(__file__).resolve().parents[4]
out = Path(sys.argv[1]).resolve()
assert not out.exists(), out
assert not subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=root, text=True).strip(), 'Commit tracked changes first'
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root).strip() + b'\n'
base = 'problems/mpopf/spar_ddp/examples/power_system/'
paths = ['AGENTS.md', '.agents/AGENTS.md', 'problems/mpopf/tadmm/review_pilot',
    'test_networks/parse_opendss.jl', 'test_networks/LoadShapeDefault.dss',
    'test_networks/LoadShapePVDefault.dss', 'test_networks/LoadShapePSubsCostDefault.dss']
paths += [base + f for f in ('export_ieee123c_data.jl', 'loadflow_start.jl',
    'terminal_soc_penalty.jl', 'control_layout.jl', 'opendss_check.jl', 'voltage_screening.jl')]
paths += ['test_networks/' + s for s in ('ieee123C_1ph', 'ieee2522C_1ph', 'large10kC_1ph')]
tracked = subprocess.check_output(['git', 'ls-files', '-z', '--', *paths], cwd=root).decode().split('\0')
tracked = [p for p in tracked if p and '/results/' not in p]
raw = subprocess.check_output(['git', '-c', 'core.autocrlf=false', 'archive', '--format=tar', 'HEAD', *tracked], cwd=root)
with tarfile.open(fileobj=io.BytesIO(raw)) as check:
    for member in check.getmembers():
        if member.isfile() and member.name.endswith('.sh'):
            assert b'\r' not in check.extractfile(member).read(), member.name
        if member.name.endswith('/worker_julia.sh'):
            assert member.mode & 0o111, 'worker_julia.sh must be executable'
buffer = io.BytesIO(raw)
with tarfile.open(fileobj=buffer, mode='a') as archive:
    info = tarfile.TarInfo('source.commit')
    info.size = len(revision)
    archive.addfile(info, io.BytesIO(revision))
out.write_bytes(gzip.compress(buffer.getvalue(), mtime=0))
print('BUNDLE', out, 'commit', revision.decode().strip(), 'bytes', out.stat().st_size)
