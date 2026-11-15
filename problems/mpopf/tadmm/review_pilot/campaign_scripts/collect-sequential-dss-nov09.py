"""Package complete discrepancy tables separately from compact audit metadata."""
from pathlib import Path
import hashlib
import json
import shutil
import tarfile
import tomli

c=Path('/home/aryan.jha/kamiak_jobs/tadmm_review_nov09_sequential_dss_v2')
dest=c/'campaign/validation_source';dest.mkdir(exist_ok=True)
for name in ('sequential_dss.jl','validate_sequential_dss.jl','run_sequential_dss.jl'):
    shutil.copy2(c/'source/problems/mpopf/tadmm/review_pilot'/name,dest/name)
shutil.copy2(c/'source/problems/mpopf/spar_ddp/examples/power_system/control_layout.jl',dest/'control_layout.jl')
shutil.copy2(c/'source/test_networks/large10kC_1ph/Storage.dss',dest/'Storage.dss')
out=c.parent;summary=[];hashes={}
for case in sorted((c/'results').iterdir()):
    if not case.is_dir():continue
    s=tomli.loads((case/'summary.toml').read_text());s['case']=case.name
    summary.append(s)
    archive=out/f'sequential-dss-{case.name}-nov09.tgz'
    with tarfile.open(archive,'w:gz') as tar:tar.add(case,arcname=case.name)
    hashes[archive.name]=hashlib.sha256(archive.read_bytes()).hexdigest()
    print(case.name,{k:s[k] for k in ('sequential_schedule_pass','dispatch_energy_execution_pass','opf_state_agreement_pass','max_power_tracking_kw','max_engine_energy_balance_kwh','max_engine_SOC_bound_kwh','max_requested_SOC_bound_kwh','max_engine_vs_opf_energy_kwh','reserve_or_capacity_clamp_records')})
(c/'campaign/replay-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
archive=out/'sequential-dss-metadata-nov09.tgz'
with tarfile.open(archive,'w:gz') as tar:
    for directory in ('campaign','logs'):
        for p in (c/directory).rglob('*'):
            if p.is_file() and p.suffix not in ('.gz','.tgz'):tar.add(p,arcname=str(p.relative_to(c)))
    for name in ('execution-source.sha256','SEQUENTIAL_GATE_VALIDATED','SEQUENTIAL_REPLAYS_RECORDED'):
        tar.add(c/name,arcname=name)
hashes[archive.name]=hashlib.sha256(archive.read_bytes()).hexdigest()
(out/'sequential-dss-archives-nov09.json').write_text(json.dumps(hashes,indent=2)+'\n')
print(json.dumps(hashes,indent=2))
