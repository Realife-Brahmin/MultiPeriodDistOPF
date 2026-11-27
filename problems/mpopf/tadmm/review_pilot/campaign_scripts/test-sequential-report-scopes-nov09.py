"""Disposable fixtures: a snapshot pass must not become a full MPOPF pass."""
from pathlib import Path
import json
import shutil
import subprocess
import sys
import tempfile

root=Path(__file__).parent
p=root/'problems/mpopf/tadmm/review_pilot'
seq=p/'results/sequential_dss_nov09'
base=p/'results/persistent_caps_nov09/completed/results/run_T24_primal_dual_31665889'
work=Path(tempfile.mkdtemp(prefix='seq-scope-',dir=root/'scratch'))
out={}
for mode in ('legacy','chronological'):
    case=work/mode;shutil.copytree(base,case)
    if mode=='chronological':
        shutil.copy2(seq/'tables/central-model/summary.toml',case/'central-sequence.toml')
        shutil.copy2(seq/'tables/cap10-model/summary.toml',case/'endpoint-sequence.toml')
        m=case/'sweep.toml';m.write_text('sequential_validation_version = 1\ncentral_sequential_schedule_pass = true\ncentral_sequential_summary = "central-sequence.toml"\n'+m.read_text())
        a=case/'rho15000/assessment.toml';a.write_text('sequential_schedule_pass = false\nsequential_summary = "endpoint-sequence.toml"\n'+a.read_text())
    with (work/f'{mode}.log').open('w') as log:
        subprocess.run([sys.executable,str(p/'validate_horizon_sweep.py'),str(case),'complete'],check=True,stdout=log)
    v=json.loads((case/'validation_complete.json').read_text());arm=v['arms']['rho15000']
    assert v['central_opendss']['scope']=='network_snapshot_only'
    assert arm['full_mpopf_validation_pass']==(False if mode=='chronological' else None)
    if mode=='chronological':
        assert v['central_sequential']['sequential_schedule_pass'] is True
        assert arm['sequential']['sequential_schedule_pass'] is False
    out[mode]={'snapshot_converged':v['central_opendss']['converged'],'full_mpopf_validation_pass':arm['full_mpopf_validation_pass'],'test_passed':True}
(seq/'report-scope-tests.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
