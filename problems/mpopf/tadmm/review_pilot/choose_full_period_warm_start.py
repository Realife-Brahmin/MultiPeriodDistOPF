"""Choose a validated implementation at the existing quality goal, then run T96.

This is implementation selection, not rho tuning. All tested outcomes and
selection cost are retained. No submission occurs without explicit --submit.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
try:
    import tomllib
except ImportError:
    import tomli as tomllib

C=Path(sys.argv[1]);jobs=[int(j) for j in sys.argv[2:] if j!='--submit'];records=[];eligible=[]
for job in jobs:
    matches=list((C/'results').glob(f'run_T48_*_{job}'))
    record=dict(job=job,eligible=False)
    if len(matches)==1:
        case=matches[0];record['case']=str(case)
        if (case/'SWEEP_VALIDATED').is_file():
            m=tomllib.loads((case/'sweep.toml').read_text())
            a=tomllib.loads((case/'rho120000/assessment.toml').read_text())
            assert m['T']==48 and m['workers']==48 and m['linear_solver']=='ma57'
            assert m['execution_policy']=='one_per_period' and m['objective_update']=='batched'
            assert m['terminal_soc_penalty']==0 and m['warm_start'] in ('primal','primal_dual')
            record.update(warm_start=m['warm_start'],core_s=a['core_s'],iterations=a['iterations'],
                          objective_gap=a['objective_gap'],original_primal=a['original_primal'],
                          final_rms_met=a['final_rms_met'],input_sha256=m['input_sha256'],C_B=m['C_B'],
                          central_initial_starts_sha256=m['central_initial_starts_sha256'],
                          validation_sha256=hashlib.sha256((case/'validation_complete.json').read_bytes()).hexdigest())
            # This is the declared residual/dual + cost goal. Preserve original
            # feasibility separately; do not label it strict physical NO.
            record['eligible']=a['final_rms_met'] and a['objective_gap']<=.001
            if record['eligible']:eligible.append(record)
    records.append(record)
decision=dict(utc=datetime.now(timezone.utc).isoformat(),candidates=records,
    selection_rule='Lowest measured core time among completed validated T48 arms at final owned/all-copy RMS<=1e-5, unique dual<=1e-2 and objective gap<=0.1%. Original physical error is recorded separately, not relabeled strict NO.',
    tuning_cost_note='Both complete warm-start implementation comparisons, all warmups and all central references count toward exploration cost.',
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
dest=C/'campaign/warm-start-selection.json'
if not eligible:
    decision['status']='No validated candidate at the declared quality goal; no T96 job submitted.'
    dest.write_text(json.dumps(decision,indent=2)+'\n')
    raise SystemExit(decision['status'])
identities={(r['input_sha256'],r['C_B'],r['central_initial_starts_sha256']) for r in records if 'input_sha256' in r}
assert len(identities)==1
winner=min(eligible,key=lambda r:r['core_s']);decision['selected_warm_start']=winner['warm_start']
cmd=['sbatch','--parsable','--account=dubey','--partition=kamiak','--nodes=2','--ntasks=2',
     '--cpus-per-task=64','--hint=nomultithread','--mem=128G','--time=08:00:00',
     '--constraint=sapphirerapids&gold-6438y+-2.00ghz&mem-1024GB&cores-64',
     '--job-name=tadmm-T96-96workers',f'--output={C}/logs/T96-96workers-%j.out',
     f'--error={C}/logs/T96-96workers-%j.err',str(C/'campaign/after-gate.sh'),str(C),'96',
     winner['warm_start'],'60000,240000']
decision['submission_command']=cmd
if '--submit' in sys.argv:
    assert not dest.exists(), 'Refuse duplicate automatic T96 submission'
    # Preserve the decision before submitting; retain raw scheduler response.
    dest.write_text(json.dumps(decision,indent=2)+'\n')
    decision['T96_job']=subprocess.check_output(cmd,text=True).strip()
dest.write_text(json.dumps(decision,indent=2)+'\n');print(json.dumps(decision,indent=2))
