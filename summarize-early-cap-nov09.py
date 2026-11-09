"""Compare a validated early-inexact trajectory with the preserved exact arm."""
import argparse
import csv
import json
from pathlib import Path
import tomllib

def rows(p):
    with p.open() as f:return list(csv.DictReader(f))
parser=argparse.ArgumentParser()
parser.add_argument('case',type=Path);parser.add_argument('audit',type=Path)
parser.add_argument('baseline',type=Path);parser.add_argument('output',type=Path)
args=parser.parse_args()
assert (args.case/'SWEEP_VALIDATED').exists()
m=tomllib.loads((args.case/'sweep.toml').read_text())
bm=tomllib.loads((args.baseline/'sweep.toml').read_text())
for k in ('system','T','C_B','input_sha256','central_initial_starts_sha256','rhos','rho_policy','rho_floor','workers','linear_solver','warm_start','objective_update'):
    assert m[k]==bm[k],(k,m[k],bm[k])
arm=args.case/'rho15000';baseline=args.baseline/'rho15000'
a=tomllib.loads((arm/'assessment.toml').read_text());ba=tomllib.loads((baseline/'assessment.toml').read_text())
audit=tomllib.loads(args.audit.read_text());e=audit['rho15000/endpoint.jls']
details=rows(arm/'trace_subproblems.csv'); trace=rows(arm/'trace.csv')
accepted=[r for r in details if r['inexact_accepted']=='true']
recovered=[r for r in details if r['cap_recovery']=='true']
attempts=[r for r in details if int(r['cap_attempt_iterations'])>0]
late=[r for r in details if int(r['iteration'])>m['cap_end_iteration']]
assert late and all(r['status'] in ('LOCALLY_SOLVED','OPTIMAL') for r in late)
checks=json.loads((args.case/'validation_complete.json').read_text())
result=dict(case=str(args.case),source=(args.case/'base_commit.txt').read_text().strip(),
    cap_screen=m.get('cap_screen','feasible'),
    cap_iterations=m['cap_iterations'],cap_end_iteration=m['cap_end_iteration'],cap_local_feasibility_tolerance=m['cap_primal_tolerance'],
    iterations=a['iterations'],baseline_iterations=ba['iterations'],core_s=a['core_s'],baseline_core_s=ba['core_s'],
    loop_s=a['loop_wall_s'],baseline_loop_s=ba['loop_wall_s'],central_s=m['central_wall_s'],baseline_central_s=bm['central_wall_s'],
    relative_core_change=a['core_s']/ba['core_s']-1,final_rms_met=a['final_rms_met'],
    strict_physical_met=a['strict_criteria_met'],cost_gap=a['objective_gap'],
    max_integrated_energy_bound_kwh=e['integrated_energy_bound_kwh'],max_SOC_bound_percentage_points=e['max_SOC_bound_percentage_points'],
    max_dynamics_error_kwh=e['max_dynamics_error_kwh'],accepted_inexact=len(accepted),recovered_exact=len(recovered),
    cap_attempts=[{k:r[k] for k in ('iteration','period','cap_attempt_iterations','cap_attempt_primal','inexact_accepted','cap_recovery','ipopt_iterations','solver_s','local_primal_gate_s')} for r in attempts],
    accepted_gate_max=max((float(r['local_primal']) for r in accepted),default=None),
    final_round_all_exact=all(r['status'] in ('LOCALLY_SOLVED','OPTIMAL') for r in details if int(r['iteration'])==a['iterations']),
    first40_core_s=float(trace[39]['algorithm_s']) if len(trace)>=40 else None,
    validation=checks,
    note='Single measured trajectory after full warmup/reset. Core counts local screening and recovery. Baseline SOC-bound error0.869624kWh; strict physical feasibility was false. No physical-quality win is inferred from RMS stopping alone.')
args.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('validation','cap_attempts')},indent=2))
print('CAP_ATTEMPTS',json.dumps(result['cap_attempts']))
