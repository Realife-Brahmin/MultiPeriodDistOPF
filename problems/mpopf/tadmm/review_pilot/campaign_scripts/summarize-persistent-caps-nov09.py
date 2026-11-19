"""Quality/timing of persistent caps; no implicit exact-final requirement."""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics
import tomllib

def rows(p):
    with p.open() as f:return list(csv.DictReader(f))
def dist(x):
    x=sorted(x)
    return dict(median=statistics.median(x),p95=x[math.ceil(.95*len(x))-1],maximum=max(x))
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('baseline',type=Path);p.add_argument('output',type=Path)
args=p.parse_args()
bm=tomllib.loads((args.baseline/'sweep.toml').read_text())
ba=tomllib.loads((args.baseline/'rho15000/assessment.toml').read_text())
result=[]
for case in sorted((args.root/'results').glob('run_*')):
    assert (case/'SWEEP_VALIDATED').exists()
    m=tomllib.loads((case/'sweep.toml').read_text());arm=case/'rho15000'
    for key in ('system','T','C_B','input_sha256','central_initial_starts_sha256','rhos','rho_policy','rho_floor','workers','linear_solver','warm_start','objective_update'):
        assert m[key]==bm[key],key
    a=tomllib.loads((arm/'assessment.toml').read_text());trace=rows(arm/'trace.csv');sub=rows(arm/'trace_subproblems.csv')
    assert m['local_solve_policy']=='persistent_cap'
    assert all(int(r['ipopt_iterations'])<=m['cap_iterations'] and r['cap_recovery']=='false' and r['cold_retry']=='false' for r in sub)
    audit=tomllib.loads((args.root/'campaign'/f"T24-{case.name.split('_')[-1]}-energy.toml").read_text())
    e=audit['rho15000/endpoint.jls'];v=json.loads((case/'validation_complete.json').read_text())
    checks=v['arms']['rho15000']['opendss'];endpoint=checks['endpoint.jls']
    last=sub[-m['T']:];tail=trace[-min(50,len(trace)):]
    result.append(dict(case=case.name,source=(case/'base_commit.txt').read_text().strip(),
        cap=m['cap_iterations'],iterations=a['iterations'],consensus_target_met=a['final_rms_met'],
        core_s=a['core_s'],loop_s=a['loop_wall_s'],central_s=m['central_wall_s'],baseline_core_s=ba['core_s'],
        relative_core_change=a['core_s']/ba['core_s']-1,
        objective_gap=a['objective_gap'],signed_objective_gap=(a['objective']-m['central_objective'])/abs(m['central_objective']),
        original_max_violation=a['original_primal'],physical_NO_met=a['original_primal']<=1e-5 and a['objective_gap']<=.001,
        strict_driver_flag=a['strict_criteria_met'],owned_rms=a['owned_r_rms'],all_copy_rms=float(trace[-1]['r_rms']),dual_rms=a['owned_s_rms'],
        max_integrated_energy_bound_kwh=e['integrated_energy_bound_kwh'],max_SOC_bound_percentage_points=e['max_SOC_bound_percentage_points'],
        max_dynamics_error_kwh=e['max_dynamics_error_kwh'],max_power_bound_kw=e['max_power_bound_kw'],normalized_cone_gap=e['normalized_cone_gap'],
        truncated_local_solves=sum(r['inexact_accepted']=='true' for r in sub),total_local_solves=len(sub),
        final_round_truncated=sum(r['inexact_accepted']=='true' for r in last),
        local_ipopt_iterations=dist([int(r['ipopt_iterations']) for r in sub]),
        local_solver_s=dist([float(r['solver_s']) for r in sub]),
        tail_last_n=len(tail),tail_all_copy_rms_min=min(float(r['r_rms']) for r in tail),tail_all_copy_rms_max=max(float(r['r_rms']) for r in tail),
        endpoint_opendss=endpoint,all_saved_network_checks_converged=all(c['converged'] for c in checks.values()),
        all_saved_network_checks_zero_voltage_violations=all(c['zero_voltage_violations'] for c in checks.values()),
        network_records=len(checks)+1,
        max_saved_network_voltage_difference_pu=max(c['max_voltage_difference_pu'] for c in [v['central_opendss'],*checks.values()]),
        timing_note='Core excludes all original-model/SOC/network diagnostics; finite-only persistent caps have no timed local feasibility report. Includes warm starts, updates, communication and coordination. Prior early-cap local reports were conservatively counted; their timing boundary differs. No exact first step, recovery or cleanup.'))
result.sort(key=lambda r:r['cap'],reverse=True)
args.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
