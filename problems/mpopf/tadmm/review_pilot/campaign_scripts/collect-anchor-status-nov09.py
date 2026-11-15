"""Collect compact immutable logs; raw dispatch arrays stay on Kamiak."""
import csv
import hashlib
import json
from pathlib import Path
import tarfile
import tomli

root = Path('/home/aryan.jha/kamiak_jobs')
campaigns = ['tadmm_review_nov09_historical_anchors', 'tadmm_review_nov09_mumps_fixed']
summary = []
archive = root / 'historical-anchors-completed-nov09.tgz'
with tarfile.open(archive, 'w:gz') as tar:
    for name in campaigns:
        c = root / name
        for p in sorted(c.rglob('*')):
            if not p.is_file() or 'julia-depot' in p.parts or (c/'source') in p.parents:
                continue
            if p.suffix in ('.jls', '.tgz', '.gz'):
                continue
            tar.add(p, arcname=str(p.relative_to(root)), recursive=False)
        for case in sorted((c/'results').glob('run_*')):
            if not (case/'SWEEP_VALIDATED').exists():
                continue
            m = tomli.loads((case/'sweep.toml').read_text())
            audit = tomli.loads((c/'campaign'/f"T{m['T']}-{case.name.split('_')[-1]}-energy.toml").read_text())
            for rho in m['rhos']:
                arm = case/f'rho{int(rho)}'
                a = tomli.loads((arm/'assessment.toml').read_text())
                marks = tomli.loads((arm/'milestones.toml').read_text())
                trace = list(csv.DictReader((arm/'trace.csv').open()))
                stops = {}
                for eps in (1e-3, 1e-4):
                    for r in trace:
                        if float(r['owned_r_rms']) <= eps and (float(r['owned_s_rms']) <= 1e-2 or r['historical_stagnant'] == 'true'):
                            stops[str(eps)] = {k:r[k] for k in ('iteration','algorithm_s','owned_r_rms','r_rms','owned_s_rms','original_primal')}
                            break
                def energy(point):
                    v = audit[f'{arm.name}/{point}']
                    return {k:v[k] for k in ('integrated_energy_bound_kwh','max_SOC_bound_percentage_points','max_dynamics_error_kwh')}
                milestone = {}
                for k in ('paper_cost_1e-3','current_cost_1e-3','all_dual_1e-3','current_cost_1e-4'):
                    if k in marks:
                        v = marks[k]
                        milestone[k] = {f:v[f] for f in ('iteration','core_s','relative_objective_gap','original_primal')}
                        milestone[k].update(energy(v['snapshot']))
                validation = json.loads((case/'validation_complete.json').read_text())
                checks = list(validation['central_opendss'] for _ in range(1)) + list(validation['arms'][arm.name]['opendss'].values())
                summary.append(dict(case=case.name,campaign=name,system=m['system'],T=m['T'],solver=m['linear_solver'],
                    iterations=a['iterations'],core_s=a['core_s'],loop_wall_s=a['loop_wall_s'],
                    central_s=m['central_wall_s'],final_rms_met=a['final_rms_met'],
                    strict_physical_met=a['strict_criteria_met'],cost_gap=a['objective_gap'],
                    final_energy=energy('endpoint.jls'),historical_stops=stops,milestones=milestone,
                    dss_checks=len(checks),all_dss_converged=all(v['converged'] for v in checks),
                    all_dss_zero_voltage_violations=all(v['zero_voltage_violations'] for v in checks),
                    max_dss_voltage_difference=max(v['max_voltage_difference_pu'] for v in checks)))
result = root/'historical-anchors-summary-nov09.json'
result.write_text(json.dumps(summary,indent=2)+'\n')
print(result.read_text())
print('ARCHIVE_SHA256',hashlib.sha256(archive.read_bytes()).hexdigest())
