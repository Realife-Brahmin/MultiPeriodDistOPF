"""Extract observed penalty decisions from completed, validated campaign cases.

This is posthoc diagnosis, never a tuning change or a replacement solver run.
"""
import csv
import hashlib
import json
from pathlib import Path
import sys
try:
    import tomllib
except ImportError:
    import tomli as tomllib

campaign = Path(sys.argv[1]).resolve()
out = Path(sys.argv[2]).resolve()
out.mkdir(parents=True, exist_ok=False)
summary_rows, events, selected = [], [], []
for case in sorted((campaign / 'results').iterdir()):
    if not case.name.startswith(('ieee2522C_1ph_', 'large10kC_1ph_')):
        continue
    if not (case / 'COMPLETE').is_file() or 'adaptive' in case.name:
        continue
    val = json.loads((case / 'validation.json').read_text())
    assert all(val.get(k) == 'pass' for k in (
        'optimization_validation', 'warmup_trace_validation',
        'opendss_validation', 'source_snapshot_validation'))
    meta = tomllib.loads((case / 'metadata.toml').read_text())
    s = tomllib.loads((case / 'summary.toml').read_text())
    trace = list(csv.DictReader((case / 'trace.csv').open()))
    assert len(trace) == s['iterations']
    ref = s['central_objective']
    near_tol = 1e-5  # Repo convention for med2522 and large10k; full solve retained.
    near = next((r for r in trace if float(r['original_primal']) <= near_tol
        and abs(float(r['objective']) - ref) / max(1, abs(ref)) <= .001), None)
    bad = [i for i, r in enumerate(trace) if float(r['original_primal']) > near_tol
        or abs(float(r['objective']) - ref) / max(1, abs(ref)) > .001]
    sustained_index = bad[-1] + 1 if bad else 0
    sustained = trace[sustained_index] if sustained_index < len(trace) else None
    changes = 0
    for prev, current in zip(trace, trace[1:]):
        old, new = float(prev['rho']), float(current['rho'])
        if old == new:
            continue
        changes += 1
        events.append(dict(case=case.name, phase='measured',
            decision_after_iteration=int(prev['iteration']),
            first_iteration_using_new_rho=int(current['iteration']),
            rho_old=old, rho_new=new, scaled_dual_rescale=old/new,
            original_primal_before=float(prev['original_primal']),
            all_copy_residual_before=float(prev['r_max']),
            dual_residual_before=float(prev['s_max']),
            relative_primal_before=float(prev['r_relative']),
            relative_dual_before=float(prev['s_relative'])))
    summary_rows.append(dict(case=case.name, system=s['system'], T=s['T'],
        role='benchmark' if '_benchmark_' in case.name else 'shared-node verification',
        backend=meta['linear_solver'], rho_policy=meta['rho_policy'],
        rho_initial=float(trace[0]['rho']), rho_final=float(trace[-1]['rho']),
        rho_min=min(float(r['rho']) for r in trace),
        rho_max=max(float(r['rho']) for r in trace), rho_changes=changes,
        iterations=s['iterations'], full_loop_wall_s=s['optimizer_loop_wall_s'],
        central_wall_s=s['central_wall_s'], objective_gap=s['relative_objective_gap'],
        original_primal=s['tadmm_original_primal'],
        first_near_optimal_iteration=int(near['iteration']) if near else '',
        first_near_optimal_loop_wall_s=float(near['loop_wall_s']) if near else '',
        near_optimal_sustained_from_iteration=int(sustained['iteration']) if sustained else '',
        near_optimal_sustained_from_loop_wall_s=float(sustained['loop_wall_s']) if sustained else '',
        near_optimal_note='First crossing: 0.1% objective and original feasibility 1e-5; not full convergence or a persistence guarantee',
        input_sha256=meta['data_sha256'],
        model_sha256=meta['model.jl_sha256'], driver_sha256=meta['run.jl_sha256']))
    selected.append(case.relative_to(campaign).as_posix())
assert summary_rows, 'No validated completed medium/large cases'
for name, rows in (('cases.csv', summary_rows), ('penalty_changes.csv', events)):
    if rows:
        with (out / name).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rows[0])
            writer.writeheader()
            writer.writerows(rows)
(out / 'selected_cases.txt').write_text('\n'.join(selected) + '\n')
(out / 'analysis-source.sha256').write_text(hashlib.sha256(Path(__file__).read_bytes()).hexdigest() + '  penalty_diagnostics.py\n')
print(json.dumps(summary_rows, indent=2))
