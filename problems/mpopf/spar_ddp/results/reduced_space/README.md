# Reduced-space MPOPF diagnostics

Raw results for the phase-1 structural test of the proposed decomposition
(outer: battery only; inner: single-period network OPF with `P_B^t` fixed).
Interpretation, caveats and the answers to the eight study questions are in
[`ddp/notes/REDUCED_SPACE_INNER_OPF_FEASIBILITY.md`](../../notes/REDUCED_SPACE_INNER_OPF_FEASIBILITY.md).

| file | contents |
|---|---|
| `inner_opf_survey_<system>_T<T>.csv` | one row per (time sample, dispatch pattern): status, feasibility, objective split, `P_Subs`/`Q_Subs`, loss, reverse-export flag, voltage envelope, `ell` magnitudes, reactive utilisation, active-bound counts, bound margin, FilterDDP cross-validation residual, SOC-box check, solve time and memory |
| `reverse_export_scan_<system>_T<T>.csv` | `all_max_discharge` at every hour — closest approach of `P_Subs` to its zero floor |
| `phi_probe_<system>_T<T>.csv` | `Phi_t` probes: base, direction, step size, value, first/second central differences, the analytic dual gradient, their relative error, and active-constraint counts |
| `infeasible_restoration_<system>_T<T>.csv` | written **only if** an inner solve is infeasible: which physical class is violated and by how much. Never a feasibility claim for the original problem. Absent for IEEE123 because nothing was infeasible. |

Scope so far: `ieee123C_1ph`, `T = 24` only. Regenerate with

```bash
julia --startup-file=no --project=envs/ddp2026 \
  ddp/examples/power_system/run_inner_opf_survey.jl ieee123C_1ph 24 24
julia --startup-file=no --project=envs/ddp2026 \
  ddp/examples/power_system/probe_reduced_value_function.jl ieee123C_1ph 24
```

Sign convention throughout: `+P_B` is discharge (injection into the bus).
