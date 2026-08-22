# FilterDDP on IEEE123C with full line losses (2026-08-21)

The nonlinear branch-flow model includes the receiving-end loss terms

```text
sum(P_jk) - P_ij + r_ij*ell_ij = P_B,j + p_D,j - p_L,j
sum(Q_jk) - Q_ij + x_ij*ell_ij = q_D,j - q_L,j
```

together with the nonlinear voltage-drop and SOCP equations.  An earlier
transcription omitted the two loss terms; all results from that transcription
are superseded.

## Corrected IEEE123C result

For `ieee123C_1ph`, `T = 3`:

- buses / branches: 128 / 127
- batteries / states (`nx`): 51
- per-stage controls (`nu`): 791
- per-stage equalities (`nc`): 562
- reduced-space dimension (`nu - nc`): 229
- centralized JuMP/Ipopt objective: `2808.924651223220`
- centralized Ipopt iterations: 34
- centralized validator: all constraints satisfied
- FilterDDP objective: `2808.924770488856`
- FilterDDP iterations: 48
- FilterDDP maximum equality residual: `6.227e-08`
- objective gap from Ipopt: `1.193e-04`
- maximum battery-energy difference from Ipopt: `1.157e-05` p.u.-h
- maximum branch-power differences: `1.269e-06` p.u. real and
  `1.716e-06` p.u. reactive

Model construction took 21.408 s, dynamic solver storage 0.524 s, and
`solve!` 70.998 s in the recorded run.

The corresponding tADMM run converged in 30 iterations with objective
`2808.58`; its lossy branch-flow validator passed every constraint.

## Implementation

The network path uses dynamic `Vector`/`Matrix` solver workspaces, LAPACK-backed
factorizations, and analytic network derivatives.  The analytic Jacobian
contains `r_ij` and `x_ij` derivatives in the real and reactive balance rows,
respectively.

Reproduction scripts:

- `envs/tadmm/root_level/run_bf.jl`
- `envs/tadmm/root_level/run_tadmm.jl`
- `ddp/examples/power_system/export_ieee123c_data.jl`
- `ddp/examples/power_system/ieee123c_filterddp.jl`
