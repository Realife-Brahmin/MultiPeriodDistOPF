# Project Development Context

## Current Status (Updated: 2025-10-29)

### Completed Features
- [x] **tADMM for LinDistFlow MPOPF** - Full implementation with temporal decomposition
- [x] **Adaptive ρ scheme** - Residual balancing for faster convergence (optional)
- [x] **Blue variable extraction** - Direct diagonal extraction from subproblem solutions
- [x] **Complete plot parity** - All BF plots have tADMM equivalents
- [x] **Convergence diagnostics** - Primal/dual residual tracking with threshold visualization
- [x] **System validation** - OpenDSS powerflow comparison and KVL verification
- [x] **PV generation profile** - Bell curve in middle 50% of horizon, zeros at edges
- [x] **All battery plotting** - Optional saving of individual plots for all batteries

### Key Implementation Details
- **System**: LinDistFlow constraints with voltage limits and battery storage
- **Penalty parameter**: ρ = 10000.0 (high for strong coordination)
- **Adaptive ρ**: Optional residual balancing (μ=10, τ=2, update every 10 iterations)
- **Convergence tolerances**: eps_pri = 1e-5, eps_dual = 1e-4
- **Battery convention**: P_B > 0 = discharging, P_B < 0 = charging
- **Consensus variable**: Battery SOC (B_j^{t0}[t]) - only variable coordinated across subproblems
- **Solution extraction**: Diagonal approach - subproblem t0 provides solution at time t0

### Current Working Directory
```
envs/tadmm/
├── tadmm.jl                      # Main tADMM implementation with LinDistFlow
├── Plotter.jl                    # Plotting utilities (input curves, batteries, voltages, PV)
├── parse_opendss.jl              # System data parser from OpenDSS
├── opendss_validator.jl          # OpenDSS powerflow validation
├── logger.jl                     # Logging utilities
└── processedData/                # Generated plots and results by system and horizon
    └── {systemName}_T{T}/
        ├── input_curves.png
        ├── battery_actions_*.png
        ├── substation_power_cost_*.png
        ├── voltage_*.png/gif
        ├── pv_power_*.png/gif
        └── convergence/
            └── tadmm_convergence.png

tex/
├── tadmm_formulation.tex         # tADMM mathematical formulation (LaTeX)
```

### Next Steps
- [ ] Analyze convergence behavior for different system sizes (IEEE 123, 730, etc.)
- [ ] Compare tADMM vs BF computational times and solution quality
- [ ] Test scalability with T=48, T=96 (larger horizons)
- [ ] Document optimal ρ selection strategies for different networks

### Development Notes
- **Subproblem structure**: Each subproblem t0 optimizes network at time t0 + full battery trajectory
- **Penalty term**: (ρ/2)||B^{t0} - B̂ + u^{t0}||² couples battery schedules across time
- **Consensus update**: B̂[t] = (1/T) × Σ_{t0} (B^{t0}[t] + u^{t0}[t]) - includes dual variables
- **Dual update**: Standard ADMM dual ascent u^{t0}[t] += (B^{t0}[t] - B̂[t])
- **High ρ**: Fast primal convergence (agreement), slow dual convergence (oscillation)
- **Adaptive ρ**: Automatically reduces when dual residual stalls, increases when primal lags
- **PV profile**: 25% zeros, 50% sin bell curve, 25% zeros - mimics solar sunrise/noon/sunset
- **Plot settings**: `showPlots=false`, `saveAllBatteryPlots=false` for production runs

### Recent Technical Decisions (Oct 27-29, 2025)
1. **Diagonal extraction**: Subproblem t provides solution at time t (no reconstruction)
2. **Full trajectory storage**: Store P_B[j][t0][t] for all times, extract diagonal
3. **Adaptive ρ enabled by default**: `adaptive_rho_tadmm = true` for faster convergence
4. **Complete plot parity**: Every BF plot has tADMM equivalent (voltages, PV, batteries, etc.)
5. **Dictionary structure**: Nested dicts for proper (battery, subproblem, time) indexing
6. **No divide-by-zero**: PV curve naturally in [0,1] without normalization
7. **Convergence threshold display**: Red dashed lines on convergence plots