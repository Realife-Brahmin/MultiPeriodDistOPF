# tADMM: Temporal ADMM for Multi-Period Optimal Power Flow

Standalone implementation of temporal ADMM for solving multi-period optimal power flow problems.

## Directory Structure

```
envs/tadmm/
├── Project.toml              # Julia environment dependencies
├── tadmm.jl                  # Main script
├── parse_opendss.jl          # Standalone OpenDSS parser
├── opendss_validator.jl      # OpenDSS validation utilities
└── README.md                 # This file
```

## Features

- **Standalone**: No dependency on MultiPeriodDistOPF package
- **Minimal Dependencies**: Only essential Julia packages
- **LinDistFlow Formulation**: Linearized power flow constraints
- **Battery Optimization**: SOC trajectory with quadratic cost
- **Network Constraints**: 
  - Nodal power balance (real & reactive)
  - KVL using LinDistFlow approximation
  - Voltage limits
  - PV reactive power limits

## Dependencies

- `JuMP` - Optimization modeling
- `Ipopt` - Nonlinear solver (default)
- `Gurobi` - Alternative solver (optional)
- `OpenDSSDirect` - OpenDSS interface
- `Parameters` - @unpack macro
- `Statistics` - Basic statistics
- Standard library: `Printf`, `Random`, `Dates`

## Installation

```julia
# Activate the environment
using Pkg
Pkg.activate("envs/tadmm")
Pkg.instantiate()
```

## Usage

### Basic Usage

```julia
# From the project root directory
julia envs/tadmm/tadmm.jl
```

### Customization

Edit the user configuration section in `tadmm.jl`:

```julia
# System and simulation parameters
systemName = "ads10A_1ph"      # System to simulate
T = 24                          # Number of time steps
delta_t_h = 1.0                 # Time step duration (hours)

# Load shapes (0-1 normalized)
LoadShapeLoad = ...             # Load multiplier per timestep
LoadShapeCost = ...             # Energy cost ($/kWh) per timestep
C_B = ...                       # Battery quadratic cost coefficient
```

## Formulation

### Objective Function
```
min Σ_t [LoadShapeCost[t] * P_Subs[t] * Δt + C_B * P_B[t]^2 * Δt]
```

### Key Constraints

1. **Power Balance** (all non-substation nodes):
   ```
   P_ij = Σ_k P_jk + p_L - p_D - P_B
   Q_ij = Σ_k Q_jk + q_L - q_D
   ```

2. **LinDistFlow (KVL)**:
   ```
   v_i - v_j = 2(r_ij * P_ij + x_ij * Q_ij)
   ```

3. **Battery SOC**:
   ```
   B[t] = B[t-1] - P_B[t-1] * Δt
   soc_min * B_R <= B[t] <= soc_max * B_R
   ```

4. **Voltage Limits**:
   ```
   V_min^2 <= v[n,t] <= V_max^2
   ```

## Output

The script prints:
- OpenDSS powerflow validation
- Optimization status
- Total cost
- Substation power summary (min/max/avg)
- Battery power and SOC
- Voltage summary

## Differences from Main Package

This implementation:
- ✅ Standalone (no MPDOPF.jl dependency)
- ✅ Simplified parser (essential data only)
- ✅ Direct script execution (no module system)
- ✅ Minimal dependencies
- ❌ No advanced features (DDP, spatial decomposition, etc.)
- ❌ No model saving/loading
- ❌ No plotting utilities

## Next Steps

Future enhancements:
1. Implement temporal decomposition (ADMM)
2. Add dual decomposition (DDP)
3. Add result plotting
4. Add more test cases
5. Benchmark against centralized solution

## Troubleshooting

### Infeasible Problem

If the optimization returns `LOCALLY_INFEASIBLE`:
- Check that all nodes in `Nm1set` have proper power balance
- Verify battery SOC limits are reasonable
- Check voltage limits are not too tight
- Ensure PV reactive power limits are feasible

### Parser Errors

If OpenDSS parsing fails:
- Verify `rawData` directory structure
- Check that `Master.dss` exists for the system
- Ensure OpenDSS can load the system

## References

- LinDistFlow: Baran & Wu (1989)
- ADMM: Boyd et al. (2011)
- MultiPeriodDistOPF: Parent repository
