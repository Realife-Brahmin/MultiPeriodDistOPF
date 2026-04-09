# ============================================================================
# SHARED CONFIGURATION - Edit this file, then run run_bf.jl or run_tadmm.jl
# ============================================================================

# System and time horizon
const SYSTEM_NAME = "large10kC_1ph"  # "ieee123A_1ph", "ieee2552_1ph", "large10k_1ph", "large10kB_1ph", "large10kC_1ph"
const T = 24                         # Time periods: 4, 6, 12, 24, 48, etc.
const DELTA_T_H = 24.0 / T           # Time step duration in hours

# Solver selection (Gurobi requires license)
const USE_GUROBI = false
const USE_GUROBI_FOR_BF = false
const USE_GUROBI_FOR_TADMM = false

# BF-specific
const BF_TIMEOUT_SEC = 0  # Kill BF after this many seconds (0 = no limit)

# Verbose output
const VERBOSE = false

# ============================================================================
# LOAD SHAPES (deterministic profiles over T periods)
# ============================================================================

# Load profile: sinusoidal variation around 0.8-1.0
const LoadShapeLoad = 0.8 .+ 0.2 .* (sin.(range(0, 2pi, length=T) .- 0.8) .+ 1) ./ 2

# Energy cost profile: time-varying $/kWh
const LoadShapeCost = 0.08 .+ 0.12 .* (sin.(range(0, 2pi, length=T)) .+ 1) ./ 2

# Solar PV profile: bell curve in middle 50% of horizon
const LoadShapePV = let
    pv = zeros(T)
    if T >= 4
        t_start = max(1, round(Int, 0.25 * T))
        t_end = min(T, round(Int, 0.75 * T))
        n_active = t_end - t_start + 1
        if n_active >= 2
            t_normalized = range(0, pi, length=n_active)
            pv[t_start:t_end] = [sin(t_norm) for t_norm in t_normalized]
        else
            pv[t_start:t_end] .= 1.0
        end
    elseif T == 3
        pv[2] = 1.0
    elseif T == 2
        pv .= 0.5
    else
        pv .= 1.0
    end
    pv
end

# Battery quadratic cost coefficient
const C_B = 1e-6 * minimum(LoadShapeCost)

# ============================================================================
# OUTPUT DIRECTORIES
# ============================================================================

const ENV_PATH = joinpath(@__DIR__, "envs", "tadmm")
const PROCESSED_DATA_DIR = joinpath(@__DIR__, "envs", "tadmm", "processedData")
const SYSTEM_DIR = joinpath(PROCESSED_DATA_DIR, "$(SYSTEM_NAME)_T$(T)")
