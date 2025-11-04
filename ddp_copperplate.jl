# ============================================================================
# DDP FOR COPPER PLATE PROBLEM
# ============================================================================
# This script solves a simplified multi-period optimal power flow (MPOPF)
# problem using the Copper Plate approximation:
# - No network constraints (single aggregated bus)
# - Only power balance constraint: P_Subs + P_B = P_L
# - Battery dynamics: B[t] = B[t-1] - P_B[t] * Δt
# - Objective: minimize energy cost + battery quadratic cost
#
# Solution methods:
# 1. Brute Force: Centralized optimization (all T time steps together)
# 2. DDP: Distributed decomposition by time period (TO BE IMPLEMENTED)
# ============================================================================

# ============================================================================
# PARALLELIZATION CONFIGURATION
# ============================================================================
function detect_cpu_cores()
    if haskey(ENV, "NUMBER_OF_PROCESSORS")
        return parse(Int, ENV["NUMBER_OF_PROCESSORS"])
    end
    sys_threads = Sys.CPU_THREADS
    if sys_threads > 1
        return sys_threads
    end
    return max(1, Threads.nthreads())
end

const AVAILABLE_CORES = detect_cpu_cores()
const USE_THREADING = true
# ============================================================================

# Activate the DDP environment
import Pkg
env_path = joinpath(@__DIR__, "envs", "ddp")
Pkg.activate(env_path)

using Revise
using LinearAlgebra
using JuMP
using Ipopt
using Gurobi
using Printf
using Statistics
using Parameters: @unpack
using Plots
using Crayons

# Display threading status at startup
println("\n" * "="^80)
println("🚀 PARALLELIZATION STATUS")
println("="^80)
println("Hardware CPU cores:  ", AVAILABLE_CORES)
println("Julia threads:       ", Threads.nthreads())
if Threads.nthreads() >= AVAILABLE_CORES
    println("Status:              ✓ Using all available cores")
elseif Threads.nthreads() > 1
    println("Status:              ⚠  Using $(Threads.nthreads())/$(AVAILABLE_CORES) cores")
else
    println("Status:              ⚠  Single-threaded mode")
end
println("="^80)

# Create shared Gurobi environment (suppresses repeated license messages)
const GUROBI_ENV = Gurobi.Env()

# Define color schemes
const COLOR_SUCCESS = Crayon(foreground = :green, bold = true)
const COLOR_WARNING = Crayon(foreground = :yellow, bold = true)
const COLOR_ERROR = Crayon(foreground = :red, bold = true)
const COLOR_INFO = Crayon(foreground = :cyan, bold = true)
const COLOR_HIGHLIGHT = Crayon(foreground = :magenta, bold = true)
const COLOR_RESET = Crayon(reset = true)

begin # entire script

# ============================================================================
# PROBLEM DATA STRUCTURE
# ============================================================================

"""
    BatteryParams

Battery parameters in per-unit system:
- E_Rated_pu: Rated energy capacity (pu on E_BASE)
- soc_min, soc_max: SOC limits (fraction of E_Rated_pu)
- P_B_R_pu: Rated power (pu on P_BASE)
- Δt: Time step duration (hours)
- B0_pu: Initial SOC (pu on E_BASE)
- B_T_target_pu: Terminal SOC constraint (optional)
"""
struct BatteryParams
    E_Rated_pu::Float64
    soc_min::Float64
    soc_max::Float64
    P_B_R_pu::Float64
    Δt::Float64
    B0_pu::Float64
    B_T_target_pu::Union{Nothing,Float64}
end

# Helper functions for battery SOC bounds
Bmin(b::BatteryParams) = b.soc_min * b.E_Rated_pu
Bmax(b::BatteryParams) = b.soc_max * b.E_Rated_pu

"""
    CopperPlateData

Complete problem instance for copper plate MPOPF:
- T: Number of time periods
- price: Energy price (\$/kWh) for each time period
- P_L_pu: Load demand (pu on P_BASE) for each time period
- bat: Battery parameters
- P_BASE: Power base (kW)
- E_BASE: Energy base (kWh)
- C_B: Battery quadratic cost coefficient (\$/kW²)
"""
struct CopperPlateData
    T::Int
    price::Vector{Float64}
    P_L_pu::Vector{Float64}
    bat::BatteryParams
    P_BASE::Float64
    E_BASE::Float64
    C_B::Float64
end

# ============================================================================
# DATA GENERATION
# ============================================================================

"""
    setup_copperplate_data(; kwargs...)

Generate synthetic copper plate problem data with realistic load/cost profiles.

Arguments:
- T: Number of time periods (default: 24)
- P_L_R_kW: Rated load power in kW (default: 1250.0)
- E_Rated_kWh: Battery energy capacity in kWh (default: 4000.0)
- P_B_R_kW: Battery power rating in kW (default: 800.0)
- soc_min, soc_max: SOC limits (default: 0.30, 0.95)
- B0_fraction: Initial SOC as fraction of E_Rated (default: 0.60)
- B_T_target_fraction: Terminal SOC constraint (default: nothing)
- Δt_h: Time step duration in hours (default: 1.0)

Returns: CopperPlateData instance
"""
function setup_copperplate_data(;
    T::Int=24,
    P_L_R_kW::Float64=1250.0,
    E_Rated_kWh::Float64=4000.0,
    P_B_R_kW::Float64=800.0,
    soc_min::Float64=0.30,
    soc_max::Float64=0.95,
    B0_fraction::Float64=0.60,
    B_T_target_fraction::Union{Nothing,Float64}=nothing,
    Δt_h::Float64=1.0
)
    # Physical bases
    KVA_B = 1000.0
    P_BASE = KVA_B  # kW
    E_BASE = P_BASE * 1.0  # kWh per 1 hour
    
    # Generate load shape: sinusoidal profile normalized to [0.6, 1.0]
    LoadShapeLoad = 0.8 .+ 0.2 .* (sin.(range(0, 2π, length=T) .- 0.8) .+ 1) ./ 2
    P_L_kW = P_L_R_kW .* LoadShapeLoad
    
    # Generate cost shape: sinusoidal profile from 8 to 20 cents/kWh
    LoadShapeCost = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T)) .+ 1) ./ 2
    
    # Battery quadratic cost coefficient
    C_B = 1e-6 * minimum(LoadShapeCost)
    
    # Convert to per-unit
    P_L_pu = P_L_kW ./ P_BASE
    E_Rated_pu = E_Rated_kWh / E_BASE
    P_B_R_pu = P_B_R_kW / P_BASE
    B0_pu = B0_fraction * E_Rated_pu
    B_T_target_pu = isnothing(B_T_target_fraction) ? nothing : (B_T_target_fraction * E_Rated_pu)
    
    bat = BatteryParams(E_Rated_pu, soc_min, soc_max, P_B_R_pu, Δt_h, B0_pu, B_T_target_pu)
    
    return CopperPlateData(T, LoadShapeCost, P_L_pu, bat, P_BASE, E_BASE, C_B)
end

# ============================================================================
# BRUTE FORCE SOLVER (CENTRALIZED)
# ============================================================================

"""
    solve_CopperPlate_BruteForced(data; solver=:gurobi)

Solve the copper plate MPOPF problem using centralized optimization.
All T time periods are solved together as a single large optimization problem.

Decision variables (T time periods):
- P_Subs[t]: Substation power at time t (≥ 0)
- P_B[t]: Battery power at time t (charge < 0, discharge > 0)
- B[t]: Battery SOC at end of time t

Constraints:
- Power balance: P_Subs[t] + P_B[t] = P_L[t]  ∀t
- Battery dynamics: B[t] = B[t-1] - P_B[t] * Δt  ∀t
- Initial SOC: B[0] = B0
- SOC bounds: B_min ≤ B[t] ≤ B_max  ∀t
- Power bounds: -P_B_R ≤ P_B[t] ≤ P_B_R  ∀t
- Terminal SOC (optional): B[T] = B_T_target

Objective: minimize Σₜ (price[t] * P_Subs[t] * Δt + C_B * P_B[t]² * Δt)

Arguments:
- data: CopperPlateData instance
- solver: :gurobi or :ipopt

Returns: Dict with solution and statistics
"""
function solve_CopperPlate_BruteForced(data::CopperPlateData; solver::Symbol=:gurobi)
    
    @unpack T, price, P_L_pu, bat, P_BASE, E_BASE, C_B = data
    Δt = bat.Δt
    
    # Create model
    model = Model()
    if solver == :gurobi
        set_optimizer(model, () -> Gurobi.Optimizer(GUROBI_ENV))
        set_silent(model)
        set_optimizer_attribute(model, "OutputFlag", 0)
        set_optimizer_attribute(model, "DualReductions", 0)
    else  # :ipopt
        set_optimizer(model, Ipopt.Optimizer)
        set_silent(model)
        set_optimizer_attribute(model, "max_iter", 3000)
    end
    
    # Decision variables
    @variable(model, P_Subs[t in 1:T] >= 0)
    @variable(model, -bat.P_B_R_pu <= P_B[t in 1:T] <= bat.P_B_R_pu)
    @variable(model, Bmin(bat) <= B[t in 1:T] <= Bmax(bat))
    
    # Constraints
    # Power balance at each time
    @constraint(model, power_balance[t in 1:T], P_Subs[t] + P_B[t] == P_L_pu[t])
    
    # Battery dynamics
    @constraint(model, B[1] == bat.B0_pu - P_B[1] * Δt)
    @constraint(model, [t in 2:T], B[t] == B[t-1] - P_B[t] * Δt)
    
    # Terminal SOC constraint (if specified)
    if !isnothing(bat.B_T_target_pu)
        @constraint(model, B[T] == bat.B_T_target_pu)
    end
    
    # Objective: energy cost + battery quadratic cost
    @objective(model, Min, 
        sum(price[t] * (P_Subs[t] * P_BASE) * Δt + 
            C_B * (P_B[t] * P_BASE)^2 * Δt 
            for t in 1:T)
    )
    
    # Solve
    wallclock_start = time()
    optimize!(model)
    wallclock_time = time() - wallclock_start
    
    # Extract results
    status = termination_status(model)
    obj_val = has_values(model) ? objective_value(model) : NaN
    solver_time = solve_time(model)
    
    # Get model statistics
    n_vars = num_variables(model)
    n_constraints = sum(num_constraints(model, F, S) 
                       for (F, S) in list_of_constraint_types(model))
    
    result = Dict(
        :status => status,
        :objective => obj_val,
        :solve_time => solver_time,
        :wallclock_time => wallclock_time,
        :solver => solver,
        :model_stats => Dict(
            :n_variables => n_vars,
            :n_constraints => n_constraints
        )
    )
    
    # Extract solution arrays if optimization succeeded
    if has_values(model)
        result[:P_Subs] = value.(P_Subs)
        result[:P_B] = value.(P_B)
        result[:B] = value.(B)
        
        # Compute objective components
        energy_cost = sum(price[t] * (value(P_Subs[t]) * P_BASE) * Δt for t in 1:T)
        battery_cost = sum(C_B * (value(P_B[t]) * P_BASE)^2 * Δt for t in 1:T)
        
        result[:energy_cost] = energy_cost
        result[:battery_cost] = battery_cost
    end
    
    return result
end

# ============================================================================
# SCENARIO CONFIGURATION
# ============================================================================

begin # scenario config
    # System parameters
    T = 24  # Number of time periods
    
    # Solver selection
    use_gurobi_for_bf = true  # Use Gurobi for brute force
    
    # Battery and load parameters (will be passed to setup function)
    P_L_R_kW = 1250.0       # Rated load power
    E_Rated_kWh = 4000.0    # Battery capacity
    P_B_R_kW = 800.0        # Battery power rating
    soc_min = 0.30
    soc_max = 0.95
    B0_fraction = 0.60      # Initial SOC (60% of capacity)
    Δt_h = 1.0              # 1-hour time steps
    
    # Plotting settings
    showPlots = false
    savePlots = true
end # scenario config

# ============================================================================
# SETUP PROBLEM DATA
# ============================================================================

println("\n" * "="^80)
println("SETTING UP COPPER PLATE PROBLEM DATA")
println("="^80)

data = setup_copperplate_data(
    T=T,
    P_L_R_kW=P_L_R_kW,
    E_Rated_kWh=E_Rated_kWh,
    P_B_R_kW=P_B_R_kW,
    soc_min=soc_min,
    soc_max=soc_max,
    B0_fraction=B0_fraction,
    Δt_h=Δt_h
)

println("\nProblem Configuration:")
println("  Time periods (T):         ", T)
println("  Time step (Δt):           ", Δt_h, " hours")
println("  Rated load (P_L_R):       ", P_L_R_kW, " kW")
println("  Battery capacity (E_R):   ", E_Rated_kWh, " kWh")
println("  Battery power (P_B_R):    ", P_B_R_kW, " kW")
println("  SOC limits:               [", soc_min*100, "%, ", soc_max*100, "%]")
println("  Initial SOC (B₀):         ", B0_fraction*100, "%")
println("  Power base (P_BASE):      ", data.P_BASE, " kW")
println("  Energy base (E_BASE):     ", data.E_BASE, " kWh")

# Print load and cost profiles summary
P_L_kW = data.P_L_pu .* data.P_BASE
price_cents = data.price .* 100
println("\nLoad Profile:")
@printf "  Min load:  %.1f kW (%.1f%% of rated)\n" minimum(P_L_kW) (minimum(P_L_kW)/P_L_R_kW*100)
@printf "  Max load:  %.1f kW (%.1f%% of rated)\n" maximum(P_L_kW) (maximum(P_L_kW)/P_L_R_kW*100)
@printf "  Avg load:  %.1f kW (%.1f%% of rated)\n" mean(P_L_kW) (mean(P_L_kW)/P_L_R_kW*100)

println("\nCost Profile:")
@printf "  Min price: %.1f cents/kWh\n" minimum(price_cents)
@printf "  Max price: %.1f cents/kWh\n" maximum(price_cents)
@printf "  Avg price: %.1f cents/kWh\n" mean(price_cents)

println("="^80)

# ============================================================================
# SOLVE USING BRUTE FORCE (CENTRALIZED)
# ============================================================================

println("\n" * "="^80)
println(COLOR_HIGHLIGHT, "SOLVING COPPER PLATE MPOPF (BRUTE-FORCED)", COLOR_RESET)
println("="^80)

solver_bf_choice = use_gurobi_for_bf ? :gurobi : :ipopt
sol_bf = solve_CopperPlate_BruteForced(data; solver=solver_bf_choice)

# Report results
println("\n--- SOLUTION STATUS ---")
print("Status: ")

if sol_bf[:status] == MOI.OPTIMAL || sol_bf[:status] == MOI.LOCALLY_SOLVED
    print(COLOR_SUCCESS)
    println(sol_bf[:status])
    println("✓ Optimization successful!")
    print(COLOR_RESET)
    
    println("\n--- OBJECTIVE VALUE ---")
    @printf "Total Cost:     \$%.4f\n" sol_bf[:objective]
    @printf "  Energy Cost:  \$%.4f (%.1f%%)\n" sol_bf[:energy_cost] (sol_bf[:energy_cost]/sol_bf[:objective]*100)
    @printf "  Battery Cost: \$%.4f (%.1f%%)\n" sol_bf[:battery_cost] (sol_bf[:battery_cost]/sol_bf[:objective]*100)
    
    println("\n--- COMPUTATION TIME ---")
    @printf "Solver time:    %.4f seconds\n" sol_bf[:solve_time]
    @printf "Wallclock time: %.4f seconds\n" sol_bf[:wallclock_time]
    
    println("\n--- MODEL STATISTICS ---")
    @printf "Variables:      %d\n" sol_bf[:model_stats][:n_variables]
    @printf "Constraints:    %d\n" sol_bf[:model_stats][:n_constraints]
    
    # Convert to physical units
    P_Subs_kW = sol_bf[:P_Subs] .* data.P_BASE
    P_B_kW = sol_bf[:P_B] .* data.P_BASE
    B_kWh = sol_bf[:B] .* data.E_BASE
    
    println("\n--- POWER SUMMARY ---")
    @printf "Substation (kW):  min=%.1f, max=%.1f, avg=%.1f\n" minimum(P_Subs_kW) maximum(P_Subs_kW) mean(P_Subs_kW)
    @printf "Battery (kW):     min=%.1f, max=%.1f, avg=%.1f\n" minimum(P_B_kW) maximum(P_B_kW) mean(P_B_kW)
    
    println("\n--- SOC SUMMARY ---")
    soc_percent = B_kWh ./ E_Rated_kWh .* 100
    @printf "SOC (%%):          min=%.1f, max=%.1f, avg=%.1f\n" minimum(soc_percent) maximum(soc_percent) mean(soc_percent)
    @printf "SOC (kWh):        min=%.1f, max=%.1f, avg=%.1f\n" minimum(B_kWh) maximum(B_kWh) mean(B_kWh)
    
    # Verify power balance
    println("\n--- POWER BALANCE VERIFICATION ---")
    P_L_kW = data.P_L_pu .* data.P_BASE
    balance_error = maximum(abs.(P_Subs_kW .+ P_B_kW .- P_L_kW))
    @printf "Max power balance error: %.2e kW\n" balance_error
    if balance_error < 1e-6
        print(COLOR_SUCCESS)
        println("✓ Power balance satisfied")
        print(COLOR_RESET)
    else
        print(COLOR_WARNING)
        println("⚠ Power balance may have numerical errors")
        print(COLOR_RESET)
    end
    
    # Verify battery dynamics
    println("\n--- BATTERY DYNAMICS VERIFICATION ---")
    B0_kWh = data.bat.B0_pu * data.E_BASE
    B_expected = [B0_kWh]
    for t in 1:T
        push!(B_expected, B_expected[end] - P_B_kW[t] * Δt_h)
    end
    dynamics_error = maximum(abs.(B_kWh .- B_expected[2:end]))
    @printf "Max battery dynamics error: %.2e kWh\n" dynamics_error
    if dynamics_error < 1e-6
        print(COLOR_SUCCESS)
        println("✓ Battery dynamics satisfied")
        print(COLOR_RESET)
    else
        print(COLOR_WARNING)
        println("⚠ Battery dynamics may have numerical errors")
        print(COLOR_RESET)
    end
    
else
    print(COLOR_ERROR)
    println("⚠ Optimization failed or did not converge to optimality")
    println("Status: ", sol_bf[:status])
    print(COLOR_RESET)
end

println("\n" * "="^80)
println(COLOR_HIGHLIGHT, "COPPER PLATE MPOPF BRUTE-FORCED SOLUTION COMPLETE", COLOR_RESET)
println("="^80)

# ============================================================================
# SAVE RESULTS TO FILE
# ============================================================================

if sol_bf[:status] == MOI.OPTIMAL || sol_bf[:status] == MOI.LOCALLY_SOLVED
    processedData_dir = joinpath(@__DIR__, "envs", "ddp", "processedData")
    system_folder = "copperplate_T$(T)"
    system_dir = joinpath(processedData_dir, system_folder)
    mkpath(system_dir)
    
    results_file = joinpath(system_dir, "results_bf.txt")
    
    open(results_file, "w") do io
        println(io, "="^80)
        println(io, "COPPER PLATE MPOPF - BRUTE FORCE SOLUTION")
        println(io, "="^80)
        println(io)
        
        println(io, "Problem Configuration:")
        println(io, "  Time periods (T):         ", T)
        println(io, "  Solver:                   ", sol_bf[:solver])
        println(io, "  Optimization status:      ", sol_bf[:status])
        println(io)
        
        println(io, "Objective Value:")
        @printf(io, "  Total Cost:     \$%.6f\n", sol_bf[:objective])
        @printf(io, "  Energy Cost:    \$%.6f (%.2f%%)\n", sol_bf[:energy_cost], (sol_bf[:energy_cost]/sol_bf[:objective]*100))
        @printf(io, "  Battery Cost:   \$%.6f (%.2f%%)\n", sol_bf[:battery_cost], (sol_bf[:battery_cost]/sol_bf[:objective]*100))
        println(io)
        
        println(io, "Computation Time:")
        @printf(io, "  Solver time:    %.6f seconds\n", sol_bf[:solve_time])
        @printf(io, "  Wallclock time: %.6f seconds\n", sol_bf[:wallclock_time])
        println(io)
        
        println(io, "Model Statistics:")
        @printf(io, "  Variables:      %d\n", sol_bf[:model_stats][:n_variables])
        @printf(io, "  Constraints:    %d\n", sol_bf[:model_stats][:n_constraints])
        println(io)
        
        P_Subs_kW = sol_bf[:P_Subs] .* data.P_BASE
        P_B_kW = sol_bf[:P_B] .* data.P_BASE
        B_kWh = sol_bf[:B] .* data.E_BASE
        
        println(io, "Solution Summary:")
        @printf(io, "  Substation (kW):  min=%.2f, max=%.2f, avg=%.2f\n", minimum(P_Subs_kW), maximum(P_Subs_kW), mean(P_Subs_kW))
        @printf(io, "  Battery (kW):     min=%.2f, max=%.2f, avg=%.2f\n", minimum(P_B_kW), maximum(P_B_kW), mean(P_B_kW))
        @printf(io, "  SOC (kWh):        min=%.2f, max=%.2f, avg=%.2f\n", minimum(B_kWh), maximum(B_kWh), mean(B_kWh))
        println(io)
        
        println(io, "="^80)
    end
    
    println(COLOR_SUCCESS, "✓ Results written to $(results_file)", COLOR_RESET)
end

end # entire script
