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

# Include plotting utilities
includet(joinpath(env_path, "Plotter.jl"))

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
    @constraint(model, B[1] - bat.B0_pu + P_B[1] * Δt == 0)
    @constraint(model, [t in 2:T], B[t] - B[t-1] + P_B[t] * Δt == 0)
    
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
# DDP SOLVER (FORWARD PASS DECOMPOSITION)
# ============================================================================

"""
    forward_pass_ddp!(B_collection, mu_collection, PB_collection, data, t; solver=:gurobi)

Solve DDP subproblem for time step t using forward pass decomposition.

Each time step optimizes only its own variables {P_Subs[t], P_B[t], B[t]} using:
- From previous step: B[t-1] (or B0 if t=1)
- From next step (last iteration): B[t+1], μ[t+1], P_B[t+1]

Objective:
    min C[t]*P_Subs[t] + C_B*P_B[t]² + μ[t+1]*(B[t+1] - B[t] + Δt*P_B[t+1])
        └─ energy cost ─┘  └ battery ┘   └────── coupling with next ──────┘

The coupling term enforces consistency with the next time step's dynamics.

Returns: Dict with solution and updates B_collection[t], PB_collection[t], mu_collection[t]
"""
function forward_pass_ddp!(B_collection, mu_collection, PB_collection, data::CopperPlateData, t::Int; solver::Symbol=:gurobi)
    
    @unpack T, price, P_L_pu, bat, P_BASE, E_BASE, C_B = data
    Δt = bat.Δt
    
    # Create model
    model = Model()
    if solver == :gurobi
        set_optimizer(model, () -> Gurobi.Optimizer(GUROBI_ENV))
        set_silent(model)
        set_optimizer_attribute(model, "OutputFlag", 0)
    else  # :ipopt
        set_optimizer(model, Ipopt.Optimizer)
        set_silent(model)
    end
    
    # Decision variables for time step t only
    @variable(model, P_Subs_t >= 0)
    @variable(model, -bat.P_B_R_pu <= P_B_t <= bat.P_B_R_pu)
    @variable(model, B_t)
    
    # SOC box constraints as NEGATIVE inequalities (≤ 0 form for KKT)
    # g_Bmin(B) = Bmin - B ≤ 0  →  B ≥ Bmin
    # g_Bmax(B) = B - Bmax ≤ 0  →  B ≤ Bmax
    @constraint(model, soc_lower, Bmin(bat) - B_t <= 0)
    @constraint(model, soc_upper, B_t - Bmax(bat) <= 0)
    
    # Get values from previous time step (or initial condition)
    B_prev = (t == 1) ? bat.B0_pu : B_collection[t-1]
    
    # Objective function components (keep in per-unit, don't scale yet)
    # 1. Energy cost at time t
    energy_cost = price[t] * P_Subs_t * Δt
    
    # 2. Battery quadratic cost at time t  
    battery_cost = C_B * P_B_t^2 * Δt
    
    # 3. Coupling term with next time step (if exists)
    if t < T
        # Get values from next time step (from last forward pass iteration)
        B_next = B_collection[t+1]
        mu_next = mu_collection[t+1]
        PB_next = PB_collection[t+1]
        
        # Coupling: +μ[t+1] * (B[t+1] - B[t] + Δt*P_B[t+1])
        # This enforces consistency with next step's dynamics
        coupling_term = mu_next * (B_next - B_t + Δt * PB_next)
        
        @objective(model, Min, energy_cost + battery_cost + coupling_term)
    else
        # For t=T (terminal step), no coupling with future
        @objective(model, Min, energy_cost + battery_cost)
    end
    
    # Constraints
    
    # 1. Real power balance at time t
    @constraint(model, power_balance, P_Subs_t + P_B_t == P_L_pu[t])
    
    # 2. Battery dynamics: B[t] = B[t-1] - Δt*P_B[t]
    @constraint(model, battery_dynamics, B_t - B_prev + Δt * P_B_t == 0)
    
    # Solve
    optimize!(model)
    
    # Extract results
    status = termination_status(model)
    
    if has_values(model)
        # Update collections (modified in-place)
        B_collection[t] = value(B_t)
        PB_collection[t] = value(P_B_t)
        
        # Extract dual variable μ[t] from battery dynamics constraint
        # This represents the marginal cost of energy at time t
        # Negate to match DDP formulation conventions
        mu_collection[t] = -dual(battery_dynamics)
        
        # Extract dual variables for SOC box constraints (negative inequality form)
        # g_Bmin: Bmin - B ≤ 0  (i.e., B ≥ Bmin), want λ_Bmin ≥ 0 when active
        # g_Bmax: B - Bmax ≤ 0  (i.e., B ≤ Bmax), want λ_Bmax ≥ 0 when active
        # JuMP returns: dual(g_Bmin) < 0 when lower bound is active (need to negate!)
        #               dual(g_Bmax) > 0 when upper bound is active (correct sign!)
        lambda_Bmin = -dual(soc_lower)  # Negate to get positive when at lower bound
        lambda_Bmax = dual(soc_upper)    # Already correct sign
        
        return Dict(
            :status => status,
            :objective => objective_value(model),
            :P_Subs => value(P_Subs_t),
            :P_B => value(P_B_t),
            :B => value(B_t),
            :mu => -dual(battery_dynamics),
            :lambda_Bmax => lambda_Bmax,
            :lambda_Bmin => lambda_Bmin,
            :solve_time => solve_time(model)
        )
    else
        error("DDP forward pass failed for time step t=$t with status: $status")
    end
end

"""
    solve_CopperPlate_DDP(data; max_iter=100, tol=1e-5, solver=:gurobi)

Solve copper plate MPOPF using DDP (Distributed Dynamic Programming) with forward pass decomposition.

Algorithm:
1. Initialize: B[t], μ[t], P_B[t] for all t
2. For each iteration k:
   - Forward pass: Solve each time step t=1,2,...,T sequentially
   - Each step uses μ[t+1], B[t+1], P_B[t+1] from previous iteration
   - Automatically updates collections
3. Check convergence: ||B^k - B^{k-1}|| < tol

Returns: Dict with solution, convergence history, and timing statistics
"""
function solve_CopperPlate_DDP(data::CopperPlateData; 
                                max_iter::Int=100, 
                                tol::Float64=1e-5,
                                solver::Symbol=:gurobi)
    
    @unpack T, price, P_L_pu, bat, P_BASE, E_BASE, C_B = data
    Δt = bat.Δt
    
    # Initialize collections
    B_collection = Dict{Int, Float64}()
    mu_collection = Dict{Int, Float64}()
    PB_collection = Dict{Int, Float64}()
    PSubs_collection = Dict{Int, Float64}()
    
    # For k=1, initialize all states with B0
    for t in 1:T
        B_collection[t] = bat.B0_pu
        mu_collection[t] = 0.0
        PB_collection[t] = 0.0
        PSubs_collection[t] = 0.0
    end
    
    # Lambda collections for SOC box constraint duals
    lambda_Bmin_collection = Dict{Int, Float64}()
    lambda_Bmax_collection = Dict{Int, Float64}()
    for t in 1:T
        lambda_Bmin_collection[t] = 0.0
        lambda_Bmax_collection[t] = 0.0
    end
    
    # History tracking
    B_history = Vector{Vector{Float64}}()
    mu_history = Vector{Vector{Float64}}()
    PB_history = Vector{Vector{Float64}}()
    lambda_Bmin_history = Vector{Vector{Float64}}()
    lambda_Bmax_history = Vector{Vector{Float64}}()
    obj_history = Float64[]
    convergence_error = Float64[]
    iteration_times = Float64[]
    
    # Store initial state
    push!(B_history, [B_collection[t] for t in 1:T])
    push!(mu_history, [mu_collection[t] for t in 1:T])
    push!(PB_history, [PB_collection[t] for t in 1:T])
    push!(lambda_Bmin_history, [lambda_Bmin_collection[t] for t in 1:T])
    push!(lambda_Bmax_history, [lambda_Bmax_collection[t] for t in 1:T])
    
    println("\n" * "="^80)
    print(COLOR_INFO)
    @printf "🎯 DDP Forward Pass: T=%d, max_iter=%d, tol=%.1e\n" T max_iter tol
    print(COLOR_RESET)
    println("="^80)
    
    converged = false
    final_iter = max_iter
    wallclock_start = time()
    
    for k in 1:max_iter
        iter_start = time()
        
        # Store previous iteration for convergence check (dual variables)
        mu_prev = copy([mu_collection[t] for t in 1:T])
        
        # FORWARD PASS: Solve all time steps sequentially
        # Each step uses B[t-1] from the CURRENT forward pass (just solved)
        # and mu[t+1], B[t+1], PB[t+1] from PREVIOUS iteration
        
        for t in 1:T
            result = forward_pass_ddp!(B_collection, mu_collection, PB_collection, data, t; solver=solver)
            PSubs_collection[t] = result[:P_Subs]
            
            # Extract dual variables for SOC box constraints
            lambda_Bmax_collection[t] = result[:lambda_Bmax]
            lambda_Bmin_collection[t] = result[:lambda_Bmin]
            
            # Note: B_collection[t] is now updated with the fresh value
            # The next time step (t+1) will use this fresh B_collection[t] as B[t]
        end
        
        # VERIFY OPTIMALITY CONDITIONS: λ_Bmax[t] - λ_Bmin[t] + μ_SOC[t] - μ_SOC[t+1] = 0
        # Also check individual lambda signs (should be ≥ 0)
        if k <= 3
            println("\n  " * "─"^70)
            println("  📊 OPTIMALITY & DUAL SIGN CHECK (Iteration $k)")
            println("  " * "─"^70)
            for t in 1:T
                B_t = B_collection[t]
                mu_t = mu_collection[t]
                mu_tp1 = (t < T) ? mu_collection[t+1] : 0.0
                lambda_Bmax_t = lambda_Bmax_collection[t]
                lambda_Bmin_t = lambda_Bmin_collection[t]
                
                # Check if constraints are active
                Bmin_val = Bmin(data.bat)
                Bmax_val = Bmax(data.bat)
                at_lower = abs(B_t - Bmin_val) < 1e-4
                at_upper = abs(B_t - Bmax_val) < 1e-4
                
                opt_residual = lambda_Bmax_t - lambda_Bmin_t + mu_t - mu_tp1
                
                println("  t=$t: B=$(round(B_t, digits=4)) ∈ [$(round(Bmin_val, digits=2)), $(round(Bmax_val, digits=2))]")
                @printf "      λ_Bmin=%+.6f %s,  λ_Bmax=%+.6f %s\n" lambda_Bmin_t (lambda_Bmin_t >= -1e-6 ? "✓" : "⚠") lambda_Bmax_t (lambda_Bmax_t >= -1e-6 ? "✓" : "⚠")
                @printf "      μ[t]=%+.4f,  μ[t+1]=%+.4f  →  Δμ=%+.4f\n" mu_t mu_tp1 (mu_t - mu_tp1)
                @printf "      Optimality residual: %.2e %s\n" opt_residual (abs(opt_residual) < 1e-3 ? "✓" : "⚠")
                if at_lower
                    println("      📍 At LOWER bound (expect λ_Bmin > 0)")
                elseif at_upper
                    println("      📍 At UPPER bound (expect λ_Bmax > 0)")
                end
                println()
            end
            println("  " * "─"^70)
        end
        
        # Compute total objective (energy + battery costs) - scale to physical units at end
        total_energy_cost_pu = sum(price[t] * PSubs_collection[t] * Δt for t in 1:T)
        total_battery_cost_pu = sum(C_B * PB_collection[t]^2 * Δt for t in 1:T)
        total_obj_pu = total_energy_cost_pu + total_battery_cost_pu
        
        # Scale to physical units for reporting
        total_energy_cost = total_energy_cost_pu * P_BASE
        total_battery_cost = total_battery_cost_pu * P_BASE
        total_obj = total_obj_pu * P_BASE
        
        push!(obj_history, total_obj)
        
        # Store iteration history
        push!(B_history, [B_collection[t] for t in 1:T])
        push!(mu_history, [mu_collection[t] for t in 1:T])
        push!(PB_history, [PB_collection[t] for t in 1:T])
        push!(lambda_Bmin_history, [lambda_Bmin_collection[t] for t in 1:T])
        push!(lambda_Bmax_history, [lambda_Bmax_collection[t] for t in 1:T])
        
        # Compute convergence metrics based on DUAL variables (μ)
        mu_curr = [mu_collection[t] for t in 1:T]
        mu_error = norm(mu_curr - mu_prev)
        
        push!(convergence_error, mu_error)
        
        iter_time = time() - iter_start
        push!(iteration_times, iter_time)
        
        # Print progress with KKT condition check
        @printf "k=%3d  obj=\$%.4f  ||Δμ||=%.2e\n" k total_obj mu_error
        
        # Show KKT condition: λ_Bmax - λ_Bmin + μ[t] - μ[t+1] = 0
        # For t < T: four terms
        # For t = T: only three terms (no future μ)
        for t in 1:T
            mu_t = mu_collection[t]
            lambda_Bmin_t = lambda_Bmin_collection[t]
            lambda_Bmax_t = lambda_Bmax_collection[t]
            
            if t < T
                # Interior time steps: λ_Bmax - λ_Bmin + μ[t] - μ[t+1] = 0
                mu_tp1 = mu_collection[t+1]
                residual = lambda_Bmax_t - lambda_Bmin_t + mu_t - mu_tp1
                kkt_satisfied = abs(residual) < 1e-3
                
                @printf "  t=%d: [%+.3f - %+.3f + %+.3f - %+.3f] = %+.2e " t lambda_Bmax_t lambda_Bmin_t mu_t mu_tp1 residual
            else
                # Terminal time step: λ_Bmax - λ_Bmin + μ[T] = 0 (no future coupling)
                residual = lambda_Bmax_t - lambda_Bmin_t + mu_t
                kkt_satisfied = abs(residual) < 1e-3
                
                @printf "  t=%d: [%+.3f - %+.3f + %+.3f] = %+.2e " t lambda_Bmax_t lambda_Bmin_t mu_t residual
            end
            
            if kkt_satisfied
                print(Crayon(foreground=:green, bold=true))
                println("✓")
            else
                print(Crayon(foreground=:red, bold=true))
                println("✗")
            end
            print(COLOR_RESET)
        end
        
        # Check convergence based on dual variables (μ)
        if mu_error < tol
            @printf "🎉 DDP converged at iteration %d (||Δμ||=%.2e < %.2e)\n" k mu_error tol
            converged = true
            final_iter = k
            break
        end
    end
    
    wallclock_time = time() - wallclock_start
    
    if !converged
        print(COLOR_WARNING)
        @printf "⚠ DDP did not converge after %d iterations (||Δμ||=%.2e >= %.2e)\n" max_iter last(convergence_error) tol
        print(COLOR_RESET)
    end
    
    # Compute final objective
    final_obj = last(obj_history)
    final_energy_cost = sum(price[t] * (PSubs_collection[t] * P_BASE) * Δt for t in 1:T)
    final_battery_cost = sum(C_B * (PB_collection[t] * P_BASE)^2 * Δt for t in 1:T)
    
    # Prepare result dictionary
    result = Dict(
        :status => converged ? MOI.OPTIMAL : MOI.ITERATION_LIMIT,
        :objective => final_obj,
        :energy_cost => final_energy_cost,
        :battery_cost => final_battery_cost,
        :converged => converged,
        :iterations => final_iter,
        :solve_time => wallclock_time,
        :wallclock_time => wallclock_time,
        :P_Subs => [PSubs_collection[t] for t in 1:T],
        :P_B => [PB_collection[t] for t in 1:T],
        :B => [B_collection[t] for t in 1:T],
        :mu => [mu_collection[t] for t in 1:T],
        :convergence_history => Dict(
            :B_history => B_history,
            :mu_history => mu_history,
            :PB_history => PB_history,
            :obj_history => obj_history,
            :convergence_error => convergence_error,
            :iteration_times => iteration_times
        )
    )
    
    return result
end

# ============================================================================
# SCENARIO CONFIGURATION
# ============================================================================

begin # scenario config
    # System parameters
    T = 4  # Number of time periods (start small for testing)
    
    # Solver selection
    use_gurobi_for_bf = true   # Use Gurobi for brute force
    use_gurobi_for_ddp = true  # Use Gurobi for DDP subproblems
    
    # DDP algorithm parameters
    max_iter_ddp = 5
    tol_ddp = 1e-5
    
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

# Report results (simplified - just status and objective)
if sol_bf[:status] == MOI.OPTIMAL || sol_bf[:status] == MOI.LOCALLY_SOLVED
    print(COLOR_SUCCESS)
    @printf "✓ Brute Force OPTIMAL: \$%.4f (%.4fs)\n" sol_bf[:objective] sol_bf[:wallclock_time]
    print(COLOR_RESET)
else
    print(COLOR_ERROR)
    @printf "⚠ Brute Force FAILED: %s\n" sol_bf[:status]
    print(COLOR_RESET)
end

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

# ============================================================================
# SOLVE USING DDP (FORWARD PASS DECOMPOSITION)
# ============================================================================

println("\n" * "="^80)
println(COLOR_HIGHLIGHT, "SOLVING COPPER PLATE MPOPF (DDP)", COLOR_RESET)
println("="^80)

solver_ddp_choice = use_gurobi_for_ddp ? :gurobi : :ipopt
sol_ddp = solve_CopperPlate_DDP(data; 
                                 max_iter=max_iter_ddp, 
                                 tol=tol_ddp,
                                 solver=solver_ddp_choice)

# Report results
println("\n--- SOLUTION STATUS ---")
print("Status: ")

if sol_ddp[:status] == MOI.OPTIMAL || sol_ddp[:converged]
    print(COLOR_SUCCESS)
    println(sol_ddp[:status])
    println("✓ DDP converged successfully!")
    print(COLOR_RESET)
    
    println("\n--- OBJECTIVE VALUE ---")
    @printf "Total Cost:     \$%.4f\n" sol_ddp[:objective]
    @printf "  Energy Cost:  \$%.4f (%.1f%%)\n" sol_ddp[:energy_cost] (sol_ddp[:energy_cost]/sol_ddp[:objective]*100)
    @printf "  Battery Cost: \$%.4f (%.1f%%)\n" sol_ddp[:battery_cost] (sol_ddp[:battery_cost]/sol_ddp[:objective]*100)
    
    # Compare with brute force
    if sol_bf[:status] == MOI.OPTIMAL || sol_bf[:status] == MOI.LOCALLY_SOLVED
        obj_diff = abs(sol_ddp[:objective] - sol_bf[:objective])
        obj_rel_error = obj_diff / sol_bf[:objective] * 100
        println("\n--- COMPARISON WITH BRUTE FORCE ---")
        @printf "BF Objective:   \$%.4f\n" sol_bf[:objective]
        @printf "DDP Objective:  \$%.4f\n" sol_ddp[:objective]
        @printf "Difference:     \$%.4f (%.2e%%)\n" obj_diff obj_rel_error
        
        if obj_rel_error < 0.1
            print(COLOR_SUCCESS)
            println("✓ DDP matches brute force solution!")
            print(COLOR_RESET)
        elseif obj_rel_error < 1.0
            print(COLOR_WARNING)
            println("⚠ DDP has small deviation from brute force")
            print(COLOR_RESET)
        else
            print(COLOR_WARNING)
            println("⚠ DDP has significant deviation from brute force")
            print(COLOR_RESET)
        end
    end
    
    println("\n--- COMPUTATION TIME ---")
    @printf "DDP time:       %.4f seconds\n" sol_ddp[:solve_time]
    @printf "Iterations:     %d\n" sol_ddp[:iterations]
    @printf "Time/iter:      %.4f seconds\n" (sol_ddp[:solve_time] / sol_ddp[:iterations])
    
    if sol_bf[:status] == MOI.OPTIMAL || sol_bf[:status] == MOI.LOCALLY_SOLVED
        speedup = sol_bf[:solve_time] / sol_ddp[:solve_time]
        @printf "\nSpeedup vs BF:  %.2fx" speedup
        if speedup > 1.0
            println(" (DDP faster)")
        else
            println(" (BF faster)")
        end
    end
    
    # Convert to physical units
    P_Subs_kW = sol_ddp[:P_Subs] .* data.P_BASE
    P_B_kW = sol_ddp[:P_B] .* data.P_BASE
    B_kWh = sol_ddp[:B] .* data.E_BASE
    
    println("\n--- POWER SUMMARY ---")
    @printf "Substation (kW):  min=%.1f, max=%.1f, avg=%.1f\n" minimum(P_Subs_kW) maximum(P_Subs_kW) mean(P_Subs_kW)
    @printf "Battery (kW):     min=%.1f, max=%.1f, avg=%.1f\n" minimum(P_B_kW) maximum(P_B_kW) mean(P_B_kW)
    
    println("\n--- SOC SUMMARY ---")
    soc_percent = B_kWh ./ E_Rated_kWh .* 100
    @printf "SOC (%%):          min=%.1f, max=%.1f, avg=%.1f\n" minimum(soc_percent) maximum(soc_percent) mean(soc_percent)
    @printf "SOC (kWh):        min=%.1f, max=%.1f, avg=%.1f\n" minimum(B_kWh) maximum(B_kWh) mean(B_kWh)
    
else
    print(COLOR_WARNING)
    println("⚠ DDP did not converge to optimality")
    println("Status: ", sol_ddp[:status])
    println("Iterations: ", sol_ddp[:iterations])
    print(COLOR_RESET)
end

println("\n" * "="^80)
println(COLOR_HIGHLIGHT, "COPPER PLATE MPOPF DDP SOLUTION COMPLETE", COLOR_RESET)
println("="^80)

# ============================================================================
# SAVE DDP RESULTS TO FILE
# ============================================================================

if sol_ddp[:converged]
    processedData_dir = joinpath(@__DIR__, "envs", "ddp", "processedData")
    system_folder = "copperplate_T$(T)"
    system_dir = joinpath(processedData_dir, system_folder)
    mkpath(system_dir)
    
    results_file = joinpath(system_dir, "results_ddp.txt")
    
    open(results_file, "w") do io
        println(io, "="^80)
        println(io, "COPPER PLATE MPOPF - DDP SOLUTION")
        println(io, "="^80)
        println(io)
        
        println(io, "Problem Configuration:")
        println(io, "  Time periods (T):         ", T)
        println(io, "  Solver:                   ", solver_ddp_choice)
        println(io, "  Optimization status:      ", sol_ddp[:status])
        println(io, "  Converged:                ", sol_ddp[:converged])
        println(io, "  Iterations:               ", sol_ddp[:iterations])
        println(io)
        
        println(io, "Objective Value:")
        @printf(io, "  Total Cost:     \$%.6f\n", sol_ddp[:objective])
        @printf(io, "  Energy Cost:    \$%.6f (%.2f%%)\n", sol_ddp[:energy_cost], (sol_ddp[:energy_cost]/sol_ddp[:objective]*100))
        @printf(io, "  Battery Cost:   \$%.6f (%.2f%%)\n", sol_ddp[:battery_cost], (sol_ddp[:battery_cost]/sol_ddp[:objective]*100))
        println(io)
        
        if sol_bf[:status] == MOI.OPTIMAL || sol_bf[:status] == MOI.LOCALLY_SOLVED
            obj_diff = abs(sol_ddp[:objective] - sol_bf[:objective])
            obj_rel_error = obj_diff / sol_bf[:objective] * 100
            println(io, "Comparison with Brute Force:")
            @printf(io, "  BF Objective:   \$%.6f\n", sol_bf[:objective])
            @printf(io, "  DDP Objective:  \$%.6f\n", sol_ddp[:objective])
            @printf(io, "  Difference:     \$%.6f (%.4e%%)\n", obj_diff, obj_rel_error)
            println(io)
        end
        
        println(io, "Computation Time:")
        @printf(io, "  Total time:     %.6f seconds\n", sol_ddp[:solve_time])
        @printf(io, "  Iterations:     %d\n", sol_ddp[:iterations])
        @printf(io, "  Time/iter:      %.6f seconds\n", sol_ddp[:solve_time] / sol_ddp[:iterations])
        println(io)
        
        P_Subs_kW = sol_ddp[:P_Subs] .* data.P_BASE
        P_B_kW = sol_ddp[:P_B] .* data.P_BASE
        B_kWh = sol_ddp[:B] .* data.E_BASE
        
        println(io, "Solution Summary:")
        @printf(io, "  Substation (kW):  min=%.2f, max=%.2f, avg=%.2f\n", minimum(P_Subs_kW), maximum(P_Subs_kW), mean(P_Subs_kW))
        @printf(io, "  Battery (kW):     min=%.2f, max=%.2f, avg=%.2f\n", minimum(P_B_kW), maximum(P_B_kW), mean(P_B_kW))
        @printf(io, "  SOC (kWh):        min=%.2f, max=%.2f, avg=%.2f\n", minimum(B_kWh), maximum(B_kWh), mean(B_kWh))
        println(io)
        
        println(io, "="^80)
    end
    
    println(COLOR_SUCCESS, "✓ DDP results written to $(results_file)", COLOR_RESET)
end

# ============================================================================
# GENERATE PLOTS
# ============================================================================

println("\n" * "="^80)
println(COLOR_INFO, "GENERATING PLOTS", COLOR_RESET)
println("="^80)

processedData_dir = joinpath(@__DIR__, "envs", "ddp", "processedData")
system_folder = "copperplate_T$(T)"
system_dir = joinpath(processedData_dir, system_folder)
mkpath(system_dir)

# Plot battery actions for both solutions
if sol_bf[:status] == MOI.OPTIMAL || sol_bf[:status] == MOI.LOCALLY_SOLVED
    println("\n📊 Plotting battery actions (Brute Force)...")
    plot_battery_actions(sol_bf, data, "Brute Force",
                        showPlots=showPlots,
                        savePlots=savePlots,
                        filename=joinpath(system_dir, "battery_actions_bf.png"))
end

if sol_ddp[:converged]
    println("📊 Plotting battery actions (DDP)...")
    plot_battery_actions(sol_ddp, data, "DDP",
                        showPlots=showPlots,
                        savePlots=savePlots,
                        filename=joinpath(system_dir, "battery_actions_ddp.png"))
    
    println("📊 Plotting DDP convergence...")
    plot_ddp_convergence(sol_ddp,
                        showPlots=showPlots,
                        savePlots=savePlots,
                        filename=joinpath(system_dir, "ddp_convergence.png"))
    
    if sol_bf[:status] == MOI.OPTIMAL || sol_bf[:status] == MOI.LOCALLY_SOLVED
        println("📊 Plotting comparison...")
        plot_comparison(sol_bf, sol_ddp, data,
                       showPlots=showPlots,
                       savePlots=savePlots,
                       filename=joinpath(system_dir, "comparison_bf_vs_ddp.png"))
    end
end

println("\n" * "="^80)
println(COLOR_SUCCESS, "PLOTTING COMPLETE", COLOR_RESET)
println("="^80)
println("📁 Results and plots saved to: ", system_dir)

end # entire script
