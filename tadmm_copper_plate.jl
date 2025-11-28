# ============================================================================
# PARALLELIZATION CONFIGURATION
# ============================================================================
# Auto-detect available CPU cores - try multiple methods for reliability
function detect_cpu_cores()
    # Try environment variable first (most reliable if set)
    if haskey(ENV, "NUMBER_OF_PROCESSORS")
        return parse(Int, ENV["NUMBER_OF_PROCESSORS"])
    end
    # Try Sys.CPU_THREADS (works on most systems)
    sys_threads = Sys.CPU_THREADS
    if sys_threads > 1
        return sys_threads
    end
    # Fallback: check what Julia is actually using
    return max(1, Threads.nthreads())
end

const AVAILABLE_CORES = detect_cpu_cores()
const USE_THREADING = true  # Enable multi-threading for time-period subproblems
# ============================================================================

# Activate the tadmm environment
import Pkg
env_path = joinpath(@__DIR__, "envs", "tadmm")
Pkg.activate(env_path)

using Revise
using JuMP
using Ipopt
using LinearAlgebra
using Printf
using Gurobi
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
    println("\nTo use all cores, restart Julia with:")
    println("  \$ export JULIA_NUM_THREADS=$(AVAILABLE_CORES)")
    println("  \$ julia tadmm_copper_plate.jl")
else
    println("Status:              ⚠  Single-threaded mode")
    println("\nTo enable parallel execution, restart Julia with:")
    println("  \$ export JULIA_NUM_THREADS=$(AVAILABLE_CORES)")
    println("  \$ julia tadmm_copper_plate.jl")
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

# ----------------------- Bases ---------------------------
systemName = "cp1"  # Copper Plate 1
# max_iter = 10000
max_iter = 1000
rho = 10.0                     # ADMM penalty parameter (initial value if adaptive)
eps_pri = 1e-5
eps_dual = 1e-4
adaptive_rho = true            # Enable adaptive ρ adjustment
# ----------------------- Scenario ---------------------------
# T = 24
T = 24
LoadShapeCost = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T)) .+ 1) ./ 2
# LoadShapeCost = 0.08*ones(T)
# LoadShapeCost = [mod(t-1, 2) == 0 ? 0.05 : 0.15 for t in 1:T]  # Square wave: 0.05 $/kWh (low) and 0.15 $/kWh (high)
# Battery quadratic cost coefficient
# C_B = 0
C_B = 1e-6 * minimum(LoadShapeCost)  # Quadratic cost coefficient for battery power: C_B * P_B^2
LoadShapeLoad = 0.8 .+ 0.2 .* (sin.(range(0, 2π, length=T) .- 0.8) .+ 1) ./ 2  # Normalized load shape [0, 1]
# LoadShapeLoad = 1 * ones(T)  # Normalized load shape [0, 1]
const KV_B = 4.16 / sqrt(3)   # kV (unused here)
const KVA_B = 1000.0           # kVA
const P_BASE = KVA_B            # kW
const E_BASE = P_BASE * 1.0     # kWh per 1 hour
# ----------------------- Data ----------------------------

struct BatteryParamsPU
    E_Rated_pu::Float64          # puh on E_BASE
    soc_min::Float64             # pu
    soc_max::Float64             # pu
    P_B_R_pu::Float64            # pu on P_BASE
    Δt::Float64                  # hours
    B0_pu::Float64               # initial energy (puh)
    B_T_target_pu::Union{Nothing,Float64}  # terminal (puh), optional
end

struct InstancePU
    T::Int
    price::Vector{Float64}       # $/kWh
    P_L_pu::Vector{Float64}      # pu on P_BASE
    bat::BatteryParamsPU
end

Bmin(b::BatteryParamsPU) = b.soc_min * b.E_Rated_pu
Bmax(b::BatteryParamsPU) = b.soc_max * b.E_Rated_pu

"""
build_instance_pu(...)
Converts physical inputs (kW/kWh) to per-unit given the fixed bases.
"""
function build_instance_pu(T::Int, price::Vector{Float64}, P_L_kW::Vector{Float64},
    E_Rated_kWh::Float64, soc_min::Float64, soc_max::Float64,
    P_B_R_kW::Float64, Δt::Float64, B0_kWh::Float64;
    B_T_target_kWh::Union{Nothing,Float64}=nothing)
    @assert length(price) == T && length(P_L_kW) == T
    P_L_pu = P_L_kW ./ P_BASE
    E_Rated_pu = E_Rated_kWh / E_BASE
    P_B_R_pu = P_B_R_kW / P_BASE
    B0_pu = B0_kWh / E_BASE
    B_T_target_pu = isnothing(B_T_target_kWh) ? nothing : (B_T_target_kWh / E_BASE)
    bat = BatteryParamsPU(E_Rated_pu, soc_min, soc_max, P_B_R_pu, Δt, B0_pu, B_T_target_pu)
    return InstancePU(T, price, P_L_pu, bat)
end

# ------------------- TADMM Update Functions --------------

"""
    primal_update_tadmm!(B_t0, Bhat, u_t0, inst, ρ, t0)

🔵 PRIMAL UPDATE for TADMM Subproblem t0 🔵

Solves the t0-th subproblem in TADMM formulation:
    min C_t0 * P_subs_t0 + (ρ/2) * ||🔵B_t0 - 🔴B̂ + 🟢u_t0||²₂

Variables being optimized (for subproblem t0):
- 🔵B_t0[t0]: Local SOC decision variable at time t0 (1 scalar)
- P_B_t0: Battery power at time t0 (1 scalar) 
- P_subs_t0: Substation power at time t0 (1 scalar)

Fixed parameters from previous iteration:
- 🔴B̂: Global consensus SOC trajectory (T-length vector)
- 🟢u_t0: Local scaled dual variable (T-length vector)

Constraints in this subproblem:
📊 CONSTRAINT COUNT:
    ✅ Equality constraints: 2
    1. SOC trajectory: 🔵B_t0[t0] = 🔴B̂[t0-1] - P_B_t0 * Δt  (1 constraint)
    2. NRPB: P_subs_t0 + P_B_t0 = P_L[t0]  (1 constraint)

    🚧 Inequality constraints: T + 1  
    1. SOC bounds: 🔵B_t0[t] ∈ [B_lower, B_upper] ∀t ∈ {1,...,T}  (T constraints)
    2. Battery power bounds: P_B_t0 ∈ [-P_B_R, P_B_R]  (1 constraint)

Arguments:
- B_t0: Local SOC variables for subproblem t0 (🔵B_t0 - modified in-place)
- Bhat: Global consensus SOC (🔴B̂ - read-only)  
- u_t0: Local scaled dual variables for subproblem t0 (🟢u_t0 - read-only)
- inst: Problem instance
- ρ: Penalty parameter
- t0: Time index for this subproblem

Returns: Dict with keys:
    - :objective => objective value for this subproblem
    - :P_B => battery power dispatch P_B_t0
    - :P_subs => substation power dispatch P_subs_t0  
    - :B_t0 => updated local SOC vector (🔵B_t0)
    - :t0 => time index of this subproblem (for reference)
"""
function primal_update_tadmm!(B_t0, Bhat, u_t0, inst::InstancePU, ρ::Float64, t0::Int)
    b = inst.bat
    
    # 🎯 Setup optimization model for subproblem t0
    # m = Model(Ipopt.Optimizer)
    m = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
    set_silent(m)
    set_optimizer_attribute(m, "OutputFlag", 0)
    set_optimizer_attribute(m, "Threads", 1)  # Each subproblem uses 1 thread
    set_optimizer_attribute(m, "DualReductions", 0)  # Better infeasibility diagnosis

    # 🔵 Decision variables for time t0 (INCREMENTAL: P_B as full vector)
    @variables(m, begin
        # Battery power for ALL time steps (T-length vector with box constraints)
        -b.P_B_R_pu <= P_B_var[1:inst.T] <= b.P_B_R_pu         
        # Substation power at time t0 (scalar)  
        P_subs_t0                                     
        # 🔵 Local SOC trajectory (T-length vector, but only B_t0[t0] is truly optimized)
        Bmin(b) <= B_t0_var[1:inst.T] <= Bmax(b)                  
    end)

    # 📌 CONSTRAINT 1: SOC Dynamics for ALL time steps (T equality constraints)
    # B[1] = B0 - P_B[1] * Δt
    # B[t] = B[t-1] - P_B[t] * Δt  for t = 2, ..., T
    @constraint(m, B_t0_var[1] == b.B0_pu - P_B_var[1] * b.Δt)
    for t in 2:inst.T
        @constraint(m, B_t0_var[t] == B_t0_var[t-1] - P_B_var[t] * b.Δt)
    end
    
    # 📌 CONSTRAINT 2: Nodal Real Power Balance ONLY for time t0 (1 equality constraint)
    # P_subs_t0 + P_B[t0] = P_L[t0]
    @constraint(m, P_subs_t0 + P_B_var[t0] == inst.P_L_pu[t0])

    # 🎯 OBJECTIVE: Economic cost + Battery quadratic cost + ADMM penalty
    # C_t0 * P_subs_t0 + C_B * P_B[t0]^2 + (ρ/2) * ||🔵B_t0 - 🔴B̂ + 🟢u_t0||²₂
    energy_cost = inst.price[t0] * (P_subs_t0 * P_BASE) * b.Δt
    battery_quad_cost = C_B * (P_B_var[t0] * P_BASE)^2 * b.Δt  # Quadratic cost for t0
    
    # ADMM penalty: (ρ/2) * ||🔵B_t0 - 🔴B̂ + 🟢u_t0||²₂
    penalty = (ρ / 2) * sum((B_t0_var[t] - Bhat[t] + u_t0[t])^2 for t in 1:inst.T)
    
    @objective(m, Min, energy_cost + battery_quad_cost + penalty)

    # 🚀 Solve subproblem
    optimize!(m)

    # 📥 Extract results and update 🔵B_t0
    # B_t0[t0] = value(B_t0_var[t0])  # Only update the t0-th component
    for t in 1:inst.T
        B_t0[t] = value(B_t0_var[t])  # Update ALL components
    end
    P_B_val = value(P_B_var[t0])    # Extract t0-th P_B value
    P_subs_val = value(P_subs_t0)
    
    #  Compute objective components using solved values (not symbolic expressions)
    energy_cost_val = inst.price[t0] * (P_subs_val * P_BASE) * b.Δt
    battery_quad_cost_val = C_B * (P_B_val * P_BASE)^2 * b.Δt
    penalty_val = (ρ / 2) * sum((value(B_t0_var[t]) - Bhat[t] + u_t0[t])^2 for t in 1:inst.T)
    
    # 📦 Return results as extensible dictionary with objective breakdown
    return Dict(
        :total_objective => objective_value(m),
        :energy_cost => energy_cost_val,
        :battery_quad_cost => battery_quad_cost_val, 
        :penalty => penalty_val,
        :P_B => P_B_val,
        :P_subs => P_subs_val,
        :B_t0 => B_t0,
        :t0 => t0,
        # Easy to add more fields later: :solver_status, :solve_time, etc.
    )
end

"""
    consensus_update_tadmm!(Bhat, B_collection, u_collection, inst, ρ)

🔴 CONSENSUS UPDATE for TADMM 🔴  

Updates global consensus variables 🔴B̂ using averaging of local solutions:
    🔴B̂[t] = (1/T) * Σ_{t0=1}^T (🔵B_t0[t] + 🟢u_t0[t])

Arguments:
- Bhat: Global consensus SOC trajectory (🔴B̂ - modified in-place)
- B_collection: Collection of all local SOC variables {🔵B_t0} for t0=1:T
- u_collection: Collection of all local scaled duals {🟢u_t0} for t0=1:T  
- inst: Problem instance
- ρ: Penalty parameter (unused in consensus update but kept for interface consistency)

Returns: Dict with keys:
    - :Bhat => updated global consensus trajectory (🔴B̂)
    - :bounds_violations => number of times clamping was needed
    - :violation_indices => time indices where bounds violations occurred
"""
function consensus_update_tadmm!(Bhat, B_collection, u_collection, inst::InstancePU, ρ::Float64)
    T = inst.T
    b = inst.bat
    violations = Int[]
    
    violation_tolerance = 1e-4  # Only warn if violation is > 1e-6

    #  Update consensus variables: B̂[t] for t = 1, 2, ..., T-1  
    # Note: B̂[t] represents SOC at END of time period t
    # B0 (initial SOC) is handled separately in primal updates
    for t in 1:T-1  # Don't update the last time step if it has terminal constraint
        # Average across all T subproblems with dual adjustments
        # 🔴B̂[t] = (1/T) * Σ_{t0=1}^T (🔵B_t0[t] + 🟢u_t0[t])
        consensus_sum = sum(B_collection[t0][t] + u_collection[t0][t] for t0 in 1:T)
        Bhat_new = consensus_sum / T
        
        # Check if projection is needed and track violations
        violation_amount = max(Bmin(b) - Bhat_new, Bhat_new - Bmax(b), 0.0)

        if violation_amount > violation_tolerance
            push!(violations, t)
            @warn "🚨 Consensus B̂[$t] = $(Bhat_new) violates bounds [$(Bmin(b)), $(Bmax(b))] by $(violation_amount). Clamping applied."
        elseif Bhat_new < Bmin(b) || Bhat_new > Bmax(b)
            # Still track minor violations for statistics, but don't warn
            push!(violations, t)
        end
        
        # Project onto feasible set
        Bhat[t] = clamp(Bhat_new, Bmin(b), Bmax(b))
    end
    
    # 📌 Terminal condition (if specified)
    if !isnothing(b.B_T_target_pu)
        Bhat[T] = b.B_T_target_pu  # B̂[T] = B_T_target (given terminal condition)
    else
        # If no terminal constraint, update the last time step too
        t = T
        consensus_sum = sum(B_collection[t0][t] + u_collection[t0][t] for t0 in 1:T)
        Bhat_new = consensus_sum / T
        
        if Bhat_new < Bmin(b) || Bhat_new > Bmax(b)
            violation_amount = max(Bmin(b) - Bhat_new, Bhat_new - Bmax(b), 0.0)

            if violation_amount > violation_tolerance
                push!(violations, t)
                @warn "🚨 Consensus B̂[$t] = $(Bhat_new) violates bounds [$(Bmin(b)), $(Bmax(b))] by $(violation_amount). Clamping applied."
            else
                push!(violations, t)
            end
        end
        
        Bhat[T] = clamp(Bhat_new, Bmin(b), Bmax(b))
    end
    
    # 📦 Return results as extensible dictionary
    return Dict(
        :Bhat => Bhat,
        :bounds_violations => length(violations),
        :violation_indices => violations,
        # Easy to add: :consensus_change, :max_violation_amount, etc.
    )
end

"""
    dual_update_tadmm!(u_collection, B_collection, Bhat, ρ)

🟢 DUAL UPDATE for TADMM 🟢

Updates scaled dual variables for each subproblem:
    🟢u_t0[t] := 🟢u_t0[t] + (🔵B_t0[t] - 🔴B̂[t])

Arguments:
- u_collection: Collection of local scaled dual variables {🟢u_t0} (modified in-place)
- B_collection: Collection of local SOC variables {🔵B_t0} (read-only)
- Bhat: Global consensus SOC trajectory (🔴B̂ - read-only)
- ρ: Penalty parameter (ρ scaling absorbed into u)

Returns: Dict with keys:
    - :u_collection => updated dual variable collection {🟢u_t0}
    - :max_dual_change => maximum absolute change in any dual variable
    - :total_updates => total number of dual variables updated (T²)
"""
function dual_update_tadmm!(u_collection, B_collection, Bhat, ρ::Float64)
    T = length(Bhat)
    max_change = 0.0
    
    # 🟢 Update scaled dual variables for each subproblem t0
    for t0 in 1:T
        for t in 1:T
            # Dual ascent step: 🟢u_t0[t] += (🔵B_t0[t] - 🔴B̂[t])
            old_u = u_collection[t0][t]
            u_collection[t0][t] += (B_collection[t0][t] - Bhat[t])
            
            # Track maximum change for diagnostics
            max_change = max(max_change, abs(u_collection[t0][t] - old_u))
        end
    end
    
    # 📦 Return results as extensible dictionary
    return Dict(
        :u_collection => u_collection,
        :max_dual_change => max_change,
        :total_updates => T * T,
        # Easy to add: :dual_norms, :convergence_metrics, etc.
    )
end

# ------------------- tADMM (T single-step blocks) --------
"""
solve_MPOPF_using_tADMM(inst; ρ=5.0, max_iter=200, eps_pri=1e-3, eps_dual=1e-3)

🎯 tADMM SOLVER using PDF formulation notation 🎯

Variables:
- 🔵B: Local SOC variables {B_t0[t]} for each subproblem t0=1:T
- 🔴B̂: Global consensus SOC trajectory  
- 🟢u: Local scaled dual variables {u_t0[t]} for each subproblem t0=1:T

Algorithm:
1. 🔵 Primal Update: Solve T subproblems in parallel
2. 🔴 Consensus Update: Average local solutions  
3. 🟢 Dual Update: Update scaled dual variables

Returns Dict(:P_B, :P_Subs, :B, :objective_history, :consensus_trajectory, :convergence_history)
"""
function solve_MPOPF_using_tADMM(inst::InstancePU; ρ::Float64=1.0,
    max_iter::Int=200, eps_pri::Float64=1e-3, eps_dual::Float64=1e-3, adaptive_rho::Bool=false)
    T = inst.T
    b = inst.bat

    # 🔴 Initialize global consensus trajectory B̂ with B0
    Bhat = fill(b.B0_pu, T)  # Initialize all time steps with B0
    # Clamp to ensure feasibility
    Bhat .= clamp.(Bhat, Bmin(b), Bmax(b))

    # 🔵 Initialize local SOC variables {B_t0} for each subproblem
    B_collection = [copy(Bhat) for t0 in 1:T]  # T copies of T-length vectors

    # 🟢 Initialize scaled dual variables {u_t0} for each subproblem  
    u_collection = [zeros(T) for t0 in 1:T]  # T copies of T-length vectors

    # 📊 Initialize power collections to store results from each subproblem
    P_B_collection = zeros(T)   # P_B[t] from subproblem t
    P_Subs_collection = zeros(T)  # P_Subs[t] from subproblem t
    B_local_collection = zeros(T)  # B[t] from subproblem t (local blue solutions)

    # 📊 History tracking
    obj_history = Float64[]  # True objective (energy + battery costs only)
    Bhat_history = Vector{Vector{Float64}}()
    B_collection_history = Vector{Vector{Vector{Float64}}}()
    u_collection_history = Vector{Vector{Vector{Float64}}}()
    r_norm_history = Float64[]
    s_norm_history = Float64[]
    rho_history = Float64[]  # Track ρ evolution if adaptive

    # Adaptive ρ parameters (simplified for copper plate - matches tadmm_socp spirit)
    μ_balance = 10.0          # Threshold for imbalance (standard Boyd parameter)
    τ_incr = 2.0              # Factor to increase ρ
    τ_decr = 2.0              # Factor to decrease ρ
    ρ_min = 0.1               # Minimum ρ value
    ρ_max = 1e6               # Maximum ρ value
    update_interval = 10      # Update ρ every N iterations (copper plate converges faster)

    # Store initial states
    push!(Bhat_history, copy(Bhat))
    push!(B_collection_history, deepcopy(B_collection))
    push!(u_collection_history, deepcopy(u_collection))

    @printf "🎯 tADMM[PDF-formulation]: T=%d, ρ_init=%.3f, adaptive=%s\n" T ρ adaptive_rho

    for k in 1:max_iter
        # 🔵 STEP 1: Primal Update - Solve T subproblems
        @printf "  🔵 Primal updates: "
        total_energy_cost = 0.0
        total_battery_cost = 0.0
        total_penalty = 0.0
        for t0 in 1:T
            result = primal_update_tadmm!(B_collection[t0], Bhat, u_collection[t0], inst, ρ, t0)
            # print(Bhat)
            # Ensure B_collection[t0] is updated with the result
            B_collection[t0] = copy(result[:B_t0])
            # Store power results from each subproblem t0
            P_B_collection[t0] = result[:P_B]
            P_Subs_collection[t0] = result[:P_subs]
            # Store local SOC from each subproblem t0 (the actual optimized SOC at time t0)
            B_local_collection[t0] = B_collection[t0][t0]
            
            # Accumulate ONLY the primary objective components (not ADMM penalty)
            total_energy_cost += result[:energy_cost]
            total_battery_cost += result[:battery_quad_cost]
            total_penalty += result[:penalty]  # Track penalty separately
            # @printf "%d " t0
        end
        @printf "\n"
        
        # True objective is energy + battery costs (penalty is just for convergence)
        true_objective = total_energy_cost + total_battery_cost
        push!(obj_history, true_objective)

        # 🔴 STEP 2: Consensus Update  
        Bhat_old = copy(Bhat)
        consensus_result = consensus_update_tadmm!(Bhat, B_collection, u_collection, inst, ρ)

        # 🟢 STEP 3: Dual Update
        dual_result = dual_update_tadmm!(u_collection, B_collection, Bhat, ρ)

        # 📊 Store iteration history
        push!(Bhat_history, copy(Bhat))
        push!(B_collection_history, deepcopy(B_collection))
        push!(u_collection_history, deepcopy(u_collection))
        push!(rho_history, ρ)  # Track current ρ value

        # 📏 STEP 4: Compute residuals using vector norms
        # Primal residual: measure consensus violation using 2-norm
        r_vectors = []
        for t0 in 1:T
            push!(r_vectors, B_collection[t0] - Bhat)
        end
        r_norm = 1/(inst.T) * norm(vcat(r_vectors...))  # Concatenate all residual vectors and take 2-norm

        # Dual residual: measure change in consensus
        s_norm = 1/(inst.T)* ρ * norm(Bhat - Bhat_old)

        push!(r_norm_history, r_norm)
        push!(s_norm_history, s_norm)

        # 🔧 Adaptive ρ adjustment (simplified Boyd et al. 2011 - better for copper plate)
        if adaptive_rho && k % update_interval == 0
            ρ_old = ρ
            
            if r_norm > μ_balance * s_norm
                # Primal residual too large -> INCREASE rho to enforce consensus
                ρ = min(ρ_max, τ_incr * ρ)
                print(COLOR_WARNING)
                @printf "  ⬆ ρ increased: %.2f → %.2f (r/s=%.1f > μ=%.1f)\n" ρ_old ρ (r_norm/s_norm) μ_balance
                print(COLOR_RESET)
                # Rescale dual variables
                for t0 in 1:inst.T
                    u_collection[t0] ./= τ_incr
                end
            elseif s_norm > μ_balance * r_norm
                # Dual residual too large -> DECREASE rho to speed convergence
                ρ = max(ρ_min, ρ / τ_decr)
                print(COLOR_WARNING)
                @printf "  ⬇ ρ decreased: %.2f → %.2f (s/r=%.1f > μ=%.1f)\n" ρ_old ρ (s_norm/r_norm) μ_balance
                print(COLOR_RESET)
                # Rescale dual variables
                for t0 in 1:inst.T
                    u_collection[t0] .*= τ_decr
                end
            end
        end

        @printf "k=%3d  obj=%.4f  ‖r‖=%.2e  ‖s‖=%.2e  ρ=%.2f\n" k true_objective r_norm s_norm ρ

        if r_norm ≤ eps_pri && s_norm ≤ eps_dual
            @printf "🎉 tADMM converged at iteration %d\n" k
            break
        end
    end

    # 📤 Return results using power values from latest subproblem optimizations
    return Dict(
        :P_B => P_B_collection,  # From latest subproblem optimizations
        :P_Subs => P_Subs_collection,  # From latest subproblem optimizations
        :B => B_local_collection,  # Use local SOC solutions (🔵 blue variables)
        :objective => last(obj_history),
        :objective_history => obj_history,
        :consensus_trajectory => Bhat,  # 🔴B̂ final trajectory
        :local_solutions => B_collection,  # All 🔵B_t0 solutions
        :dual_variables => u_collection,  # All 🟢u_t0 variables
        :convergence_history => Dict(
            :Bhat_history => Bhat_history,
            :B_collection_history => B_collection_history,
            :u_collection_history => u_collection_history,
            :r_norm_history => r_norm_history,
            :s_norm_history => s_norm_history,
            :obj_history => obj_history,
            :rho_history => rho_history
        )
    )
end

# ---------------- Brute-force (monolithic, pu) -----------
function solve_MPOPF_using_BruteForce(inst::InstancePU)
    T = inst.T
    b = inst.bat
    m = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
    set_silent(m)
    set_optimizer_attribute(m, "OutputFlag", 0)
    set_optimizer_attribute(m, "DualReductions", 0)  # Better infeasibility diagnosis

    @variables(m, begin
        -b.P_B_R_pu <= P_B[1:T] <= b.P_B_R_pu
        P_Subs[1:T]
        Bmin(b) <= B[1:T] <= Bmax(b)
    end)

    @constraint(m, [t = 1:T], P_Subs[t] + P_B[t] == inst.P_L_pu[t])
    @constraint(m, B[1] == b.B0_pu - P_B[1] * b.Δt)
    @constraint(m, [t = 2:T], B[t] == B[t-1] - P_B[t] * b.Δt)
    if !isnothing(b.B_T_target_pu)
        @constraint(m, B[T] == b.B_T_target_pu)
    end

    @objective(m, Min, sum(inst.price[t] * (P_Subs[t] * P_BASE) * b.Δt + C_B * (P_B[t] * P_BASE)^2 * b.Δt for t in 1:T))
    optimize!(m)

    return Dict(:P_B => value.(P_B), :P_Subs => value.(P_Subs), :B => value.(B),
        :objective => objective_value(m))
end


# ----------------------- Example -------------------------


# Load profile definition
P_L_R_kW = 1250.0  # Rated load power (kW) - maximum possible load
P_L_kW = P_L_R_kW .* LoadShapeLoad  # Actual load profile (kW)

E_Rated_kWh = 4000.0
soc_min, soc_max = 0.30, 0.95
P_B_R_kW = 800.0
Δt_h = 1.0
B0_kWh = 0.60 * E_Rated_kWh
# B_T_target_kWh = 0.60 * E_Rated_kWh

inst = build_instance_pu(T, LoadShapeCost, P_L_kW,
    E_Rated_kWh, soc_min, soc_max,
    P_B_R_kW, Δt_h, B0_kWh;
    B_T_target_kWh=nothing)

println("\n" * "="^80)
println(COLOR_HIGHLIGHT, "SOLVING MPOPF WITH BRUTE FORCE (MONOLITHIC)", COLOR_RESET)
println("="^80)

# Brute-force baseline (unchanged)
sol_bf = solve_MPOPF_using_BruteForce(inst)

println("\n--- BRUTE FORCE RESULTS ---")
print(COLOR_SUCCESS)
@printf "Objective: \$%.4f\n" sol_bf[:objective]
print(COLOR_RESET)
@printf "B[min,max]=[%.4f, %.4f] pu  | P_B[min,max]=[%.4f, %.4f] pu\n" minimum(sol_bf[:B]) maximum(sol_bf[:B]) minimum(sol_bf[:P_B]) maximum(sol_bf[:P_B])
println("="^80)

println("\n" * "="^80)
println(COLOR_HIGHLIGHT, "SOLVING MPOPF WITH tADMM (TEMPORAL DECOMPOSITION)", COLOR_RESET)
println("="^80)

# tADMM with PDF formulation variable names 
sol_tadmm = solve_MPOPF_using_tADMM(inst; ρ=rho, max_iter=max_iter, eps_pri=eps_pri, eps_dual=eps_dual, adaptive_rho=adaptive_rho)

println("\n--- tADMM RESULTS ---")
print(COLOR_SUCCESS)
@printf "Objective: \$%.4f\n" sol_tadmm[:objective]
print(COLOR_RESET)
@printf "B[min,max]=[%.4f, %.4f] pu  | P_B[min,max]=[%.4f, %.4f] pu\n" minimum(sol_tadmm[:B]) maximum(sol_tadmm[:B]) minimum(sol_tadmm[:P_B]) maximum(sol_tadmm[:P_B])

# Convert tADMM results back to physical for sanity checks
P_B_kW = sol_tadmm[:P_B] .* P_BASE
PSubs_kW = sol_tadmm[:P_Subs] .* P_BASE
B_kWh = sol_tadmm[:B] .* E_BASE
@printf "P_B[kW] in [%.1f, %.1f], B[kWh] in [%.1f, %.1f]\n" minimum(P_B_kW) maximum(P_B_kW) minimum(B_kWh) maximum(B_kWh)
println("="^80)

# ----------------------- Plotting Functions --------------
using Plots

"""
    downsample_for_plotting(data_vectors...; max_points::Int=25)

Downsample multiple data vectors to at most max_points for cleaner plotting.
Always includes first and last points, then evenly spaces the rest.
Returns tuple of (downsampled_indices, downsampled_vectors...)
"""
function downsample_for_plotting(data_vectors...; max_points::Int=25)
    n_points = length(data_vectors[1])
    
    # If we have fewer points than max_points, return all data
    if n_points <= max_points
        return (1:n_points, data_vectors...)
    end
    
    # Create indices that include first, last, and evenly spaced points in between
    # Always include iteration 1 and final iteration
    indices = [1]  # Always include first point
    
    if max_points > 2
        # Add evenly spaced intermediate points
        intermediate_count = max_points - 2
        step = (n_points - 1) / (intermediate_count + 1)
        for i in 1:intermediate_count
            idx = round(Int, 1 + i * step)
            if idx != 1 && idx != n_points  # Avoid duplicates
                push!(indices, idx)
            end
        end
    end
    
    push!(indices, n_points)  # Always include last point
    indices = sort(unique(indices))  # Remove any duplicates and sort
    
    # Downsample all data vectors
    downsampled_vectors = tuple([vec[indices] for vec in data_vectors]...)
    
    return (indices, downsampled_vectors...)
end

"""
    plot_load_and_cost_curves(; showPlots::Bool=true, savePlots::Bool=false, filename::String="load_cost_curves.png")

Plot LoadShape and cost curves using dark yellow and dark green themes on the same plot.
Uses the global variables T, LoadShapeLoad, and LoadShapeCost from the current file.
"""
function plot_load_and_cost_curves(; showPlots::Bool=true, savePlots::Bool=false, filename::String="load_cost_curves.png")
    # Prepare data for plotting
    time_steps = 1:T
    load_cost_cents = LoadShapeCost .* 100  # Convert from $/kWh to cents/kWh
    
    # Calculate y-axis limits
    left_min = -0.05
    left_max = 1.05 * maximum(LoadShapeLoad)
    right_min = floor(0.95 * minimum(load_cost_cents))
    right_max = ceil(1.05 * maximum(load_cost_cents))

    # Set theme and backend
    gr()
    theme(:mute)
    
    # Create main plot with LoadShape (dark yellow theme)
    p = plot(
        time_steps, LoadShapeLoad,
        dpi=600,
        label="Loading Factor (λᵗ)",
        xlabel="Time Period (t)",
        ylabel="Loading Factor [Dimensionless]",
        legend=:bottomright,
        lw=4,
        color=:darkgoldenrod2,  # Dark yellow theme
        markershape=:square,
        markersize=6,
        markerstrokecolor=:black,
        markerstrokewidth=2.0,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        ylims=(left_min, left_max),
        xticks=1:T,
        title="Load Shape and Cost Curves",
        titlefont=font(14, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern",
        right_margin=10Plots.mm
    )

    # Add secondary y-axis for cost curve (dark green theme)
    ax2 = twinx()
    plot!(
        ax2, time_steps, load_cost_cents,
        label="Substation Power Cost (Cᵗ)",
        lw=4,
        color=:darkgreen,  # Dark green theme
        linestyle=:solid,
        markershape=:diamond,
        markersize=7,
        markerstrokecolor=:black,
        markerstrokewidth=2.0,
        ylabel="Cost [cents/kWh]",
        ylims=(right_min, right_max),
        legend=:bottomleft,
        titlefont=font(14, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )

    # Show the plot if requested
    if showPlots
        display(p)
    end

    # Save the plot if requested
    if savePlots
        @printf "Saving plot to: %s\n" filename
        savefig(p, filename)
    end
    
    return p
end

"""
    plot_battery_actions_single(solution, inst, method_name; showPlots::Bool=true, savePlots::Bool=false, filename_prefix::String="battery_actions", rho_val::Union{Nothing,Float64}=nothing)

Plot battery charging/discharging power and SOC for a single solution method.
Creates separate charging/discharging bars and SOC plot in the same style as Plotter.jl.
If rho_val is provided and method_name contains "tADMM", includes rho value in the title.
"""
function plot_battery_actions_single(solution, inst::InstancePU, method_name::String; 
    showPlots::Bool=true, savePlots::Bool=false, filename_prefix::String="battery_actions", rho_val::Union{Nothing,Float64}=nothing)
    
    T = inst.T
    b = inst.bat
    time_steps = 1:T
    
    # Convert to physical units for plotting
    P_B_kW = solution[:P_B] .* P_BASE
    B_kWh = [b.B0_pu * E_BASE; solution[:B] .* E_BASE]  # Include initial SOC
    
    # Separate charging and discharging (P_B > 0 is discharging, P_B < 0 is charging)
    charging_power_kW = -min.(P_B_kW, 0.0)  # P_B < 0 made positive for green bars
    discharging_power_kW = max.(P_B_kW, 0.0)  # P_B > 0 kept positive for wine red bars
    
    # Calculate SOC percentages
    E_Rated_kWh = b.E_Rated_pu * E_BASE
    soc_percent = B_kWh ./ E_Rated_kWh .* 100
    
    # Physical limits
    P_B_R_kW = b.P_B_R_pu * P_BASE
    ylimit = (-P_B_R_kW, P_B_R_kW)
    
    # Set theme
    gr()
    theme(:mute)
    
    # Create title with rho value if applicable
    title_text = if !isnothing(rho_val) && occursin("tADMM", method_name)
        "Battery Actions - $(method_name) (ρ=$(rho_val))\nCharging and Discharging"
    else
        "Battery Actions - $(method_name)\nCharging and Discharging"
    end
    
    # Create charging/discharging bar plot with exact positioning
    # P_B bars: positioned at interval boundaries [0.5, 1.5, 2.5, 3.5] for intervals [0,1), [1,2), [2,3), [3,4)
    power_x_positions = collect(0.5:1.0:(T-0.5))  # [0.5, 1.5, 2.5, 3.5] for T=4
    
    # Set common x-axis limits for both plots
    x_min, x_max = -1.0, T + 0.5
    
    charging_discharge_plot = bar(
        power_x_positions, charging_power_kW,
        dpi=600,
        label="Charging (P_c)",
        color=:green,
        bar_width=1.0,  # Full width bars to exactly touch
        legend=:topleft,
        legendfontsize=8,
        legend_background_color=RGBA(1,1,1,0.7),  # Semi-transparent white background
        xlabel="Time Interval (t)",
        ylabel="P_c/P_d [kW]",
        ylim=ylimit,
        xlim=(x_min, x_max),  # Common x-axis range
        xticks=0:T,  # Show ticks at integer positions
        gridstyle=:solid,
        gridlinewidth=1.0,
        gridalpha=0.2,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        title=title_text,
        titlefont=font(12, "Computer Modern"),
        guidefont=font(15, "Computer Modern"),
        tickfontfamily="Computer Modern",
        top_margin=5Plots.mm
    )
    
    bar!(charging_discharge_plot,
        power_x_positions, -discharging_power_kW,
        dpi=600,
        label="Discharging (P_d)", 
        color=:maroon,  # Wine red color
        bar_width=1.0  # Full width bars to exactly touch
    )
    
    hline!(charging_discharge_plot, [0], color=:black, lw=2, label=false)
    
    # Create SOC bar plot with exact positioning
    # B bars: each covering exactly [t-0.5, t+0.5] for t = 0,1,2,3,4
    # B₀ at x=0: covers [-0.5, 0.5], B₁ at x=1: covers [0.5, 1.5], etc.
    soc_x_positions = [0.0; collect(1.0:T)]  # [0, 1, 2, 3, 4] for T=4
    
    # Plot B0 (initial SOC) with dark plum color to distinguish it as fixed
    soc_plot = bar(
        [soc_x_positions[1]], [soc_percent[1]],  # B0 at x = 0
        dpi=600,
        label="Initial SOC B₀ (Fixed)",
        color=:indigo,  # Much darker color for B₀
        alpha=0.9,
        bar_width=1.0,  # Full width bars to exactly touch
        linewidth=0,
        legend=:top,
        legendfontsize=8,
        legend_background_color=RGBA(1,1,1,0.7),  # Semi-transparent white background
        xlabel="Time Interval (t)",
        ylabel="SOC [%]",
        ylim=(b.soc_min * 100 * 0.95, b.soc_max * 100 * 1.10),
        xlim=(x_min, x_max),  # Same x-axis range as power plot
        xticks=0:T,  # Show ticks at integer positions
        yticks=5*(div(b.soc_min * 100 * 0.95, 5)-1):10:5*(div(b.soc_max * 100 * 1.05, 5)+1),
        gridstyle=:solid,
        gridlinewidth=1.0,
        gridalpha=0.2,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.05,
        title="SOC",
        titlefont=font(12, "Computer Modern"),
        guidefont=font(15, "Computer Modern"),
        tickfontfamily="Computer Modern",
        top_margin=0Plots.mm
    )
    
    # Add B1, B2, B3, B4 (optimized SOC values) with normal styling and exact positioning
    if length(soc_percent) > 1
        bar!(soc_plot,
            soc_x_positions[2:end], soc_percent[2:end],  # B1, B2, B3, B4 at x = 1, 2, 3, 4
            dpi=600,
            label="SOC (B)",
            color=:purple,
            alpha=1.0,
            bar_width=1.0  # Full width bars to exactly touch
        )
    end
    
    # Combine the two plots in vertical layout
    plot_combined = plot(charging_discharge_plot, soc_plot, layout=(2,1))
    
    # Show the plot if requested
    if showPlots
        display(plot_combined)
    end
    
    # Save the plot if requested
    if savePlots
        filename = "$(filename_prefix)_$(lowercase(replace(method_name, " " => "_"))).png"
        @printf "Saving %s battery actions plot to: %s\n" method_name filename
        savefig(plot_combined, filename)
    end
    
    return plot_combined
end

"""
    plot_battery_actions_both(sol_bf, sol_tadmm, inst; showPlots::Bool=true, savePlots::Bool=false, rho_val::Union{Nothing,Float64}=nothing)

Plot battery actions for both Brute Force and tADMM solutions as separate plots.
Uses the same styling as the original Plotter.jl battery actions function.
Includes rho value in tADMM plot title if rho_val is provided.
"""
function plot_battery_actions_both(sol_bf, sol_tadmm, inst::InstancePU; 
    showPlots::Bool=true, savePlots::Bool=false, rho_val::Union{Nothing,Float64}=nothing,
    filename::String="battery_actions.png")
    
    # Extract directory and base filename
    dir_path = dirname(filename)
    base_name = splitext(basename(filename))[1]
    
    # Plot Brute Force solution
    filename_bf = isempty(dir_path) ? "$(base_name)_bf.png" : joinpath(dir_path, "$(base_name)_bf.png")
    plot_bf = plot_battery_actions_single(sol_bf, inst, "Brute Force", 
        showPlots=showPlots, savePlots=savePlots, filename_prefix=filename_bf)
    
    # Plot tADMM solution with rho value in title
    filename_tadmm = isempty(dir_path) ? "$(base_name)_tadmm.png" : joinpath(dir_path, "$(base_name)_tadmm.png")
    plot_tadmm = plot_battery_actions_single(sol_tadmm, inst, "tADMM",
        showPlots=showPlots, savePlots=savePlots, filename_prefix=filename_tadmm, rho_val=rho_val)
    
    return plot_bf, plot_tadmm
end

"""
    plot_power_balance_verification(sol_bf, sol_tadmm, inst; showPlots::Bool=true, savePlots::Bool=false, filename::String="power_balance_verification.png")

Plot power balance verification: P_B + P_Subs = P_L for both Brute Force and tADMM solutions.
Shows all three power components (P_B, P_Subs, P_L) in kW across all time steps for both models.
"""
function plot_power_balance_verification(sol_bf, sol_tadmm, inst::InstancePU; 
    showPlots::Bool=true, savePlots::Bool=false, filename::String="power_balance_verification.png")
    
    T = inst.T
    time_steps = 1:T
    
    # Convert to physical units (kW)
    P_B_bf_kW = sol_bf[:P_B] .* P_BASE
    P_Subs_bf_kW = sol_bf[:P_Subs] .* P_BASE
    P_L_kW = inst.P_L_pu .* P_BASE
    
    P_B_tadmm_kW = sol_tadmm[:P_B] .* P_BASE
    P_Subs_tadmm_kW = sol_tadmm[:P_Subs] .* P_BASE
    
    # Calculate power balance check: P_B + P_Subs (should equal P_L)
    P_total_bf_kW = P_B_bf_kW .+ P_Subs_bf_kW
    P_total_tadmm_kW = P_B_tadmm_kW .+ P_Subs_tadmm_kW
    
    # Set theme and colors
    gr()
    theme(:mute)
    line_colour_bf = :darkorange2
    line_colour_tadmm = :darkorchid1
    line_colour_load = :darkgoldenrod2  # Dark yellow (consistent with input curve plotter)
    line_colour_balance = :crimson      # Complementary red to dark yellow
    
    # Create Brute Force power balance plot
    bf_plot = plot(
        time_steps, P_B_bf_kW,
        dpi=600,
        label="P_B (Battery)",
        xlabel="Time Period (t)",
        ylabel="Power [kW]",
        legend=:bottomleft,
        lw=3,
        color=:purple,
        markershape=:circle,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        xlims=(0.5, T + 0.5),
        xticks=1:T,
        title="Power Balance Verification - Brute Force\nP_B + P_Subs = P_L",
        titlefont=font(12, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
    
    plot!(bf_plot, time_steps, P_Subs_bf_kW,
        label="P_Subs (Substation)",
        lw=3,
        color=line_colour_bf,
        markershape=:square,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5
    )
    
    plot!(bf_plot, time_steps, P_L_kW,
        label="P_L (Load)",
        lw=3,
        color=line_colour_load,
        markershape=:diamond,
        markersize=5,
        markerstrokecolor=:black,
        markerstrokewidth=1.5
    )
    
    plot!(bf_plot, time_steps, P_total_bf_kW,
        label="P_B + P_Subs (Total)",
        lw=2,
        color=line_colour_balance,
        linestyle=:dash,
        markershape=:utriangle,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5
    )
    
    # Create tADMM power balance plot
    tadmm_plot = plot(
        time_steps, P_B_tadmm_kW,
        dpi=600,
        label="P_B (Battery)",
        xlabel="Time Period (t)",
        ylabel="Power [kW]",
        legend=:bottomleft,
        lw=3,
        color=:purple,
        markershape=:circle,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        xlims=(0.5, T + 0.5),
        xticks=1:T,
        title="Power Balance Verification - tADMM (ρ=$(rho))\nP_B + P_Subs = P_L",
        titlefont=font(12, "Computer Modern"),
        guidefont=font(12, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
    
    plot!(tadmm_plot, time_steps, P_Subs_tadmm_kW,
        label="P_Subs (Substation)",
        lw=3,
        color=line_colour_tadmm,
        markershape=:square,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5
    )
    
    plot!(tadmm_plot, time_steps, P_L_kW,
        label="P_L (Load)",
        lw=3,
        color=line_colour_load,
        markershape=:diamond,
        markersize=5,
        markerstrokecolor=:black,
        markerstrokewidth=1.5
    )
    
    plot!(tadmm_plot, time_steps, P_total_tadmm_kW,
        label="P_B + P_Subs (Total)",
        lw=2,
        color=line_colour_balance,
        linestyle=:dash,
        markershape=:utriangle,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5
    )
    
    # Combine plots in vertical layout
    combined_plot = plot(bf_plot, tadmm_plot, layout=(2,1), size=(900, 700))
    
    # Calculate and print power balance errors
    bf_error = maximum(abs.(P_total_bf_kW .- P_L_kW))
    tadmm_error = maximum(abs.(P_total_tadmm_kW .- P_L_kW))
    
    @printf "Power Balance Verification:\n"
    @printf "  Brute Force max error: %.6f kW\n" bf_error
    @printf "  tADMM max error:       %.6f kW\n" tadmm_error
    
    # Show the plot if requested
    if showPlots
        display(combined_plot)
    end
    
    # Save the plot if requested
    if savePlots
        @printf "Saving power balance verification plot to: %s\n" filename
        savefig(combined_plot, filename)
    end
    
    return combined_plot
end

"""
    plot_tadmm_convergence(sol_tadmm, sol_bf; showPlots::Bool=true, savePlots::Bool=false, filename::String="tadmm_convergence.png")

Plot tADMM convergence showing objective function value, primal residual (r_norm), and dual residual (s_norm) across iterations.
Includes horizontal orange line showing Brute Force objective value for comparison.
Uses convergence history from the tADMM solution dictionary.
"""
function plot_tadmm_convergence(sol_tadmm, sol_bf; showPlots::Bool=true, savePlots::Bool=false, filename::String="tadmm_convergence.png")
    
    # Extract convergence history
    conv_hist = sol_tadmm[:convergence_history]
    obj_history = conv_hist[:obj_history]
    r_norm_history = conv_hist[:r_norm_history]
    s_norm_history = conv_hist[:s_norm_history]
    rho_history = conv_hist[:rho_history]
    
    # Downsample data for cleaner plotting (keep all data, just plot subset)
    plot_indices, obj_plot, r_norm_plot, s_norm_plot, rho_plot = downsample_for_plotting(
        obj_history, r_norm_history, s_norm_history, rho_history; max_points=25)
    
    iterations = 1:length(obj_history)  # Full range for x-axis limits
    plot_iterations = plot_indices      # Downsampled points for plotting
    
    # Set theme and colors
    gr()
    theme(:mute)
    line_colour_obj = :dodgerblue       # Blue for tADMM (magenta reserved for DDP)
    line_colour_primal = :darkgreen     # Green for primal residual
    line_colour_dual = :darkorange2     # Orange for dual residual
    
    # Create objective function plot (using downsampled data)
    obj_plot = plot(
        plot_iterations, obj_plot,
        dpi=600,
        label="Objective Value",
        xlabel="Iteration (k)",
        ylabel="Objective Function [\$]",
        legend=:bottomleft,
        lw=3,
        color=line_colour_obj,
        markershape=:circle,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        xlims=(0.5, length(obj_history) + 0.5),
        xticks=1:max(1, div(length(obj_history), 5)):length(obj_history),  # Smart x-ticks
        title="tADMM Convergence - Objective Function",
        titlefont=font(11, "Computer Modern"),
        guidefont=font(11, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
    
    # Add horizontal line for Brute Force objective value
    hline!(obj_plot, [sol_bf[:objective]], color=:darkorange, lw=2, linestyle=:dash, alpha=0.8, label="BF Objective = \$ $(round(sol_bf[:objective], digits=4))")
    
    # Create primal residual plot (log scale)
    primal_plot = plot(
        plot_iterations, r_norm_plot,
        dpi=600,
        label="Primal Residual (‖r‖)",
        xlabel="Iteration (k)",
        ylabel="Primal Residual ‖r‖ [log scale]",
        legend=:bottomleft,
        lw=3,
        color=line_colour_primal,
        markershape=:square,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        xlims=(0.5, length(r_norm_history) + 0.5),
        xticks=1:max(1, div(length(r_norm_history), 5)):length(r_norm_history),
        yscale=:log10,
        title="Primal Residual Convergence",
        titlefont=font(11, "Computer Modern"),
        guidefont=font(11, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
    
    # Add convergence threshold line for primal residual
    hline!(primal_plot, [eps_pri], color=:red, lw=2, linestyle=:dash, alpha=0.7, label="ε_pri = $(eps_pri)")
    
    # Create dual residual plot (log scale)
    dual_plot = plot(
        plot_iterations, s_norm_plot,
        dpi=600,
        label="Dual Residual (‖s‖)",
        xlabel="Iteration (k)",
        ylabel="Dual Residual ‖s‖ [log scale]",
        legend=:bottomleft,
        lw=3,
        color=line_colour_dual,
        markershape=:diamond,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        xlims=(0.5, length(s_norm_history) + 0.5),
        xticks=1:max(1, div(length(s_norm_history), 5)):length(s_norm_history),
        yscale=:log10,
        title="Dual Residual Convergence",
        titlefont=font(11, "Computer Modern"),
        guidefont=font(11, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
    
    # Add convergence threshold line for dual residual
    hline!(dual_plot, [eps_dual], color=:red, lw=2, linestyle=:dash, alpha=0.7, label="ε_dual = $(eps_dual)")
    
    # Create rho evolution plot
    rho_plot_fig = plot(
        plot_iterations, rho_plot,
        dpi=600,
        label="Penalty Parameter ρ",
        xlabel="Iteration (k)",
        ylabel="ρ Value",
        legend=:topright,
        lw=3,
        color=:purple,
        markershape=:circle,
        markersize=4,
        markerstrokecolor=:black,
        markerstrokewidth=1.5,
        gridstyle=:solid,
        gridalpha=0.3,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.15,
        xlims=(0.5, length(rho_history) + 0.5),
        xticks=1:max(1, div(length(rho_history), 5)):length(rho_history),
        title="Penalty Parameter Evolution",
        titlefont=font(11, "Computer Modern"),
        guidefont=font(11, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
    
    # Use log scale if rho varies significantly
    if maximum(rho_history) / minimum(rho_history) > 10
        plot!(rho_plot_fig, yscale=:log10)
    end
    
    # Combine plots in vertical layout (4x1)
    combined_plot = plot(obj_plot, primal_plot, dual_plot, rho_plot_fig, layout=(4,1), size=(800, 1200))
    
    # Print convergence summary
    final_obj = last(obj_history)
    final_r_norm = last(r_norm_history)
    final_s_norm = last(s_norm_history)
    converged = (final_r_norm ≤ eps_pri) && (final_s_norm ≤ eps_dual)
    
    @printf "tADMM Convergence Summary:\n"
    @printf "  Total iterations: %d\n" length(obj_history)
    @printf "  Final objective:  \$%.6f \n" final_obj
    @printf "  Final ‖r‖:        %.2e (threshold: %.1e)\n" final_r_norm eps_pri
    @printf "  Final ‖s‖:        %.2e (threshold: %.1e)\n" final_s_norm eps_dual
    @printf "  Converged:        %s\n" converged ? "✅ YES" : "❌ NO"
    
    # Show the plot if requested
    if showPlots
        display(combined_plot)
    end
    
    # Save the plot if requested
    if savePlots
        @printf "Saving tADMM convergence plot to: %s\n" filename
        savefig(combined_plot, filename)
    end
    
    return combined_plot
end

"""
    plot_objective_breakdown(sol_tadmm; showPlots::Bool=true, savePlots::Bool=false, filename::String="objective_breakdown.png")

Plot breakdown of tADMM objective function showing energy costs, battery costs, and ADMM penalty across iterations.
This helps identify if penalty terms are dominating the primary objective.
"""
# Create output directories for plots
processedData_dir = joinpath(@__DIR__, "envs", "tadmm", "processedData")
system_folder = "$(systemName)_T$(T)"
system_dir = joinpath(processedData_dir, system_folder)
mkpath(system_dir)

println("\n" * "="^80)
println(COLOR_INFO, "GENERATING PLOTS", COLOR_RESET)
println("="^80)
println("Saving plots to: $(system_dir)")
println("="^80)

# Create and display the plots
plot_load_and_cost_curves(showPlots=false, savePlots=true, filename=joinpath(system_dir, "load_and_cost_curves.png"))
plot_battery_actions_both(sol_bf, sol_tadmm, inst, showPlots=false, savePlots=true, rho_val=rho, filename=joinpath(system_dir, "battery_actions.png"))
# plot_power_balance_verification(sol_bf, sol_tadmm, inst, showPlots=false, savePlots=true, filename=joinpath(system_dir, "power_balance_verification.png"))
plot_tadmm_convergence(sol_tadmm, sol_bf, showPlots=false, savePlots=true, filename=joinpath(system_dir, "tadmm_convergence.png"))

println(COLOR_SUCCESS, "✓ All plots saved to $(system_dir)", COLOR_RESET)
println("="^80)

