############################################################
# Copper-Plate MPOPF in per-unit
# - ADMM with exactly T single-step blocks (no K)
# - Brute-force (monolithic) baseline
#
# Bases:
#   kV_B = 4.16/sqrt(3)  (phase-to-neutral; not used in power-only)
#   kVA_B = 1000  => P_BASE = 1000 kW, E_BASE = 1000 kWh (per 1h)
#
# Variables (optimization; all per-unit):
#   P_Subs[t], P_B[t] on P_BASE;  B[t] on E_BASE  (puh)
# Constraints:
#   B[t] = B[t-1] - P_B[t]*Δt
#   P_Subs[t] + P_B[t] = P_L[t]
#   soc_min*E_Rated_pu ≤ B[t] ≤ soc_max*E_Rated_pu
# Objective ($):  sum_t price[t] * (P_Subs_pu[t]*P_BASE) * Δt
#
# Requires: Julia ≥ 1.9, JuMP, Ipopt
############################################################

using JuMP
using Ipopt
using LinearAlgebra
using Printf

# ----------------------- Bases ---------------------------
max_iter = 1000
rho = 5.0                       # ADMM penalty parameter
eps_pri = 1e-3
eps_dual = 1e-3
# ----------------------- Scenario ---------------------------
T = 24
LoadShapeCost = 0.08 .+ 0.12 .* (sin.(range(0, 2π, length=T)) .+ 1) ./ 2
LoadShapeLoad = 0.8 .+ 0.2 .* (sin.(range(0, 2π, length=T) .- 0.8) .+ 1) ./ 2  # Normalized load shape [0, 1]
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
    m = Model(Ipopt.Optimizer)
    set_silent(m)

    # 🔵 Decision variables for time t0
    @variables(m, begin
        # Battery power at time t0 (scalar)
        -b.P_B_R_pu <= P_B_t0 <= b.P_B_R_pu         
        # Substation power at time t0 (scalar)  
        P_subs_t0                                     
        # 🔵 Local SOC trajectory (T-length vector, but only B_t0[t0] is truly optimized)
        Bmin(b) <= B_t0_var[1:inst.T] <= Bmax(b)                  
    end)

    # 📌 CONSTRAINT 1: SOC Trajectory (1 equality constraint)
    # 🔵B_t0[t0] = 🔴B̂[t0-1] - P_B_t0 * Δt
    if t0 == 1
        # Special case: use B0 for t0-1
        @constraint(m, B_t0_var[t0] == b.B0_pu - P_B_t0 * b.Δt)
    else
        @constraint(m, B_t0_var[t0] == Bhat[t0-1] - P_B_t0 * b.Δt)
    end
    
    # 📌 CONSTRAINT 2: Nodal Real Power Balance (1 equality constraint)
    # P_subs_t0 + P_B_t0 = P_L[t0]
    @constraint(m, P_subs_t0 + P_B_t0 == inst.P_L_pu[t0])

    # 🎯 OBJECTIVE: Economic cost + ADMM penalty
    # C_t0 * P_subs_t0 + (ρ/2) * ||🔵B_t0 - 🔴B̂ + 🟢u_t0||²₂
    energy_cost = inst.price[t0] * (P_subs_t0 * P_BASE) * b.Δt
    
    # ADMM penalty: (ρ/2) * ||🔵B_t0 - 🔴B̂ + 🟢u_t0||²₂
    penalty = (ρ / 2) * sum((B_t0_var[t] - Bhat[t] + u_t0[t])^2 for t in 1:inst.T)
    
    @objective(m, Min, energy_cost + penalty)

    # 🚀 Solve subproblem
    optimize!(m)

    # 📥 Extract results and update 🔵B_t0
    # B_t0[t0] = value(B_t0_var[t0])  # Only update the t0-th component
    for t in 1:inst.T
        B_t0[t] = value(B_t0_var[t])  # Update ALL components
    end
    P_B_val = value(P_B_t0)
    P_subs_val = value(P_subs_t0)
    
    # 📦 Return results as extensible dictionary
    return Dict(
        :objective => objective_value(m),
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

    # � Update consensus variables: B̂[t] for t = 1, 2, ..., T-1  
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
    max_iter::Int=200, eps_pri::Float64=1e-3, eps_dual::Float64=1e-3)
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

    # 📊 History tracking
    obj_history = Float64[]
    Bhat_history = Vector{Vector{Float64}}()
    B_collection_history = Vector{Vector{Vector{Float64}}}()
    u_collection_history = Vector{Vector{Vector{Float64}}}()
    r_norm_history = Float64[]
    s_norm_history = Float64[]

    # Store initial states
    push!(Bhat_history, copy(Bhat))
    push!(B_collection_history, deepcopy(B_collection))
    push!(u_collection_history, deepcopy(u_collection))

    @printf "🎯 tADMM[PDF-formulation]: T=%d, ρ=%.3f\n" T ρ

    for k in 1:max_iter
        # 🔵 STEP 1: Primal Update - Solve T subproblems
        @printf "  🔵 Primal updates: "
        total_obj = 0.0
        for t0 in 1:T
            result = primal_update_tadmm!(B_collection[t0], Bhat, u_collection[t0], inst, ρ, t0)
            total_obj += result[:objective]
            # @printf "%d " t0
        end
        @printf "\n"
        push!(obj_history, total_obj)

        # 🔴 STEP 2: Consensus Update  
        Bhat_old = copy(Bhat)
        consensus_result = consensus_update_tadmm!(Bhat, B_collection, u_collection, inst, ρ)

        # 🟢 STEP 3: Dual Update
        dual_result = dual_update_tadmm!(u_collection, B_collection, Bhat, ρ)

        # 📊 Store iteration history
        push!(Bhat_history, copy(Bhat))
        push!(B_collection_history, deepcopy(B_collection))
        push!(u_collection_history, deepcopy(u_collection))

        # 📏 STEP 4: Compute residuals using vector norms
        # Primal residual: measure consensus violation using 2-norm
        r_vectors = []
        for t0 in 1:T
            push!(r_vectors, B_collection[t0] - Bhat)
        end
        r_norm = norm(vcat(r_vectors...))  # Concatenate all residual vectors and take 2-norm

        # Dual residual: measure change in consensus
        s_norm = ρ * norm(Bhat - Bhat_old)

        push!(r_norm_history, r_norm)
        push!(s_norm_history, s_norm)

        @printf "k=%3d  obj=%10.4f  ‖r‖=%.2e  ‖s‖=%.2e\n" k total_obj r_norm s_norm

        if r_norm ≤ eps_pri && s_norm ≤ eps_dual
            @printf "🎉 tADMM converged at iteration %d\n" k
            break
        end
    end

    # 📤 Compute final P_B and P_Subs from consensus solution
    # Reconstruct power schedules from final consensus trajectory
    P_B_final = zeros(T)
    P_Subs_final = zeros(T)

    for t in 1:T
        if t == 1
            P_B_final[t] = (b.B0_pu - Bhat[t]) / b.Δt
        else
            P_B_final[t] = (Bhat[t-1] - Bhat[t]) / b.Δt
        end
        P_Subs_final[t] = inst.P_L_pu[t] - P_B_final[t]
    end

    return Dict(
        :P_B => P_B_final,
        :P_Subs => P_Subs_final,
        :B => Bhat,  # Use consensus SOC trajectory
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
            :obj_history => obj_history
        )
    )
end

# ---------------- Brute-force (monolithic, pu) -----------
function solve_MPOPF_using_BruteForce(inst::InstancePU)
    T = inst.T
    b = inst.bat
    m = Model(Ipopt.Optimizer)
    set_silent(m)

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

    @objective(m, Min, sum(inst.price[t] * (P_Subs[t] * P_BASE) * b.Δt for t in 1:T))
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

# Brute-force baseline (unchanged)
sol_bf = solve_MPOPF_using_BruteForce(inst)
@printf "\nBrute-force objective: %.4f\n" sol_bf[:objective]

# tADMM with PDF formulation variable names 
sol_tadmm = solve_MPOPF_using_tADMM(inst; ρ=rho, max_iter=max_iter, eps_pri=eps_pri, eps_dual=eps_dual)
@printf "tADMM objective: %.4f\n" sol_tadmm[:objective]

# Quick pu checks - compare all methods
@printf "BruteForced     B[min,max]=[%.4f, %.4f] pu  | PB[min,max]=[%.4f, %.4f] pu\n" minimum(sol_bf[:B]) maximum(sol_bf[:B]) minimum(sol_bf[:P_B]) maximum(sol_bf[:P_B])

@printf "tADMM  B[min,max]=[%.4f, %.4f] pu  | PB[min,max]=[%.4f, %.4f] pu\n" minimum(sol_tadmm[:B]) maximum(sol_tadmm[:B]) minimum(sol_tadmm[:P_B]) maximum(sol_tadmm[:P_B])

# Convert tADMM results back to physical for sanity checks
P_B_kW = sol_tadmm[:P_B] .* P_BASE
PSubs_kW = sol_tadmm[:P_Subs] .* P_BASE
B_kWh = sol_tadmm[:B] .* E_BASE
@printf "(tADMM) P_B[kW] in [%.1f, %.1f], B[kWh] in [%.1f, %.1f]\n" minimum(P_B_kW) maximum(P_B_kW) minimum(B_kWh) maximum(B_kWh)

# ----------------------- Plotting Functions --------------
using Plots

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
        legend=:topleft,
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
        legend=:topright,
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

# Create and display the plot
plot_load_and_cost_curves(showPlots=true, savePlots=true, filename="load_and_cost_curves.png")

