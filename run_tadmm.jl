# ============================================================================
# run_tadmm.jl — Temporal ADMM decomposed SOCP solve for MPOPF
# ============================================================================
# Usage: edit config.jl, then run this script
#   julia run_tadmm.jl
#   OR in VS Code: Ctrl+Shift+Enter
# ============================================================================

# Fix GR backend display issues on headless/Windows environments
ENV["GKSwstype"] = "png"

# Load shared configuration
include("config.jl")

# Activate the tadmm environment
import Pkg
Pkg.activate(ENV_PATH)

using Revise
using LinearAlgebra
using JuMP
using Ipopt
if USE_GUROBI
    using Gurobi
end
using OpenDSSDirect
using Printf
using Statistics
using Parameters: @unpack
using Crayons
using Serialization

const MOI = JuMP.MOI

# Threading status
println("\n" * "="^80)
println("PARALLELIZATION STATUS")
println("="^80)
println("Julia threads:       ", Threads.nthreads())
if Threads.nthreads() > 1
    println("Status:              Using $(Threads.nthreads()) threads")
else
    println("Status:              Single-threaded mode")
    println("  To enable parallel: JULIA_NUM_THREADS=16 julia run_tadmm.jl")
end
println("="^80)

const USE_THREADING = Threads.nthreads() > 1

# Create shared Gurobi environment
const GUROBI_ENV = USE_GUROBI ? Gurobi.Env() : nothing

# Color scheme
const COLOR_SUCCESS = Crayon(foreground = :green, bold = true)
const COLOR_WARNING = Crayon(foreground = :yellow, bold = true)
const COLOR_ERROR = Crayon(foreground = :red, bold = true)
const COLOR_INFO = Crayon(foreground = :cyan, bold = true)
const COLOR_HIGHLIGHT = Crayon(foreground = :magenta, bold = true)
const COLOR_RESET = Crayon(reset = true)

# Include utilities
includet(joinpath(ENV_PATH, "parse_opendss.jl"))
includet(joinpath(ENV_PATH, "opendss_validator.jl"))
includet(joinpath(ENV_PATH, "solution_validator.jl"))
includet(joinpath(ENV_PATH, "logger.jl"))

# ============================================================================
# tADMM PARAMETERS
# ============================================================================

rho_base = 3000.0
rho_scaling_with_T = true
rho_tadmm = rho_scaling_with_T ? rho_base * sqrt(T / 24.0) : rho_base
max_iter_tadmm = 500
eps_pri_tadmm = 1e-3
eps_dual_tadmm = 1e-2
adaptive_rho_tadmm = true
stagnation_window = 5
stagnation_threshold = 1e-3

enable_warm_start = true
track_subproblem_details = true

# ============================================================================
# PARSE SYSTEM DATA
# ============================================================================

println("\n" * "="^80)
println("PARSING SYSTEM DATA")
println("="^80)

data = parse_system_from_dss(
    SYSTEM_NAME, T;
    LoadShapeLoad=LoadShapeLoad,
    LoadShapeCost=LoadShapeCost,
    LoadShapePV=LoadShapePV,
    C_B=C_B,
    delta_t_h=DELTA_T_H
)

print_powerflow_summary(data)

# ============================================================================
# PRIMAL UPDATE (Step 1: solve per-time-period subproblems)
# ============================================================================

function primal_update_tadmm_socp!(B_local, Bhat, u_local, data, rho::Float64, t0::Int;
                                   solver::Symbol=:ipopt,
                                   warm_start::Bool=false,
                                   prev_solution::Union{Nothing,Dict}=nothing,
                                   track_subproblem_details::Bool=false,
                                   lindistflow::Bool=false)
    @unpack Nset, Lset, L1set, Nm1set, NLset, Dset, Bset, Tset = data
    @unpack substationBus, parent, children = data
    @unpack rdict_pu, xdict_pu, Vminpu, Vmaxpu = data
    @unpack p_L_pu, q_L_pu, p_D_pu, S_D_R = data
    @unpack B0_pu, B_R_pu, P_B_R_pu, soc_min, soc_max = data
    @unpack LoadShapeCost, C_B, delta_t_h, kVA_B = data

    j1 = substationBus
    dt = delta_t_h
    P_BASE = kVA_B

    # Create model for subproblem t0
    model = Model()
    if solver == :gurobi
        set_optimizer(model, () -> Gurobi.Optimizer(GUROBI_ENV))
        set_silent(model)
        set_optimizer_attribute(model, "NonConvex", 2)
        set_optimizer_attribute(model, "OutputFlag", 0)
        set_optimizer_attribute(model, "Threads", 1)
        set_optimizer_attribute(model, "DualReductions", 0)
    else
        set_optimizer(model, Ipopt.Optimizer)
        set_silent(model)
        set_optimizer_attribute(model, "print_level", 0)
        set_optimizer_attribute(model, "max_iter", 30)
        set_optimizer_attribute(model, "tol", 1e-6)
        set_optimizer_attribute(model, "acceptable_tol", 1e-4)
        set_optimizer_attribute(model, "mu_strategy", "adaptive")
        set_optimizer_attribute(model, "linear_solver", "mumps")
        set_optimizer_attribute(model, "nlp_scaling_method", "gradient-based")
        set_optimizer_attribute(model, "bound_relax_factor", 1e-8)
    end

    # ===== NETWORK VARIABLES (time t0 ONLY) =====
    @variable(model, P_Subs_t0 >= 0)
    @variable(model, Q_Subs_t0)
    @variable(model, P_t0[(i, j) in Lset])
    @variable(model, Q_t0[(i, j) in Lset])
    @variable(model, v_t0[j in Nset])
    if !lindistflow
        @variable(model, l_t0[(i, j) in Lset] >= 1e-6)
    end
    @variable(model, q_D_t0[j in Dset])

    # ===== BATTERY VARIABLES (LOCAL COUPLING - NEIGHBORS ONLY) =====
    if t0 == 1
        local_times = [1, 2]
    elseif t0 == length(Tset)
        local_times = [t0-1, t0]
    else
        local_times = [t0-1, t0, t0+1]
    end

    @variable(model, P_B_var[j in Bset, t in Tset])
    @variable(model, B_var[j in Bset, t in local_times])

    # ===== OBJECTIVE =====
    energy_cost = LoadShapeCost[t0] * P_Subs_t0 * P_BASE * dt
    battery_cost = sum(C_B * (P_B_var[j, t0] * P_BASE)^2 * dt for j in Bset)
    penalty = (rho / 2) * sum(sum((B_var[j, t] - Bhat[j][t] + u_local[j][t])^2 for t in local_times) for j in Bset)
    @objective(model, Min, energy_cost + battery_cost + penalty)

    # ===== SPATIAL CONSTRAINTS (time t0 ONLY) =====

    # Real power balance
    @constraint(model, P_Subs_t0 - sum(P_t0[(j1, j)] for (j1, j) in L1set) == 0)
    for j in Nm1set
        i = parent[j]
        sum_Pjk = isempty(children[j]) ? 0.0 : sum(P_t0[(j, k)] for k in children[j])
        p_L_j = (j in NLset) ? p_L_pu[j, t0] : 0.0
        p_D_j = (j in Dset) ? p_D_pu[j, t0] : 0.0
        P_B_j = (j in Bset) ? P_B_var[j, t0] : 0.0
        @constraint(model, sum_Pjk - P_t0[(i, j)] == P_B_j + p_D_j - p_L_j)
    end

    # Reactive power balance
    @constraint(model, Q_Subs_t0 - sum(Q_t0[(j1, j)] for (j1, j) in L1set) == 0)
    for j in Nm1set
        i = parent[j]
        sum_Qjk = isempty(children[j]) ? 0.0 : sum(Q_t0[(j, k)] for k in children[j])
        q_L_j = (j in NLset) ? q_L_pu[j, t0] : 0.0
        q_D_j = (j in Dset) ? q_D_t0[j] : 0.0
        @constraint(model, sum_Qjk - Q_t0[(i, j)] == q_D_j - q_L_j)
    end

    # Voltage drop
    if lindistflow
        for (i, j) in Lset
            r_ij = rdict_pu[(i, j)]
            x_ij = xdict_pu[(i, j)]
            @constraint(model, v_t0[j] == v_t0[i] - 2*(r_ij*P_t0[(i,j)] + x_ij*Q_t0[(i,j)]))
        end
    else
        for (i, j) in Lset
            r_ij = rdict_pu[(i, j)]
            x_ij = xdict_pu[(i, j)]
            @constraint(model, v_t0[j] == v_t0[i] - 2*(r_ij*P_t0[(i,j)] + x_ij*Q_t0[(i,j)]) + (r_ij^2 + x_ij^2)*l_t0[(i,j)])
        end
        for (i, j) in Lset
            @constraint(model, P_t0[(i,j)]^2 + Q_t0[(i,j)]^2 <= v_t0[i] * l_t0[(i,j)])
        end
    end

    # Fixed substation voltage
    @constraint(model, v_t0[j1] == 1.05^2)

    # Voltage limits
    for j in Nset
        @constraint(model, Vminpu[j]^2 <= v_t0[j] <= Vmaxpu[j]^2)
    end

    # PV reactive power limits
    for j in Dset
        p_D_val = p_D_pu[j, t0]
        S_D_R_val = S_D_R[j]
        q_max = sqrt(S_D_R_val^2 - p_D_val^2)
        @constraint(model, -q_max <= q_D_t0[j] <= q_max)
    end

    # ===== TEMPORAL CONSTRAINTS (LOCAL COUPLING) =====
    for j in Bset
        if t0 == 1
            @constraint(model, B_var[j, 1] == B0_pu[j] - P_B_var[j, 1] * dt)
            @constraint(model, B_var[j, 2] == B_var[j, 1] - P_B_var[j, 2] * dt)
        elseif t0 == length(Tset)
            @constraint(model, B_var[j, t0] == B_var[j, t0-1] - P_B_var[j, t0] * dt)
        else
            @constraint(model, B_var[j, t0] == B_var[j, t0-1] - P_B_var[j, t0] * dt)
            @constraint(model, B_var[j, t0+1] == B_var[j, t0] - P_B_var[j, t0+1] * dt)
        end
        for t in local_times
            @constraint(model, soc_min[j] * B_R_pu[j] <= B_var[j, t] <= soc_max[j] * B_R_pu[j])
        end
        for t in local_times
            @constraint(model, -P_B_R_pu[j] <= P_B_var[j, t] <= P_B_R_pu[j])
        end
    end

    # ===== WARM-START =====
    if warm_start && !isnothing(prev_solution)
        if haskey(prev_solution, :P_Subs)
            set_start_value(P_Subs_t0, prev_solution[:P_Subs])
        end
        if haskey(prev_solution, :Q_Subs)
            set_start_value(Q_Subs_t0, prev_solution[:Q_Subs])
        end
        if haskey(prev_solution, :P)
            for (i, j) in Lset
                if haskey(prev_solution[:P], (i, j))
                    set_start_value(P_t0[(i, j)], prev_solution[:P][(i, j)])
                end
            end
        end
        if haskey(prev_solution, :Q)
            for (i, j) in Lset
                if haskey(prev_solution[:Q], (i, j))
                    set_start_value(Q_t0[(i, j)], prev_solution[:Q][(i, j)])
                end
            end
        end
        if haskey(prev_solution, :v)
            for j in Nset
                if haskey(prev_solution[:v], j)
                    set_start_value(v_t0[j], prev_solution[:v][j])
                end
            end
        end
        if !lindistflow && haskey(prev_solution, :l)
            for (i, j) in Lset
                if haskey(prev_solution[:l], (i, j))
                    set_start_value(l_t0[(i, j)], prev_solution[:l][(i, j)])
                end
            end
        end
        if haskey(prev_solution, :q_D)
            for j in Dset
                if haskey(prev_solution[:q_D], j)
                    set_start_value(q_D_t0[j], prev_solution[:q_D][j])
                end
            end
        end
        if haskey(prev_solution, :P_B)
            for j in Bset
                if haskey(prev_solution[:P_B], j)
                    for t in Tset
                        if haskey(prev_solution[:P_B][j], t)
                            set_start_value(P_B_var[j, t], prev_solution[:P_B][j][t])
                        end
                    end
                end
            end
        end
        if haskey(prev_solution, :B_local)
            for j in Bset
                if haskey(prev_solution[:B_local], j)
                    for t in local_times
                        if haskey(prev_solution[:B_local][j], t)
                            set_start_value(B_var[j, t], prev_solution[:B_local][j][t])
                        end
                    end
                end
            end
        end
    end

    # Model stats (only for first subproblem)
    subproblem_stats = (t0 == 1) ? get_model_size_statistics(model) : nothing

    # Solve
    optimize!(model)

    status = termination_status(model)
    acceptable = (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.ALMOST_LOCALLY_SOLVED || status == MOI.ITERATION_LIMIT)
    if !acceptable
        error("Subproblem t0=$t0 failed with status: $status")
    end

    solver_time_t0 = solve_time(model)

    ipopt_iters = 0
    if solver == :ipopt && track_subproblem_details
        try
            ipopt_iters = MOI.get(model, MOI.BarrierIterations())
        catch
            ipopt_iters = 0
        end
    end

    # Extract results
    for j in Bset
        for t in local_times
            B_local[j][t] = value(B_var[j, t])
        end
    end

    P_Subs_val = value(P_Subs_t0)
    Q_Subs_val = value(Q_Subs_t0)
    P_vals = Dict((i, j) => value(P_t0[(i, j)]) for (i, j) in Lset)
    Q_vals = Dict((i, j) => value(Q_t0[(i, j)]) for (i, j) in Lset)
    v_vals = Dict(j => value(v_t0[j]) for j in Nset)
    if lindistflow
        l_vals = Dict((i, j) => (P_vals[(i,j)]^2 + Q_vals[(i,j)]^2) / v_vals[i] for (i, j) in Lset)
    else
        l_vals = Dict((i, j) => value(l_t0[(i, j)]) for (i, j) in Lset)
    end
    q_D_vals = Dict(j => value(q_D_t0[j]) for j in Dset)
    P_B_vals_all_times = Dict(j => Dict(t => value(P_B_var[j, t]) for t in Tset) for j in Bset)
    P_B_val_t0 = Dict(j => value(P_B_var[j, t0]) for j in Bset)

    energy_cost_val = LoadShapeCost[t0] * P_Subs_val * P_BASE * dt
    battery_cost_val = sum(C_B * (P_B_val_t0[j] * P_BASE)^2 * dt for j in Bset)
    penalty_val = (rho / 2) * sum(sum((value(B_var[j, t]) - Bhat[j][t] + u_local[j][t])^2 for t in local_times) for j in Bset)

    return Dict(
        :total_objective => objective_value(model),
        :energy_cost => energy_cost_val,
        :battery_cost => battery_cost_val,
        :penalty => penalty_val,
        :solve_time => solver_time_t0,
        :ipopt_iters => ipopt_iters,
        :P_Subs => P_Subs_val,
        :Q_Subs => Q_Subs_val,
        :P => P_vals,
        :Q => Q_vals,
        :v => v_vals,
        :l => l_vals,
        :q_D => q_D_vals,
        :P_B => P_B_vals_all_times,
        :B_local => B_local,
        :t0 => t0,
        :model_stats => subproblem_stats
    )
end

# ============================================================================
# CONSENSUS UPDATE (Step 2: average local SOC across overlapping subproblems)
# ============================================================================

function consensus_update_tadmm_socp!(Bhat, B_collection, u_collection, data, rho::Float64)
    @unpack Bset, Tset, B_R_pu, soc_min, soc_max = data

    for j in Bset
        for t in Tset
            if t == 1
                subproblems_with_t = [1, 2]
            elseif t == length(Tset)
                subproblems_with_t = [length(Tset)-1, length(Tset)]
            else
                subproblems_with_t = [t-1, t, t+1]
            end

            consensus_sum = 0.0
            num_contributors = 0
            for t0 in subproblems_with_t
                if haskey(B_collection[t0][j], t) && haskey(u_collection[t0][j], t)
                    consensus_sum += B_collection[t0][j][t] + u_collection[t0][j][t]
                    num_contributors += 1
                end
            end
            Bhat_new = num_contributors > 0 ? consensus_sum / num_contributors : 0.0

            B_min = soc_min[j] * B_R_pu[j]
            B_max = soc_max[j] * B_R_pu[j]
            Bhat[j][t] = clamp(Bhat_new, B_min, B_max)
        end
    end
end

# ============================================================================
# DUAL UPDATE (Step 3: update scaled dual variables)
# ============================================================================

function dual_update_tadmm_socp!(u_collection, B_collection, Bhat, rho::Float64, data)
    @unpack Bset, Tset = data

    for t0 in Tset
        for j in Bset
            for t in keys(u_collection[t0][j])
                u_collection[t0][j][t] += (B_collection[t0][j][t] - Bhat[j][t])
            end
        end
    end
end

# ============================================================================
# MAIN tADMM LOOP
# ============================================================================

function solve_MPOPF_SOCP_tADMM(data; rho::Float64=1.0,
                                    max_iter::Int=1000, eps_pri::Float64=1e-5, eps_dual::Float64=1e-4,
                                    adaptive_rho::Bool=false, solver::Symbol=:ipopt,
                                    enable_warm_start::Bool=false, track_subproblem_details::Bool=false)
    @unpack Bset, Tset, B0_pu, B_R_pu, soc_min, soc_max = data

    # Initialize consensus variables
    Bhat = Dict(j => fill(B0_pu[j], length(Tset)) for j in Bset)
    for j in Bset
        Bhat[j] .= clamp.(Bhat[j], soc_min[j] * B_R_pu[j], soc_max[j] * B_R_pu[j])
    end

    # Initialize local SOC and dual variables per subproblem
    B_collection = Dict{Int,Dict{Int,Dict{Int,Float64}}}()
    u_collection = Dict{Int,Dict{Int,Dict{Int,Float64}}}()

    for t0 in Tset
        if t0 == 1
            local_times = [1, 2]
        elseif t0 == length(Tset)
            local_times = [t0-1, t0]
        else
            local_times = [t0-1, t0, t0+1]
        end
        B_collection[t0] = Dict(j => Dict(t => Bhat[j][t] for t in local_times) for j in Bset)
        u_collection[t0] = Dict(j => Dict(t => 0.0 for t in local_times) for j in Bset)
    end

    # Storage for decision variables
    P_Subs_collection = Dict{Int,Float64}()
    Q_Subs_collection = Dict{Int,Float64}()
    P_collection = Dict{Tuple{Int,Int},Dict{Int,Float64}}()
    Q_collection = Dict{Tuple{Int,Int},Dict{Int,Float64}}()
    v_collection = Dict{Int,Dict{Int,Float64}}()
    l_collection = Dict{Tuple{Int,Int},Dict{Int,Float64}}()
    q_D_collection = Dict{Int,Dict{Int,Float64}}()
    P_B_collection = Dict{Int,Dict{Int,Dict{Int,Float64}}}()

    for (i, j) in data[:Lset]
        P_collection[(i, j)] = Dict{Int,Float64}()
        Q_collection[(i, j)] = Dict{Int,Float64}()
        l_collection[(i, j)] = Dict{Int,Float64}()
    end
    for j in data[:Nset]
        v_collection[j] = Dict{Int,Float64}()
    end
    for j in data[:Dset]
        q_D_collection[j] = Dict{Int,Float64}()
    end
    for j in Bset
        P_B_collection[j] = Dict{Int,Dict{Int,Float64}}()
    end

    # History tracking
    obj_history = Float64[]
    r_norm_history = Float64[]
    s_norm_history = Float64[]
    rho_history = Float64[]
    subproblem_times_history = Vector{Vector{Float64}}()
    iteration_effective_times = Float64[]

    # Detailed per-subproblem tracking
    n_T = length(Tset)
    subproblem_times_matrix = zeros(Float64, max_iter, n_T)
    subproblem_iters_matrix = zeros(Int, max_iter, n_T)

    # Adaptive rho parameters
    rho_current = rho
    mu_balance = 5.0
    tau_incr = 2.0
    tau_decr = 2.0
    rho_min = 1.0
    rho_max = 1e6
    update_interval = 5

    # Two-phase adaptive rho
    primal_converged_once = false
    phase1_only_increase = true

    # Primal residual watchdog
    enable_rnorm_watchdog = true
    rnorm_watchdog_window = 15
    rnorm_watchdog_factor = 2.5
    rnorm_high_counter = 0

    # Output control
    progress_interval = 5

    println("\n" * "="^80)
    print(COLOR_INFO)
    @printf "tADMM[SOCP]: T=%d, rho_init=%.1f, adaptive=%s, warm_start=%s, |Bset|=%d, |Nset|=%d, |Lset|=%d\n" length(Tset) rho adaptive_rho enable_warm_start length(Bset) length(data[:Nset]) length(data[:Lset])
    print(COLOR_RESET)
    println("="^80)

    converged = false
    final_iter = max_iter

    tadmm_subproblem_stats = nothing
    prev_solutions = Dict{Int, Union{Nothing,Dict}}()
    for t0 in Tset
        prev_solutions[t0] = nothing
    end

    tadmm_wallclock_start = time()

    try
        tadmm_start_time = time()
        for k in 1:max_iter
            # STEP 1: Primal Update
            total_energy_cost = 0.0
            total_battery_cost = 0.0
            total_penalty = 0.0
            subproblem_times_k = Float64[]

            if USE_THREADING
                results_collection = Vector{Any}(undef, length(Tset))
                Threads.@threads for i in 1:length(Tset)
                    t0 = Tset[i]
                    results_collection[i] = primal_update_tadmm_socp!(
                        B_collection[t0], Bhat, u_collection[t0], data, rho_current, t0;
                        solver=solver,
                        warm_start=enable_warm_start && k > 1,
                        prev_solution=prev_solutions[t0],
                        track_subproblem_details=track_subproblem_details
                    )
                end

                for i in 1:length(Tset)
                    t0 = Tset[i]
                    result = results_collection[i]
                    push!(subproblem_times_k, result[:solve_time])

                    if track_subproblem_details
                        subproblem_idx = findfirst(==(t0), Tset)
                        subproblem_times_matrix[k, subproblem_idx] = result[:solve_time]
                        if haskey(result, :ipopt_iters)
                            subproblem_iters_matrix[k, subproblem_idx] = result[:ipopt_iters]
                        end
                    end

                    if k == 1 && t0 == 1 && !isnothing(result[:model_stats])
                        tadmm_subproblem_stats = result[:model_stats]
                    end
                    if enable_warm_start
                        prev_solutions[t0] = result
                    end

                    B_collection[t0] = result[:B_local]
                    P_Subs_collection[t0] = result[:P_Subs]
                    Q_Subs_collection[t0] = result[:Q_Subs]
                    for (i, j) in data[:Lset]
                        P_collection[(i, j)][t0] = result[:P][(i, j)]
                        Q_collection[(i, j)][t0] = result[:Q][(i, j)]
                        l_collection[(i, j)][t0] = result[:l][(i, j)]
                    end
                    for j in data[:Nset]
                        v_collection[j][t0] = result[:v][j]
                    end
                    for j in data[:Dset]
                        q_D_collection[j][t0] = result[:q_D][j]
                    end
                    for j in Bset
                        if !haskey(P_B_collection[j], t0)
                            P_B_collection[j][t0] = Dict{Int,Float64}()
                        end
                        for t in Tset
                            P_B_collection[j][t0][t] = result[:P_B][j][t]
                        end
                    end

                    total_energy_cost += result[:energy_cost]
                    total_battery_cost += result[:battery_cost]
                    total_penalty += result[:penalty]
                end
            else
                for t0 in Tset
                    result = primal_update_tadmm_socp!(
                        B_collection[t0], Bhat, u_collection[t0], data, rho_current, t0;
                        solver=solver,
                        warm_start=enable_warm_start && k > 1,
                        prev_solution=prev_solutions[t0],
                        track_subproblem_details=track_subproblem_details
                    )
                    push!(subproblem_times_k, result[:solve_time])

                    if track_subproblem_details
                        subproblem_idx = findfirst(==(t0), Tset)
                        subproblem_times_matrix[k, subproblem_idx] = result[:solve_time]
                        if haskey(result, :ipopt_iters)
                            subproblem_iters_matrix[k, subproblem_idx] = result[:ipopt_iters]
                        end
                    end

                    if k == 1 && t0 == 1 && !isnothing(result[:model_stats])
                        tadmm_subproblem_stats = result[:model_stats]
                    end
                    if enable_warm_start
                        prev_solutions[t0] = result
                    end

                    B_collection[t0] = result[:B_local]
                    P_Subs_collection[t0] = result[:P_Subs]
                    Q_Subs_collection[t0] = result[:Q_Subs]
                    for (i, j) in data[:Lset]
                        P_collection[(i, j)][t0] = result[:P][(i, j)]
                        Q_collection[(i, j)][t0] = result[:Q][(i, j)]
                        l_collection[(i, j)][t0] = result[:l][(i, j)]
                    end
                    for j in data[:Nset]
                        v_collection[j][t0] = result[:v][j]
                    end
                    for j in data[:Dset]
                        q_D_collection[j][t0] = result[:q_D][j]
                    end
                    for j in Bset
                        P_B_collection[j][t0] = result[:P_B][j]
                    end

                    total_energy_cost += result[:energy_cost]
                    total_battery_cost += result[:battery_cost]
                    total_penalty += result[:penalty]
                end
            end

            true_objective = total_energy_cost + total_battery_cost
            push!(obj_history, true_objective)
            push!(subproblem_times_history, subproblem_times_k)
            push!(iteration_effective_times, maximum(subproblem_times_k))

            # STEP 2: Consensus Update
            Bhat_old = Dict(j => copy(Bhat[j]) for j in Bset)
            consensus_update_tadmm_socp!(Bhat, B_collection, u_collection, data, rho_current)

            # STEP 3: Dual Update
            dual_update_tadmm_socp!(u_collection, B_collection, Bhat, rho_current, data)

            # STEP 4: Compute residuals
            r_values = Float64[]
            for t0 in Tset
                for j in Bset
                    if haskey(B_collection[t0][j], t0)
                        push!(r_values, B_collection[t0][j][t0] - Bhat[j][t0])
                    end
                end
            end
            r_norm = norm(r_values) / sqrt(length(r_values))

            s_vectors = []
            for j in Bset
                push!(s_vectors, Bhat[j] - Bhat_old[j])
            end
            all_s_values = vcat(s_vectors...)
            s_norm = rho_current * norm(all_s_values) / sqrt(length(all_s_values))

            push!(r_norm_history, r_norm)
            push!(s_norm_history, s_norm)
            push!(rho_history, rho_current)

            # Watchdog: monitor primal residual trend
            if enable_rnorm_watchdog
                if r_norm > 2.0 * eps_pri
                    rnorm_high_counter += 1
                elseif r_norm < 0.5 * eps_pri
                    rnorm_high_counter = 0
                end
                if rnorm_high_counter >= rnorm_watchdog_window
                    rho_old = rho_current
                    rho_current = min(rho_max, rnorm_watchdog_factor * rho_current)
                    @printf "  [k=%3d] RNORM WATCHDOG: rho %.1f -> %.1f (r=%.2e > %.1e for %d iters)\n" k rho_old rho_current r_norm eps_pri rnorm_high_counter
                    scale_factor = rho_old / rho_current
                    for t0 in Tset, j in Bset
                        for t_key in keys(u_collection[t0][j])
                            u_collection[t0][j][t_key] *= scale_factor
                        end
                    end
                    rnorm_high_counter = 0
                end
            end

            # STEP 5: Adaptive rho (two-phase)
            if adaptive_rho && k % update_interval == 0 && k < max_iter - 50
                rho_old = rho_current
                rho_changed = false

                if r_norm <= eps_pri && !primal_converged_once
                    primal_converged_once = true
                    phase1_only_increase = false
                    @printf "  [k=%3d] Phase 1->2: primal converged\n" k
                end

                if phase1_only_increase
                    if r_norm > mu_balance * s_norm
                        rho_current = min(rho_max, tau_incr * rho_current)
                        @printf "  [PHASE1-UP] rho: %.1f -> %.1f (r/s=%.1f)\n" rho_old rho_current (r_norm/s_norm)
                        rho_changed = true
                    end
                else
                    if r_norm > mu_balance * s_norm
                        rho_current = min(rho_max, tau_incr * rho_current)
                        @printf "  [k=%3d] rho %.1f -> %.1f (r/s=%.1f)\n" k rho_old rho_current (r_norm/s_norm)
                        rho_changed = true
                    elseif s_norm > mu_balance * r_norm && r_norm <= eps_pri
                        rho_current = max(rho_min, rho_current / tau_decr)
                        @printf "  [k=%3d] rho %.1f -> %.1f (s/r=%.1f, r converged)\n" k rho_old rho_current (s_norm/r_norm)
                        rho_changed = true
                    end
                end

                if rho_current != rho_old
                    scale_factor = rho_old / rho_current
                    for t0 in Tset, j in Bset
                        for t in keys(u_collection[t0][j])
                            u_collection[t0][j][t] *= scale_factor
                        end
                    end
                end
            end

            # Print progress
            show_this_iter = (k % progress_interval == 0) || (k == 1)
            if show_this_iter
                times = subproblem_times_k
                @printf "k=%3d  obj=\$%.2f  r=%.2e  s=%.2e  rho=%.1f  sp_max=%.3fs  sp_med=%.3fs\n" k true_objective r_norm s_norm rho_current maximum(times) median(times)
                if mod(k, 10) == 0
                    elapsed = time() - tadmm_start_time
                    @printf "      [Progress] %.1fs elapsed, %d iters/min\n" elapsed round(Int, k / (elapsed/60))
                end
            end

            # Convergence check
            stagnation_detected = false
            if length(obj_history) >= stagnation_window
                recent_objs = obj_history[end-stagnation_window+1:end]
                obj_min = minimum(recent_objs)
                obj_max = maximum(recent_objs)
                if obj_max > 0
                    improvement_ratio = (obj_max - obj_min) / abs(obj_max)
                    stagnation_detected = improvement_ratio < stagnation_threshold
                end
            end

            r_converged = r_norm <= eps_pri
            s_converged = s_norm <= eps_dual
            either_exit = s_converged || stagnation_detected

            if r_converged && either_exit
                converged = true
                print(COLOR_SUCCESS)
                if s_converged
                    @printf "tADMM converged at k=%d (r=%.2e <= %.2e, s=%.2e <= %.2e)\n" k r_norm eps_pri s_norm eps_dual
                else
                    @printf "tADMM converged at k=%d (r=%.2e <= %.2e, stagnation)\n" k r_norm eps_pri
                end
                @printf "  Final objective: \$%.2f\n" true_objective
                print(COLOR_RESET)
                break
            end

            final_iter = k
        end
    catch e
        if isa(e, InterruptException)
            @printf "\ntADMM interrupted at iteration %d\n" final_iter
        else
            rethrow(e)
        end
    end

    if !converged
        print(COLOR_WARNING)
        @printf "tADMM did NOT converge after %d iterations (r=%.2e, s=%.2e)\n" max_iter r_norm_history[end] s_norm_history[end]
        print(COLOR_RESET)
    end

    # Extract final solution (diagonal: subproblem t0 provides time t0)
    tadmm_wallclock_time = time() - tadmm_wallclock_start

    P_Subs_final = Dict(t => P_Subs_collection[t] for t in Tset)
    Q_Subs_final = Dict(t => Q_Subs_collection[t] for t in Tset)
    P_final = Dict((ij, t) => P_collection[ij][t] for ij in keys(P_collection), t in Tset)
    Q_final = Dict((ij, t) => Q_collection[ij][t] for ij in keys(Q_collection), t in Tset)
    l_final = Dict((ij, t) => l_collection[ij][t] for ij in keys(l_collection), t in Tset)
    v_final = Dict((j, t) => v_collection[j][t] for j in keys(v_collection), t in Tset)
    q_D_final = Dict((j, t) => q_D_collection[j][t] for j in keys(q_D_collection), t in Tset)
    P_B_final = Dict((j, t) => P_B_collection[j][t][t] for j in Bset, t in Tset)
    B_final = Dict((j, t) => B_collection[t][j][t] for j in Bset, t in Tset)

    # Print performance summary
    println("\n" * "="^80)
    println(COLOR_INFO, "SOLVER PERFORMANCE SUMMARY", COLOR_RESET)
    println("="^80)
    @printf "Iterations:       %d\n" final_iter
    total_effective = sum(iteration_effective_times)
    @printf "Effective time:   %.2f seconds\n" total_effective
    @printf "Wall-clock time:  %.2f seconds\n" tadmm_wallclock_time
    @printf "Overhead:         %.1f%%\n" 100*(tadmm_wallclock_time - total_effective)/tadmm_wallclock_time
    println("="^80)

    return Dict(
        :status => converged ? MOI.OPTIMAL : MOI.ITERATION_LIMIT,
        :P_Subs => P_Subs_final,
        :Q_Subs => Q_Subs_final,
        :P => P_final,
        :Q => Q_final,
        :v => v_final,
        :l => l_final,
        :q_D => q_D_final,
        :P_B => P_B_final,
        :B => B_final,
        :objective => last(obj_history),
        :convergence_history => Dict(
            :obj_history => obj_history,
            :r_norm_history => r_norm_history,
            :s_norm_history => s_norm_history,
            :rho_history => rho_history,
        ),
        :timing => Dict(
            :subproblem_times_history => subproblem_times_history,
            :iteration_effective_times => iteration_effective_times,
            :total_effective_time => total_effective,
            :total_sequential_time => sum(sum.(subproblem_times_history)),
            :wallclock_time => tadmm_wallclock_time,
            :subproblem_times_matrix => subproblem_times_matrix,
            :subproblem_iters_matrix => subproblem_iters_matrix,
        ),
        :subproblem_model_stats => tadmm_subproblem_stats,
    )
end

# ============================================================================
# SAVE tADMM RESULTS
# ============================================================================

function save_tadmm_results(sol, data)
    mkpath(SYSTEM_DIR)

    hist = sol[:convergence_history]
    obj_h = hist[:obj_history]
    r_h = hist[:r_norm_history]
    s_h = hist[:s_norm_history]
    rho_h = hist[:rho_history]
    eff_times = sol[:timing][:iteration_effective_times]
    n_iter = length(obj_h)

    # --- 1. Convergence CSV ---
    conv_csv = joinpath(SYSTEM_DIR, "convergence_data.csv")
    open(conv_csv, "w") do io
        println(io, "iteration,objective,r_norm,s_norm,rho,eff_time_this_k,cum_eff_time,max_sp_time_this_k,mean_ipopt_iters_this_k")
        cum_eff_time = 0.0
        for k in 1:n_iter
            rho_val = k <= length(rho_h) ? rho_h[k] : NaN
            eff_time_k = k <= length(eff_times) ? eff_times[k] : NaN
            cum_eff_time += isnan(eff_time_k) ? 0.0 : eff_time_k

            sp_times_k = sol[:timing][:subproblem_times_matrix][k, :]
            max_sp_time_k = maximum(sp_times_k)
            ipopt_iters_k = sol[:timing][:subproblem_iters_matrix][k, :]
            ipopt_nz = ipopt_iters_k[ipopt_iters_k .> 0]
            mean_ipopt_k = isempty(ipopt_nz) ? 0.0 : mean(ipopt_nz)

            @printf(io, "%d,%.6f,%.10e,%.10e,%.2f,%.4f,%.4f,%.4f,%.2f\n",
                    k, obj_h[k], r_h[k], s_h[k], rho_val, eff_time_k, cum_eff_time, max_sp_time_k, mean_ipopt_k)
        end
    end
    println(COLOR_SUCCESS, "Convergence CSV: $(conv_csv)", COLOR_RESET)

    # --- 2. Subproblem timing CSV ---
    csv_file = joinpath(SYSTEM_DIR, "subproblem_timing_details.csv")
    open(csv_file, "w") do io
        println(io, "iteration,subproblem,solve_time_sec,ipopt_iters")
        times_matrix = sol[:timing][:subproblem_times_matrix][1:n_iter, :]
        iters_matrix = sol[:timing][:subproblem_iters_matrix][1:n_iter, :]
        Tset_vec = collect(data[:Tset])
        for k in 1:n_iter
            for (idx, t0) in enumerate(Tset_vec)
                @printf(io, "%d,%d,%.6f,%d\n", k, t0, times_matrix[k, idx], iters_matrix[k, idx])
            end
        end
    end
    println(COLOR_SUCCESS, "Subproblem timing CSV: $(csv_file)", COLOR_RESET)

    # --- 3. Load BF results for comparison (if available) ---
    bf_objective = NaN
    bf_wallclock = NaN
    bf_solver_time = NaN
    bf_results_file = joinpath(SYSTEM_DIR, "results_socp_bf.txt")
    if isfile(bf_results_file)
        for line in eachline(bf_results_file)
            m = match(r"Total Cost: \$([0-9.]+)", line)
            if !isnothing(m)
                bf_objective = parse(Float64, m.captures[1])
            end
            m = match(r"Wall-clock time: ([0-9.]+) seconds", line)
            if !isnothing(m)
                bf_wallclock = parse(Float64, m.captures[1])
            end
            m = match(r"Solver time: ([0-9.]+) seconds", line)
            if !isnothing(m)
                bf_solver_time = parse(Float64, m.captures[1])
            end
        end
        if !isnan(bf_objective)
            println(COLOR_INFO, "Loaded BF results from $(bf_results_file)", COLOR_RESET)
        end
    else
        println(COLOR_WARNING, "BF results not found at $(bf_results_file), skipping comparison", COLOR_RESET)
    end

    # Validate solution
    tadmm_validation = validate_branch_flow_equations(sol, data; tol=1e-4, verbose=false)

    # --- 4. Results text file ---
    results_file = joinpath(SYSTEM_DIR, "results_socp_tadmm.txt")
    open(results_file, "w") do io
        println(io, "="^80)
        println(io, "SOCP tADMM OPTIMIZATION RESULTS")
        println(io, "="^80)
        println(io, "System: $(SYSTEM_NAME)")
        println(io, "Time horizon: T=$(T) periods")
        println(io, "Time step: $(DELTA_T_H) hours")
        println(io, "Number of buses: $(length(data[:Nset]))")
        println(io, "Number of branches: $(length(data[:Lset]))")
        println(io, "Number of batteries: $(length(data[:Bset]))")
        println(io, "Number of PV units: $(length(data[:Dset]))")
        if !isnothing(sol[:subproblem_model_stats])
            stats = sol[:subproblem_model_stats]
            println(io, "\n--- SUBPROBLEM SIZE (one time period) ---")
            @printf(io, "Number of variables: %d\n", stats[:n_variables])
            @printf(io, "Linear constraints: %d\n", stats[:n_linear_constraints])
            @printf(io, "Quadratic constraints (SOCP): %d\n", stats[:n_quadratic_constraints])
            @printf(io, "Nonlinear constraints: %d\n", stats[:n_nonlinear_constraints])
            @printf(io, "Total constraints: %d\n", stats[:total_constraints])
        end
        println(io, "\n--- tADMM PARAMETERS ---")
        @printf(io, "Initial rho: %.1f\n", rho_tadmm)
        println(io, "Adaptive rho: $(adaptive_rho_tadmm)")
        @printf(io, "Max iterations: %d\n", max_iter_tadmm)
        @printf(io, "Primal tolerance: %.1e\n", eps_pri_tadmm)
        @printf(io, "Dual tolerance: %.1e\n", eps_dual_tadmm)
        println(io, "\n--- CONVERGENCE ---")
        @printf(io, "Iterations: %d\n", n_iter)
        @printf(io, "Final primal residual: %.2e\n", r_h[end])
        @printf(io, "Final dual residual: %.2e\n", s_h[end])
        println(io, "\n--- OBJECTIVE VALUE ---")
        @printf(io, "Total Cost: \$%.4f\n", sol[:objective])
        if !isnan(bf_objective)
            obj_diff = abs(sol[:objective] - bf_objective)
            obj_rel = obj_diff / bf_objective * 100
            @printf(io, "Brute Force objective: \$%.4f\n", bf_objective)
            @printf(io, "Absolute difference: \$%.4f\n", obj_diff)
            @printf(io, "Relative difference: %.4f%%\n", obj_rel)
        end
        println(io, "\n--- COMPUTATION TIME ---")
        if !isnan(bf_wallclock)
            println(io, "Brute Force:")
            @printf(io, "  Wall-clock time: %.4f seconds\n", bf_wallclock)
            @printf(io, "  Solver time: %.4f seconds\n", bf_solver_time)
        end
        println(io, "tADMM:")
        @printf(io, "  Wall-clock time: %.4f seconds\n", sol[:timing][:wallclock_time])
        @printf(io, "  Effective time: %.4f seconds\n", sol[:timing][:total_effective_time])
        @printf(io, "  Sequential time: %.4f seconds\n", sol[:timing][:total_sequential_time])
        if !isnan(bf_wallclock)
            println(io, "Speedup Metrics:")
            @printf(io, "  tADMM vs BF (wall-clock): %.2fx\n", bf_wallclock / sol[:timing][:wallclock_time])
        end
        println(io, "\n--- SOLUTION FEASIBILITY ---")
        if tadmm_validation[:feasible]
            println(io, "Status: FEASIBLE - All constraints satisfied")
        else
            println(io, "Status: INFEASIBLE - $(tadmm_validation[:total_violations]) constraint violations")
            for (key, count) in tadmm_validation[:violations]
                if count > 0
                    max_viol = tadmm_validation[:max_violations][key]
                    @printf(io, "  %-30s: %6d violations (max: %.2e)\n", string(key), count, max_viol)
                end
            end
        end
        println(io, "\n--- SOLVER ---")
        solver_name = USE_GUROBI_FOR_TADMM ? "Gurobi" : "Ipopt"
        println(io, "Subproblem solver: $(solver_name)")
        println(io, "Formulation: SOCP (BFM-NL)")
        println(io, "Decomposition: Temporal ADMM")
        println(io, "="^80)
    end
    println(COLOR_SUCCESS, "Results: $(results_file)", COLOR_RESET)

    # --- 5. Serialized solution (.jls, gitignored) ---
    sol_file = joinpath(SYSTEM_DIR, "sol_socp_tadmm.jls")
    open(sol_file, "w") do io
        serialize(io, sol)
    end
    println(COLOR_SUCCESS, "Solution saved: $(sol_file)", COLOR_RESET)
end

# ============================================================================
# RUN
# ============================================================================

if isempty(data[:Bset])
    println(COLOR_WARNING, "No batteries in system — skipping tADMM", COLOR_RESET)
else
    println("\n" * "="^80)
    println(COLOR_HIGHLIGHT, "SOLVING MPOPF WITH SOCP (tADMM)", COLOR_RESET)
    println("="^80)

    if rho_scaling_with_T
        println(COLOR_INFO, "rho scaled with T: $(rho_base) * sqrt($(T)/24) = $(round(rho_tadmm, digits=1))", COLOR_RESET)
    else
        println(COLOR_INFO, "Using fixed rho = $(rho_tadmm)", COLOR_RESET)
    end

    solver_choice = USE_GUROBI_FOR_TADMM ? :gurobi : :ipopt
    sol_socp_tadmm = solve_MPOPF_SOCP_tADMM(data;
        rho=rho_tadmm,
        max_iter=max_iter_tadmm,
        eps_pri=eps_pri_tadmm,
        eps_dual=eps_dual_tadmm,
        adaptive_rho=adaptive_rho_tadmm,
        solver=solver_choice,
        enable_warm_start=enable_warm_start,
        track_subproblem_details=track_subproblem_details
    )

    # Report
    final_r = sol_socp_tadmm[:convergence_history][:r_norm_history][end]
    final_s = sol_socp_tadmm[:convergence_history][:s_norm_history][end]
    println("\n--- tADMM RESULT ---")
    print(COLOR_SUCCESS)
    @printf "Objective: \$%.2f\n" sol_socp_tadmm[:objective]
    @printf "Wall-clock: %.2f seconds\n" sol_socp_tadmm[:timing][:wallclock_time]
    @printf "Effective:  %.2f seconds (%d iterations)\n" sol_socp_tadmm[:timing][:total_effective_time] length(sol_socp_tadmm[:timing][:iteration_effective_times])
    print(COLOR_RESET)

    # Validate
    println("\n" * "="^80)
    println(COLOR_INFO, "VALIDATING tADMM SOLUTION", COLOR_RESET)
    println("="^80)
    tadmm_validation = validate_branch_flow_equations(sol_socp_tadmm, data; tol=1e-4, verbose=false)
    print_validation_summary(tadmm_validation; solution_name="tADMM SOCP")

    save_tadmm_results(sol_socp_tadmm, data)

    println("\n" * "="^80)
    println(COLOR_HIGHLIGHT, "tADMM PIPELINE COMPLETE", COLOR_RESET)
    println("="^80)
end
