# Import JuMP and solver
import JuMP
import Gurobi
import Ipopt

# solve_MPOPF_with_LDF_BruteForced.jl
"""
    solve_MPOPF_with_LDF_BruteForced(data; Tset=nothing, verbose=false)

Solve the Multi-Period Optimal Power Flow (MPOPF) problem using a brute-force approach with the LinDistFlow (LDF) model.
This function is self-contained and does not depend on ModelBuilder or other modules. It builds and solves the optimization
problem for the given data and returns the solution.

# Arguments
- `data::Dict`: Dictionary containing all system and simulation data.
- `Tset::Union{Nothing, Vector{Int}}`: Optional vector of time steps to consider. Defaults to all time steps in data.
- `verbose::Bool`: Print solver output if true.

# Returns
- `result::Dict`: Dictionary containing the solution (optimal values for variables, objective value, status, etc.)

# Notes
- This function assumes a single-phase, linearized (LDF) model.
- The optimization is performed for all time steps in Tset (brute-forced, not decomposed).
- Only standard JuMP and solver packages are used.
"""
function solve_MPOPF_with_LDF_BruteForced(data; Tset=nothing, verbose=false, solver=:gurobi)

    # Get time set
    if Tset === nothing
        Tset = data[:Tset]
    end
    N = data[:N]
    m = data[:m]
    NLset = data[:NLset]
    Dset = data[:Dset]
    Bset = data[:Bset]
    Lset = data[:Lset]
    parent = data[:parent]
    children = data[:children]
    p_L = data[:p_L]
    q_L = data[:q_L]
    p_D = data[:p_D]
    B0 = data[:B0]
    B_R = data[:B_R]
    P_B_R = data[:P_B_R]
    eta_C = data[:eta_C]
    eta_D = data[:eta_D]
    soc_min = data[:soc_min]
    soc_max = data[:soc_max]
    rdict = data[:rdict]
    xdict = data[:xdict]
    Vminpu = data[:Vminpu]
    Vmaxpu = data[:Vmaxpu]
    substationBus = data[:substationBus]
    delta_t = data[:delta_t]

    # Create JuMP model
    if solver == :gurobi
        model = JuMP.Model(Gurobi.Optimizer)
    else
        model = JuMP.Model(Ipopt.Optimizer)
    end
    JuMP.set_silent(model)
    if verbose
        JuMP.set_silent(model, false)
    end

    # Variables
    @JuMP.variable(model, V[bus in 1:N, t in Tset] >= 0)
    @JuMP.variable(model, P[branch in Lset, t in Tset])
    @JuMP.variable(model, Q[branch in Lset, t in Tset])
    @JuMP.variable(model, P_B[bus in Bset, t in Tset])
    @JuMP.variable(model, soc[bus in Bset, t in Tset] >= 0)


    # Initial SOC
    for bus in Bset
        JuMP.@constraint(model, soc[bus, Tset[1]] == B0[bus])
    end

    # SOC dynamics (P_B: positive = discharge, negative = charge)
    for bus in Bset, (i, t) in enumerate(Tset[2:end])
        prev_t = Tset[i]
        JuMP.@constraint(model, soc[bus, t] == soc[bus, prev_t] - P_B[bus, prev_t]*delta_t)
    end

    # SOC limits
    for bus in Bset, t in Tset
        JuMP.@constraint(model, soc_min[bus]*B_R[bus] <= soc[bus, t] <= soc_max[bus]*B_R[bus])
    end

    # Battery power limits
    for bus in Bset, t in Tset
        JuMP.@constraint(model, -P_B_R[bus] <= P_B[bus, t] <= P_B_R[bus])
    end

    # LinDistFlow constraints (simplified)
    for (i, j) in Lset, t in Tset
        # Power flow equations
        JuMP.@constraint(model, P[(i,j), t] == sum(p_L[(k, t)] for k in children[j] if k in NLset) + sum(p_D[(k, t)] for k in children[j] if k in Dset) + sum(P_B[k, t] for k in children[j] if k in Bset) + sum(P[(j, k2), t] for k2 in children[j] if (j, k2) in Lset))
        JuMP.@constraint(model, Q[(i,j), t] == sum(q_L[(k, t)] for k in children[j] if k in NLset) + sum(Q[(j, k2), t] for k2 in children[j] if (j, k2) in Lset))
        # Voltage drop
        JuMP.@constraint(model, V[j, t] == V[i, t] - 2*rdict[(i,j)]*P[(i,j), t] - 2*xdict[(i,j)]*Q[(i,j), t])
    end

    # Voltage limits
    for bus in 1:N, t in Tset
        JuMP.@constraint(model, Vminpu[bus] <= V[bus, t] <= Vmaxpu[bus])
    end

    # Substation voltage fixed
    for t in Tset
        JuMP.@constraint(model, V[substationBus, t] == 1.0)
    end

    # Objective: minimize total substation power cost + battery quadratic cost over all time
    LoadShapeCost = get(data, :LoadShapeCost, nothing)
    C_B = get(data, :C_B, 0.0)
    obj = 0.0
    for t in Tset
        # Substation power cost (if available)
        if LoadShapeCost !== nothing
            for (substationBus, k) in Lset
                if substationBus == 1
                    obj += LoadShapeCost[t] * P[(substationBus, k), t] * delta_t
                end
            end
        else
            for (substationBus, k) in Lset
                if substationBus == 1
                    obj += P[(substationBus, k), t]
                end
            end
        end
        # Battery quadratic cost
        for bus in Bset
            obj += C_B * (P_B[bus, t])^2 * delta_t
        end
    end
    JuMP.@objective(model, Min, obj)

    JuMP.optimize!(model)

    status = JuMP.termination_status(model)
    result = Dict(
        :status => status,
        :objective => JuMP.objective_value(model),
        :V => JuMP.value.(V),
        :P => JuMP.value.(P),
        :Q => JuMP.value.(Q),
        :P_B => JuMP.value.(P_B),
        :soc => JuMP.value.(soc)
    )
    return result
end
