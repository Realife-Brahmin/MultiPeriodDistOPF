module ModelCopier

export
    copy_modelVals,
    create_variable_dict,
    ModelVals,
    store_FS_t_k_decvar_values

using JuMP
using Ipopt
using Parameters

#region copy_modelVals
"""
    copy_modelVals(modelDict, model_Tset; Tset=nothing)

Copy the values of decision variables from the optimization model to modelVals.

This function extracts the values of decision variables from `model_Tset` and stores them in `modelVals` within `modelDict`.
It handles different sets of variables and updates the objective value, termination status, and solve time.
"""
function copy_modelVals(modelDict, model_Tset;
    Tset=nothing, solverCallType=nothing) # modelDict could be ddpModel or modelDict (temporallyBruteforced)

    @unpack modelVals, data = modelDict
    if Tset === nothing
        Tset = modelDict[:data][:Tset]
    end
    # Extract necessary sets from data
    @unpack Bset, Dset, Lset, Nset, linearizedModel, warmStart_mu = data

    if isnothing(solverCallType)
        solverCallType = "nonlinear"
    end

    # Retrieve variables from the model
    P_Subs_model = model_Tset[:P_Subs]
    P_model = model_Tset[:P]
    Q_model = model_Tset[:Q]
    v_model = model_Tset[:v]
    q_D_model = model_Tset[:q_D]
    q_B_model = model_Tset[:q_B]
    P_c_model = model_Tset[:P_c]
    P_d_model = model_Tset[:P_d]
    B_model = model_Tset[:B]

    # Store values into modelVals using the indices from data

    # P_Subs[t]
    for t in Tset
        modelVals[:P_Subs][t] = value(P_Subs_model[t])
    end

    # P[(i,j), t], Q[(i,j), t], l[(i,j), t] for (i,j) in Lset
    for (i, j) in Lset, t in Tset
        modelVals[:P][(i, j), t] = value(P_model[(i, j), t])
        modelVals[:Q][(i, j), t] = value(Q_model[(i, j), t])
        if !linearizedModel && solverCallType == "nonlinear" # neither warm-starting DDP nor doing linear model BF
            l_model = model_Tset[:l]
            modelVals[:l][(i, j), t] = value(l_model[(i, j), t])
        elseif linearizedModel || solverCallType == "linear" # either we're warm-starting DDP with a linear model or doing linear model BF
            modelVals[:l][(i, j), t] = (value(P_model[(i, j), t])^2 + value(Q_model[(i, j), t])^2) / value(v_model[i, t])
        else
            @error "Invalid solverCallType: $solverCallType"
            return
        end
    end

    # v[j, t] for j in Nset
    for j in Nset, t in Tset
        modelVals[:v][j, t] = value(v_model[j, t])
    end

    # q_D[j, t] for j in Dset
    for j in Dset, t in Tset
        modelVals[:q_D][j, t] = value(q_D_model[j, t])
    end

    # q_B[j, t], P_c[j, t], P_d[j, t], B[j, t] for j in Bset
    for j in Bset, t in Tset
        modelVals[:q_B][j, t] = value(q_B_model[j, t])
        modelVals[:P_c][j, t] = value(P_c_model[j, t])
        modelVals[:P_d][j, t] = value(P_d_model[j, t])
        modelVals[:B][j, t] = value(B_model[j, t])
    end

    @unpack T = data
    if length(Tset) == 1
        for t ∈ Tset
            modelVals[:objective_value_vs_t][t] = objective_value(model_Tset)
            modelVals[:termination_status_vs_t][t] = termination_status(model_Tset)
            modelVals[:solve_time_vs_t][t] = solve_time(model_Tset)
        end
    elseif length(Tset) == T
        modelVals[:objective_value] = objective_value(model_Tset)
        modelVals[:termination_status] = termination_status(model_Tset)
        modelVals[:solve_time] = solve_time(model_Tset)
    else
        @error "Invalid length of Tset: $(length(Tset))"
        return
    end

    @pack! modelDict = modelVals
    return modelDict
end
#endregion
#region create_variable_dict
"""
    create_variable_dict(model)

Create a dictionary of variable names and their values from the optimization model.

This function iterates over all variables in the given `model` and creates a dictionary mapping variable names to their values.
"""
function create_variable_dict(model)
    var_dict = Dict{Symbol,Float64}()
    for v in all_variables(model)
        var_name = Symbol(name(v))
        var_dict[var_name] = value(v)
    end
    return var_dict
end
#endregion

#region ModelVals
"""
    ModelVals(data)

Initialize a dictionary to store the values of decision variables.

This function creates and initializes a dictionary `modelVals` to store the values of various decision variables and other related information.
"""
function ModelVals(data)
    modelVals = Dict{Symbol,Any}()

    # Extract necessary sets from data
    Tset = data[:Tset]
    Bset = data[:Bset]
    Dset = data[:Dset]
    Lset = data[:Lset]
    Nset = data[:Nset]

    # Initialize containers for each variable
    modelVals[:P_Subs] = Dict{Int,Float64}()  # t => value
    modelVals[:P] = Dict{Tuple{Tuple{Int,Int},Int},Float64}()  # (i, j), t => value
    modelVals[:Q] = Dict{Tuple{Tuple{Int,Int},Int},Float64}()
    modelVals[:l] = Dict{Tuple{Tuple{Int,Int},Int},Float64}()
    modelVals[:v] = Dict{Tuple{Int,Int},Float64}() 
    modelVals[:q_D] = Dict{Tuple{Int,Int},Float64}()
    modelVals[:q_B] = Dict{Tuple{Int,Int},Float64}()
    modelVals[:P_c] = Dict{Tuple{Int,Int},Float64}()
    modelVals[:P_d] = Dict{Tuple{Int,Int},Float64}()
    modelVals[:B] = Dict{Tuple{Int,Int},Float64}()

    modelVals[:termination_status_vs_t] = Dict{Int,Any}()
    modelVals[:solve_time_vs_t] = Dict{Int,Float64}()
    modelVals[:objective_value_vs_t] = Dict{Int,Float64}()

    return modelVals
end
#endregion

#region store_FS_decvar_values_t_in_Tset
"""
    Store current forward step decision variable values
"""
function store_FS_t_k_decvar_values(ddpModel, Tset=nothing;
    optModel="Linear", verbose=false)
    
    if !optModel in ["Linear", "Nonlinear"]
        @error "Invalid optModel: $optModel. Must be 'Linear' or 'Nonlinear'."
        return
    elseif optModel == "Nonlinear"
        @error "Currently Nonlinear model is not supported."
        return
    elseif isnothing(Tset) || length(Tset) != 1
        @error "Tset must be a single time step."
        return
    end

    t_ddp = Tset[1]
    @unpack data, k_ddp = ddpModel

    @unpack modelVals_ddp_vs_t_vs_k, models_ddp_vs_t_vs_k = ddpModel
    modelVals_t0_k0 = modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp]
    model_t0_k0 = models_ddp_vs_t_vs_k[t_ddp, k_ddp]

    # The actual copying

    # Copy P_Subs[t]
    modelVals_t0_k0[:P_Subs][t_ddp] = JuMP.value(model_t0_k0[:P_Subs][t_ddp])

    # Copy P[(i,j), t], Q[(i,j), t] for (i,j) in Lset
    @unpack Lset = data
    for (i, j) in Lset
        modelVals_t0_k0[:P][(i, j), t_ddp] = JuMP.value(model_t0_k0[:P][(i, j), t_ddp])
        modelVals_t0_k0[:Q][(i, j), t_ddp] = JuMP.value(model_t0_k0[:Q][(i, j), t_ddp])
    end

    # Copy v[j, t] for j in Nset
    @unpack Nset = data
    for j in Nset
        modelVals_t0_k0[:v][j, t_ddp] = JuMP.value(model_t0_k0[:v][j, t_ddp])
    end

    if optModel == "Nonlinear"
        # Copy l[(i,j), t] for (i,j) in Lset
        for (i, j) in Lset
            modelVals_t0_k0[:l][(i, j), t_ddp] = JuMP.value(model_t0_k0[:l][(i, j), t_ddp])
        end
    elseif optModel == "Linear"
        # Compute (approximately) and copy l[(i,j), t] for (i,j) in Lset
        for (i, j) in Lset
            modelVals_t0_k0[:l][(i, j), t_ddp] = (modelVals_t0_k0[:P][(i, j), t_ddp]^2 + modelVals_t0_k0[:Q][(i, j), t_ddp]^2) / modelVals_t0_k0[:v][i, t_ddp]
        end
    else
        @error "Invalid optModel: $optModel. Must be 'Linear' or 'Nonlinear'."
        return
    end

    # Copy q_D[j, t] for j in Dset
    @unpack Dset = data
    for j in Dset
        modelVals_t0_k0[:q_D][j, t_ddp] = JuMP.value(model_t0_k0[:q_D][j, t_ddp])
    end

    # Copy q_B[j, t], P_c[j, t], P_d[j, t], B[j, t] for j in Bset
    @unpack Bset = data
    for j in Bset
        modelVals_t0_k0[:q_B][j, t_ddp] = JuMP.value(model_t0_k0[:q_B][j, t_ddp])
        modelVals_t0_k0[:P_c][j, t_ddp] = JuMP.value(model_t0_k0[:P_c][j, t_ddp])
        modelVals_t0_k0[:P_d][j, t_ddp] = JuMP.value(model_t0_k0[:P_d][j, t_ddp])
        modelVals_t0_k0[:B][j, t_ddp] = JuMP.value(model_t0_k0[:B][j, t_ddp])
    end

    # Copy termination status, solve time, and objective value
    modelVals_t0_k0[:termination_status] = JuMP.termination_status(model_t0_k0)
    modelVals_t0_k0[:solve_time] = JuMP.solve_time(model_t0_k0)
    modelVals_t0_k0[:objective_value] = JuMP.objective_value(model_t0_k0)

    modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp] = modelVals_t0_k0
    @pack! ddpModel = modelVals_ddp_vs_t_vs_k
    return ddpModel
end
#endregion

end