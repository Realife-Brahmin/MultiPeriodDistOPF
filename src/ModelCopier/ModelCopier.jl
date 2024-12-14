module ModelCopier

export
    copy_modelVals,
    create_variable_dict,
    ModelVals

using JuMP
using Ipopt
using Parameters

function copy_modelVals(modelDict, model_Tset;
    Tset=nothing) # modelDict could be ddpModel or modelDict (temporallyBruteforced)

    @unpack modelVals, data = modelDict;
    if Tset === nothing
        Tset = modelDict[:data][:Tset]
    end
    # Extract necessary sets from data
    @unpack Bset, Dset, Lset, Nset = data;

    # Retrieve variables from the model
    P_Subs_model = model_Tset[:P_Subs]
    P_model = model_Tset[:P]
    Q_model = model_Tset[:Q]
    l_model = model_Tset[:l]
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
        modelVals[:l][(i, j), t] = value(l_model[(i, j), t])
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

    @unpack T = data;
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
    
    @pack! modelDict = modelVals;
    return modelDict
end

function create_variable_dict(model)
    var_dict = Dict{Symbol,Float64}()
    for v in all_variables(model)
        var_name = Symbol(name(v))
        var_dict[var_name] = value(v)
    end
    return var_dict
end

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

end