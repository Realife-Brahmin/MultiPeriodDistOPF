# optimizer.jl
module Playbook_of_MPOPF

export
    get_soc_dual_variables_fullMPOPF,
    optimize_MPOPF_1ph_NL_TemporallyBruteforced,
    print_mu

include("./ModelBuilder/ModelBuilder.jl")
import .ModelBuilder as MB

include("./ModelCopier/ModelCopier.jl")
import .ModelCopier as MC

include("./helperFunctions.jl")
using .helperFunctions

include("./DDP/DDP.jl")
using .DDP

include("./SolverArranger/SolverArranger.jl")
import .SolverArranger as SolverArranger

include("./Exporter.jl")
import .Exporter as Exporter

using Crayons
using JuMP
using EAGO
using Gurobi
using Ipopt
using Juniper
using MadNLP
using Parameters: @unpack, @pack!

function optimize_MPOPF_1ph_NL_TemporallyBruteforced(data)

    Tset = data[:Tset] # In this case, Tset = [1, 2, ... T]
    modelDict = MB.build_MPOPF_1ph_NL_model_t_in_Tset(data, Tset=Tset)

    @unpack model, data = modelDict
    optimize!(model)
    
    @pack! modelDict = model

    # modelDict = generate_1ph_NL_model_decvar_value_dict(modelDict)
    modelDict = MC.copy_modelVals(modelDict, model, Tset=Tset)

    # Extract SOC dual variables and update modelDict
    mu = get_soc_dual_variables_fullMPOPF(modelDict, Tset=Tset)
    modelDict[:mu] = mu

    # Check solver status and retrieve results
    # Define crayons for green and red text
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    crayon_red = Crayon(foreground=:red, bold=true)

    @unpack modelVals = modelDict
    termination_status = modelVals[:termination_status]

    crayon_final_green = Crayon(foreground=:white,background=:green, bold=true)
    if termination_status == LOCALLY_SOLVED
        println(crayon_final_green("Optimal solution found."))
    else
        println(crayon_red("Optimization did not find an optimal solution."))
    end

    # optimal_obj_value = objective_value(model)
    optimal_obj_value = modelVals[:objective_value]
    println("Optimal objective function value: ", optimal_obj_value)
    
    Exporter.export_optimization_model(modelDict, verbose=false)
    
    return modelDict

end

function get_soc_dual_variables_fullMPOPF(modelDict; Tset=nothing)
    @unpack model, data = modelDict
    if Tset === nothing
        Tset = data[:Tset]
    end
    @unpack Bset, T = data
    
    mu = Dict{Tuple{Int,Int},Float64}()
    for t in Tset
        for j in Bset
            if t == 1
                constraint_name = "h_SOC_j^{t=1}_Initial_SOC_Node_j_$(j)_t1"
            elseif 2 <= t <= T
                constraint_name = "h_SOC_j^{t=2toT}_SOC_Trajectory_Node_j_$(j)_t_$(t)"
            else
                @error "Invalid value of t: $t"
                continue
            end
            constraint_j_t = constraint_by_name(model, constraint_name)
            mu[(j, t)] = dual(constraint_j_t)
        end
    end
    return mu
end

function print_mu(modelDict)
    crayon_header = Crayon(foreground=:white, background=:blue, bold=true)
    crayon_value = Crayon(foreground=:light_green, bold=true)
    crayon_error = Crayon(foreground=:red, bold=true)

    println(crayon_header("Dual Variables (mu) for SOC Constraints:"))

    @unpack mu, data = modelDict
    @unpack Tset, Bset, temporal_decmp = data

    # Limit to the first 2 batteries if there are more than 2
    # Bset_to_print = Bset
    Bset_to_print = length(Bset) > 2 ? [Bset[1], Bset[end]] : Bset

    if temporal_decmp == false
        for t in Tset
            for j in Bset_to_print
                if haskey(mu, (j, t))
                    println(crayon_value("mu[$j, $t] = $(mu[(j, t)])"))
                else
                    println(crayon_error("mu[$j, $t] not found"))
                end
            end
        end
    elseif temporal_decmp == true
        @unpack k_ddp = modelDict
        for t in Tset
            for j in Bset_to_print
                if haskey(mu, (j, t, k_ddp - 1))
                    println(crayon_value("mu[$j, $t, $(k_ddp-1)] = $(mu[(j, t, k_ddp-1)])"))
                else
                    println(crayon_error("mu[$j, $t, $(k_ddp-1)] not found"))
                end
            end
        end
    else
        println(crayon_error("Invalid value for temporal_decmp: $temporal_decmp"))
    end
end

end # module Playbook_of_MPOPF
