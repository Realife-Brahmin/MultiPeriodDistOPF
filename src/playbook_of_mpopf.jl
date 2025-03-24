# optimizer.jl
module Playbook_of_MPOPF

export
    get_soc_dual_variables_fullMPOPF,
    optimize_MPOPF_1ph_NL_TemporallyBruteforced,
    optimize_MPOPF_1ph_L,
    get_dual_variables_state_fullMPOPF

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

    # Exporter.export_optimization_model(modelDict, verbose=false)

    @unpack model, data = modelDict
    set_optimizer_attribute(model, "print_level", 3)
    # set_optimizer_attribute(model, "print_frequency_iter", 100)


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
    println("Optimal objective function value:  $(optimal_obj_value)")
        
    return modelDict

end

function optimize_MPOPF_1ph_L(data)

    Tset = data[:Tset] # In this case, Tset = [1, 2, ... T]
    modelDict = MB.build_MPOPF_1ph_L_model_t_in_Tset(data, Tset=Tset)

    @unpack model, data = modelDict
    optimize!(model)

    @pack! modelDict = model

    # modelDict = generate_1ph_L_model_decvar_value_dict(modelDict)
    modelDict = MC.copy_modelVals(modelDict, model, Tset=Tset, solverCallType="linear")

    # Extract SOC dual variables and update modelDict
    mu = get_soc_dual_variables_fullMPOPF(modelDict, Tset=Tset)
    modelDict[:mu] = mu

    # Check solver status and retrieve results
    # Define crayons for green and red text
    crayon_light_green = Crayon(foreground=:light_green, bold=true)
    crayon_red = Crayon(foreground=:red, bold=true)

    @unpack modelVals = modelDict
    termination_status = modelVals[:termination_status]

    crayon_final_green = Crayon(foreground=:white, background=:green, bold=true)
    if termination_status == LOCALLY_SOLVED
        println(crayon_final_green("Optimal solution found."))
    else
        println(crayon_red("Optimization did not find an optimal solution."))
    end

    # optimal_obj_value = objective_value(model)
    optimal_obj_value = modelVals[:objective_value]
    println("Optimal objective function value: $(optimal_obj_value)")

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
            mu[(j, t)] = -dual(constraint_j_t)
        end
    end
    return mu
end

function get_dual_variables_state_fullMPOPF(ddpModel;
    tolKKT=1e-6,
    verbose=false)

    kktBalanceDict = Dict()
    kktMetricDict = 0.0
    muDict = Dict()
    lambda_lb_Dict = Dict()
    lambda_ub_Dict = Dict()

    crayon_header = Crayon(foreground=:white, background=:blue, bold=true)
    crayon_error = Crayon(foreground=:red, bold=true)

    # Define colors for each battery for non-temporal decomposition
    battery_colors_non_temporal = Dict(
        1 => Crayon(foreground=:red, bold=true),
        2 => Crayon(foreground=:green, bold=true),
        3 => Crayon(foreground=:blue, bold=true)
    )

    # Define negative colors for each battery for temporal decomposition
    battery_colors_temporal = Dict(
        1 => Crayon(foreground=:light_red, bold=true),
        2 => Crayon(foreground=:light_green, bold=true),
        3 => Crayon(foreground=:light_blue, bold=true)
    )

    myprintln(verbose, crayon_header("Dual Variables (mu) for SOC Constraints:"))

    @unpack mu, data = ddpModel
    @unpack Tset, Bset, temporal_decmp = data

    # Limit to the first and last batteries if there are more than 2
    Bset_to_print = length(Bset) > 2 ? [Bset[1], Bset[end]] : Bset

    if temporal_decmp == false
        @unpack model = ddpModel
        for (j_B, j) in enumerate(Bset_to_print)
            battery_color = battery_colors_non_temporal[j_B]
            for t in Tset
                if haskey(mu, (j, t))
                    muDict[(j, t)] = mu[(j, t)]
                    mu_j_t_str = trim_number_for_printing(mu[(j, t)], sigdigits=4)
                    myprintln(verbose, battery_color("μ[$j, $t] = $(mu_j_t_str)"))
                else
                    myprintln(verbose, crayon_error("μ[$j, $t] not found"))
                end
            end
        end
    elseif temporal_decmp == true
        @unpack models_ddp_vs_t_vs_k, k_ddp = ddpModel
        for (j_B, j) in enumerate(Bset_to_print)
            battery_color = battery_colors_temporal[j_B]
            for t in Tset
                if haskey(mu, (j, t, k_ddp - 1))
                    muDict[(j, t)] = mu[(j, t, k_ddp - 1)]
                    mu_j_t_k_str = trim_number_for_printing(mu[(j, t, k_ddp - 1)], sigdigits=4)
                    myprintln(verbose, battery_color("μ[$j, $t, $(k_ddp-1)] = $(mu_j_t_k_str)"))
                else
                    myprintln(verbose, crayon_error("μ[$j, $t, $(k_ddp-1)] not found"))
                end
            end
        end
    end

    myprintln(verbose, crayon_header("Dual Variables (lambda) for SOC Limits:"))

    for (j_B, j) in enumerate(Bset_to_print)
        battery_color = temporal_decmp ? battery_colors_temporal[j_B] : battery_colors_non_temporal[j_B]
        for t in Tset
            if temporal_decmp == false
                @unpack model = ddpModel
            else
                model = models_ddp_vs_t_vs_k[t, k_ddp-1]
            end
            lambda_lower_name = "g_11_j^t_MinSOC_Node_j_$(j)_t_$(t)"
            lambda_upper_name = "g_12_j^t_MaxSOC_Node_j_$(j)_t_$(t)"
            lambda_lower = -dual(constraint_by_name(model, lambda_lower_name))
            lambda_lb_Dict[(j, t)] = lambda_lower
            lambda_lower_str = trim_number_for_printing(lambda_lower, sigdigits=4)
            lambda_upper = -dual(constraint_by_name(model, lambda_upper_name))
            lambda_ub_Dict[(j, t)] = lambda_upper
            lambda_upper_str = trim_number_for_printing(lambda_upper, sigdigits=4)
            myprintln(verbose, battery_color("λ_lb[$j, $t] = $(lambda_lower_str) | λ_ub[$j, $t] = $(lambda_upper_str)"))
        end
    end

    myprintln(verbose, crayon_header("Checking for KKT Necessary Condition (∇L_{B_j^t} = 0 ∀ j ∈ B, t ∈ τ):"))

    for (j_B, j) in enumerate(Bset_to_print)
        battery_color = temporal_decmp ? battery_colors_temporal[j_B] : battery_colors_non_temporal[j_B]
        kkt_balance_total = 0.0
        for t in Tset
            if temporal_decmp == false
                @unpack model = ddpModel
            else
                model = models_ddp_vs_t_vs_k[t, k_ddp-1]
            end
            lambda_lower_name = "g_11_j^t_MinSOC_Node_j_$(j)_t_$(t)"
            lambda_upper_name = "g_12_j^t_MaxSOC_Node_j_$(j)_t_$(t)"
            lambda_lower = -dual(constraint_by_name(model, lambda_lower_name))
            lambda_upper = -dual(constraint_by_name(model, lambda_upper_name))
            if temporal_decmp == false
                mu_current = mu[(j, t)]
                if t < maximum(Tset)
                    mu_next = mu[(j, t + 1)]
                    balance = -lambda_lower + lambda_upper + mu_current - mu_next
                else
                    balance = -lambda_lower + lambda_upper + mu_current
                end
                kktBalanceDict[(j, t)] = balance
            else
                mu_current = mu[(j, t, k_ddp - 1)]
                if t < maximum(Tset)
                    mu_next = mu[(j, t + 1, k_ddp - 1)]
                    balance = -lambda_lower + lambda_upper + mu_current - mu_next
                else
                    balance = -lambda_lower + lambda_upper + mu_current
                end
                kktBalanceDict[(j, t)] = balance
            end
            kkt_balance_total += abs(balance)
            balance_str = trim_number_for_printing(balance, sigdigits=4)
            myprintln(verbose, battery_color("∇L_{B_j^t} for [$j, $t]: $balance_str"))
            kktMetricDict = balance
        end
        @unpack T = data
        kkt_balance_avg_str = trim_number_for_printing(kkt_balance_total/T, sigdigits=4)
        myprintln(verbose, "**************")
        myprintln(verbose, battery_color("Total KKT balance for B_$j: $(kkt_balance_avg_str)"))
        myprintln(verbose, "**************")

        dualStateDict = Dict(
            :muDict => muDict,
            :lambda_lb_Dict => lambda_lb_Dict,
            :lambda_ub_Dict => lambda_ub_Dict,
            :kktBalanceDict => kktBalanceDict,
            :kktMetricDict => kktMetricDict
        )

        return dualStateDict
    end
end

end # module Playbook_of_MPOPF
