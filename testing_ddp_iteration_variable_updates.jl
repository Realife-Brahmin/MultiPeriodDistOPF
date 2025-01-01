# everything is already loaded, just create the for loop of comparison

@unpack data = modelDict
@unpack temporal_decmp = data

if !temporal_decmp
    @error "Cannot test for DDP as it wasn't run"
    return
end

ddpModel = modelDict
@unpack data, modelVals, mu,modelVals_ddp_vs_t_vs_k = ddpModel
@unpack Bset, Compset = data

using Crayons
crayon_light_green = Crayon(foreground=:light_green, bold=true)
crayon_red = Crayon(foreground=:red, bold=true)

μ = mu

# @unpack Tset = data;
# for t_ddp ∈ Tset
#     println("t_ddp = $t_ddp")
#     for j in Bset
#         values_mu = []
#         values_B = []
#         for k_ddp in 1:3
#             B = modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp][:B]
#             push!(values_B, trim_number_for_printing(B[j, t_ddp], sigdigits=2))
#             push!(values_mu, trim_number_for_printing(μ[j, t_ddp, k_ddp], sigdigits=2))
#         end
#         println(crayon_light_green("j = $j: μ = ", join(values_mu, ", ")))
#         println("j = $j: B = ", join(values_B, ", "))
#     end
# end

using Crayons

@unpack modelVals_ddp_vs_t_vs_k = ddpModel;

crayon_red = Crayon(foreground=:red, bold=true)
crayon_red_negative = Crayon(foreground=:red, bold=true, negative=true)
crayon_light_green = Crayon(foreground=:light_green, bold=true)
crayon_blue = Crayon(foreground=:blue, bold=true)
crayon_magenta = Crayon(foreground=:magenta, bold=true)
crayon_cyan = Crayon(foreground=:cyan, bold=true)
crayon_yellow = Crayon(foreground=:yellow, bold=true)
crayon_light_red = Crayon(foreground=:light_red, bold=true)
crayon_light_blue = Crayon(foreground=:light_blue, bold=true)
crayon_light_magenta = Crayon(foreground=:light_magenta, bold=true)
crayon_light_cyan = Crayon(foreground=:light_cyan, bold=true)
crayon_light_yellow = Crayon(foreground=:light_yellow, bold=true)

@unpack Tset, Bset, Nset, Lset, Compset, Dset = data;

# Compare variable values for all t and k
for t_ddp in Tset
    println(crayon_red_negative("Comparing variable values for t_ddp = $t_ddp"))

    # Compare B variables
    for j in Bset
        values_B = []
        for k_ddp in 1:3
            B = modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp][:B]
            push!(values_B, trim_number_for_printing(B[j, t_ddp], sigdigits=2))
        end
        println(crayon_light_green("B[$j, $t_ddp] = ", join(values_B, ", ")))
    end

    # Compare PSubs variables
    values_P_Subs = []
    for k_ddp in 1:3
        P_Subs = modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp][:P_Subs]
        push!(values_P_Subs, trim_number_for_printing(P_Subs[t_ddp], sigdigits=2))
    end
    println(crayon_blue("P_Subs[$t_ddp] = ", join(values_P_Subs, ", ")))

    # Compare P and Q variables
    for (i, j) in Lset
        values_P = []
        values_Q = []
        for k_ddp in 1:3
            P = modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp][:P]
            Q = modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp][:Q]
            push!(values_P, trim_number_for_printing(P[(i, j), t_ddp], sigdigits=4))
            push!(values_Q, trim_number_for_printing(Q[(i, j), t_ddp], sigdigits=4))
        end
        println(crayon_magenta("P[($i, $j), $t_ddp] = ", join(values_P, ", ")))
        println(crayon_cyan("Q[($i, $j), $t_ddp] = ", join(values_Q, ", ")))
    end

    # Compare v variables
    for j in Compset
        values_v = []
        for k_ddp in 1:3
            v = modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp][:v]
            push!(values_v, trim_number_for_printing(v[j, t_ddp], sigdigits=4))
        end
        println(crayon_yellow("v[$j, $t_ddp] = ", join(values_v, ", ")))
    end

    # Compare q_D variables
    for j in Dset
        values_q_D = []
        for k_ddp in 1:3
            q_D = modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp][:q_D]
            push!(values_q_D, trim_number_for_printing(q_D[j, t_ddp], sigdigits=2))
        end
        println(crayon_light_red("q_D[$j, $t_ddp] = ", join(values_q_D, ", ")))
    end

    # Compare q_B variables
    for j in Bset
        values_q_B = []
        for k_ddp in 1:3
            q_B = modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp][:q_B]
            push!(values_q_B, trim_number_for_printing(q_B[j, t_ddp], sigdigits=2))
        end
        println(crayon_light_blue("q_B[$j, $t_ddp] = ", join(values_q_B, ", ")))
    end

    # Compare P_c variables
    for j in Bset
        values_P_c = []
        for k_ddp in 1:3
            P_c = modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp][:P_c]
            push!(values_P_c, trim_number_for_printing(P_c[j, t_ddp], sigdigits=2))
        end
        println(crayon_light_magenta("P_c[$j, $t_ddp] = ", join(values_P_c, ", ")))
    end

    # Compare P_d variables
    for j in Bset
        values_P_d = []
        for k_ddp in 1:3
            P_d = modelVals_ddp_vs_t_vs_k[t_ddp, k_ddp][:P_d]
            push!(values_P_d, trim_number_for_printing(P_d[j, t_ddp], sigdigits=2))
        end
        println(crayon_light_cyan("P_d[$j, $t_ddp] = ", join(values_P_d, ", ")))
    end
end