# Build the centralized loss-aware BFM-SOCP model from already-exported network
# data and ask Ipopt to report its global NLP/KKT statistics. This is a short
# diagnostic, not a benchmark driver.

using Ipopt
using JuMP
using Serialization

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const MOI = JuMP.MOI

system = length(ARGS) >= 1 ? ARGS[1] : "ieee123C_1ph"
T = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3
datafile = joinpath(REPO, "ddp", "results", "network_filterddp",
                    "network_data_$(system)_T$(T).jls")
data = deserialize(datafile)

Nset, Lset, Bset, Dset = data[:Nset], data[:Lset], data[:Bset], data[:Dset]
Tset = 1:T
root, nonroot = data[:substationBus], data[:Nm1set]
dt, pbase = data[:delta_t_h], data[:kVA_B]

model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 5)
set_optimizer_attribute(model, "print_timing_statistics", "yes")
set_optimizer_attribute(model, "output_file", joinpath(REPO, "ddp", "results",
    "network_filterddp", "$(system)_ipopt_T$(T)_matrix.log"))

@variable(model, P_Subs[Tset] >= 0)
@variable(model, Q_Subs[Tset])
@variable(model, P[Lset, Tset])
@variable(model, Q[Lset, Tset])
@variable(model, v[Nset, Tset])
@variable(model, ell[Lset, Tset] >= 0)
@variable(model, P_B[Bset, Tset])
@variable(model, B[Bset, Tset])
@variable(model, q_D[Dset, Tset])

@objective(model, Min,
    sum(data[:LoadShapeCost][t] * pbase * dt * P_Subs[t] for t in Tset) +
    sum(data[:C_B] * pbase^2 * dt * P_B[j,t]^2 for j in Bset, t in Tset))

for t in Tset
    @constraint(model, P_Subs[t] == sum(P[e,t] for e in data[:L1set]))
    @constraint(model, Q_Subs[t] == sum(Q[e,t] for e in data[:L1set]))
    for j in nonroot
        incoming = (data[:parent][j], j)
        outgoing_p = sum(P[(j,k),t] for k in data[:children][j]; init=0.0)
        outgoing_q = sum(Q[(j,k),t] for k in data[:children][j]; init=0.0)
        pL = j in data[:NLset] ? data[:p_L_pu][j,t] : 0.0
        qL = j in data[:NLset] ? data[:q_L_pu][j,t] : 0.0
        pD = j in Dset ? data[:p_D_pu][j,t] : 0.0
        pb = j in Bset ? P_B[j,t] : 0.0
        qd = j in Dset ? q_D[j,t] : 0.0
        @constraint(model, outgoing_p - P[incoming,t] +
            data[:rdict_pu][incoming] * ell[incoming,t] == pb + pD - pL)
        @constraint(model, outgoing_q - Q[incoming,t] +
            data[:xdict_pu][incoming] * ell[incoming,t] == qd - qL)
    end
    for e in Lset
        i, j = e
        r, x = data[:rdict_pu][e], data[:xdict_pu][e]
        @constraint(model, v[j,t] == v[i,t] - 2(r*P[e,t] + x*Q[e,t]) +
            (r^2+x^2)*ell[e,t])
        @constraint(model, P[e,t]^2 + Q[e,t]^2 <= v[i,t]*ell[e,t])
    end
    @constraint(model, v[root,t] == 1.05^2)
    for j in Nset
        set_lower_bound(v[j,t], data[:Vminpu][j]^2)
        set_upper_bound(v[j,t], data[:Vmaxpu][j]^2)
    end
    for j in Dset
        qmax = sqrt(max(0.0, data[:S_D_R][j]^2 - data[:p_D_pu][j,t]^2))
        set_lower_bound(q_D[j,t], -qmax)
        set_upper_bound(q_D[j,t], qmax)
    end
    for j in Bset
        set_lower_bound(P_B[j,t], -data[:P_B_R_pu][j])
        set_upper_bound(P_B[j,t], data[:P_B_R_pu][j])
        set_lower_bound(B[j,t], data[:soc_min][j] * data[:B_R_pu][j])
        set_upper_bound(B[j,t], data[:soc_max][j] * data[:B_R_pu][j])
        if t == 1
            @constraint(model, B[j,t] == data[:B0_pu][j] - dt*P_B[j,t])
        else
            @constraint(model, B[j,t] == B[j,t-1] - dt*P_B[j,t])
        end
    end
end

println("IPOPT_DIAGNOSTIC system=$system T=$T JuMP_variables=$(num_variables(model))")
optimize!(model)
println("IPOPT_DIAGNOSTIC status=$(termination_status(model)) objective=$(objective_value(model))")
println("IPOPT_DIAGNOSTIC solver=$(MOI.get(model, MOI.SolverName()))")
