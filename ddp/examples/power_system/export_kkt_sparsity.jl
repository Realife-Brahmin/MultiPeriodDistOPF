# Export the KKT sparsity pattern of the centralized MPOPF for spy plots.
#
# Builds the JuMP/Ipopt model (same as centralized_ipopt_matched.jl) WITHOUT
# solving, extracts the Jacobian and Hessian sparsity via MOI, assembles the
# full KKT pattern, and writes CSV files for external plotting.
#
#   PROFILE_PERIODIC=1 REDUCED_CB=1e-3 TERMINAL_SOC_SOFT=1 \
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/export_kkt_sparsity.jl <system> <T>

using Ipopt
using JuMP
using LinearAlgebra
using Printf
using Serialization
using SparseArrays

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const MOI = JuMP.MOI
include(joinpath(@__DIR__, "terminal_soc_penalty.jl"))

system = length(ARGS) >= 1 ? ARGS[1] : "ieee123C_1ph"
T = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
ptag = haskey(ENV, "REDUCED_PROFILE") ? "_" * ENV["REDUCED_PROFILE"] : ""
datafile = joinpath(REPO, "ddp", "results", "network_filterddp",
                    "network_data_$(system)_T$(T)$(ptag).jls")

data = deserialize(datafile)
haskey(ENV, "REDUCED_CB") && (data[:C_B] = parse(Float64, ENV["REDUCED_CB"]))

Nset, Lset, Bset, Dset = data[:Nset], data[:Lset], data[:Bset], data[:Dset]
Tset = 1:T
root, nonroot = data[:substationBus], data[:Nm1set]
dt, pbase = data[:delta_t_h], data[:kVA_B]

model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 0)

@variable(model, P_Subs[Tset] >= 0)
@variable(model, Q_Subs[Tset])
@variable(model, P[Lset, Tset])
@variable(model, Q[Lset, Tset])
@variable(model, v[Nset, Tset])
@variable(model, ell[Lset, Tset] >= 0)
@variable(model, P_B[Bset, Tset])
@variable(model, B[Bset, Tset])
@variable(model, q_D[Dset, Tset])

gammaT = terminal_soc_soft() ? gamma_terminal(system) : 0.0
@objective(model, Min,
    sum(data[:LoadShapeCost][t] * pbase * dt * P_Subs[t] for t in Tset) +
    sum(data[:C_B] * pbase^2 * dt * P_B[j,t]^2 for j in Bset, t in Tset) +
    gammaT * sum((B[j,T] - data[:B0_pu][j])^2 for j in Bset))

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

nvar = num_variables(model)
ncon = num_constraints(model; count_variable_in_set_constraints=false)
@printf("system=%s T=%d nvar=%d ncon=%d\n", system, T, nvar, ncon)

outdir = joinpath(REPO, "ddp", "results", "sparsity")
mkpath(outdir)
_varidx(v::VariableRef) = index(v).value

# Build the Jacobian as a sparse matrix by iterating over constraints
function extract_jacobian_sparsity(model)
    rows = Int[]
    cols = Int[]
    row = 0
    for (F, S) in list_of_constraint_types(model)
        F == VariableRef && continue
        for con in all_constraints(model, F, S)
            row += 1
            obj = constraint_object(con)
            func = obj.func
            _add_vars!(rows, cols, row, func, model)
        end
    end
    return rows, cols, row
end

function _add_vars!(rows, cols, row, func::AffExpr, model)
    for (var, _) in func.terms
        push!(rows, row)
        push!(cols, _varidx(var))
    end
end

function _add_vars!(rows, cols, row, func::QuadExpr, model)
    _add_vars!(rows, cols, row, func.aff, model)
    for (pair, _) in func.terms
        push!(rows, row)
        push!(cols, _varidx(pair.a))
        push!(rows, row)
        push!(cols, _varidx(pair.b))
    end
end

function _add_vars!(rows, cols, row, func, model)
    @warn "unhandled function type: $(typeof(func))"
end

@printf("extracting Jacobian sparsity... ")
flush(stdout)
jr, jc, m = extract_jacobian_sparsity(model)
J = sparse(jr, jc, ones(length(jr)), m, nvar)
dropzeros!(J)
@printf("done: %d x %d, nnz=%d\n", m, nvar, nnz(J))

# Hessian: only the objective contributes nonlinear terms (P_B^2 and SOC cone)
# For the objective, the Hessian is diagonal in P_B and B variables.
# For the SOC constraints (P^2 + Q^2 <= v*ell), the Hessian involves P, Q, v, ell.
# Extract by constraint type.
function extract_hessian_sparsity(model)
    rows = Int[]
    cols = Int[]
    obj = objective_function(model)
    if obj isa QuadExpr
        for (pair, _) in obj.terms
            i = _varidx(pair.a)
            j = _varidx(pair.b)
            push!(rows, i); push!(cols, j)
            if i != j
                push!(rows, j); push!(cols, i)
            end
        end
    end
    for (F, S) in list_of_constraint_types(model)
        F == VariableRef && continue
        F <: AffExpr && continue
        for con in all_constraints(model, F, S)
            obj_c = constraint_object(con)
            func = obj_c.func
            if func isa QuadExpr
                for (pair, _) in func.terms
                    i = _varidx(pair.a)
                    j = _varidx(pair.b)
                    push!(rows, i); push!(cols, j)
                    if i != j
                        push!(rows, j); push!(cols, i)
                    end
                end
            end
        end
    end
    return rows, cols
end

@printf("extracting Hessian sparsity... ")
flush(stdout)
hr, hc = extract_hessian_sparsity(model)
H = sparse(hr, hc, ones(length(hr)), nvar, nvar)
dropzeros!(H)
@printf("done: %d x %d, nnz=%d\n", nvar, nvar, nnz(H))

# Build KKT = [H J'; J 0]
n = nvar
KKT_dim = n + m
@printf("assembling KKT... ")
flush(stdout)
KKT = blockdiag(H, spzeros(m, m))
JI, JJ, _ = findnz(J)
for k in eachindex(JI)
    KKT[n + JI[k], JJ[k]] = 1.0
    KKT[JJ[k], n + JI[k]] = 1.0
end
dropzeros!(KKT)
@printf("done: %d x %d, nnz=%d\n", KKT_dim, KKT_dim, nnz(KKT))

# Write CSV: row, col for each nonzero (1-indexed)
function write_spy_csv(filename, S)
    open(filename, "w") do io
        println(io, "row,col")
        I, J, _ = findnz(S)
        for k in eachindex(I)
            println(io, I[k], ",", J[k])
        end
    end
end

prefix = "$(system)_T$(T)"
write_spy_csv(joinpath(outdir, "$(prefix)_jacobian.csv"), J)
write_spy_csv(joinpath(outdir, "$(prefix)_hessian.csv"), H)
write_spy_csv(joinpath(outdir, "$(prefix)_kkt.csv"), KKT)

# Also write dimensions
open(joinpath(outdir, "$(prefix)_dims.txt"), "w") do io
    @printf(io, "system=%s T=%d nvar=%d ncon=%d nnz_J=%d nnz_H=%d nnz_KKT=%d\n",
            system, T, nvar, m, nnz(J), nnz(H), nnz(KKT))
end

@printf("wrote %s/{jacobian,hessian,kkt}.csv + dims.txt\n", prefix)
