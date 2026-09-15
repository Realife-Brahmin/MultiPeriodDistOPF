# Persistent inner OPF: build the JuMP model ONCE per stage, then re-solve it
# at new battery dispatches by editing the balance right-hand sides in place.
#
# Motivation, measured rather than assumed. Warm-starting alone cut the inner
# solve from ~35 interior-point iterations to ~4-5, but wall time fell only
# 3.8-4.7x against a 7-9.4x iteration drop -- the residue is JuMP model
# CONSTRUCTION, which `inner_opf` repeats on every single call. P_B enters the
# real-power balance rows purely as a constant, so it can be moved with
# `set_normalized_rhs` and the model reused.
#
# `inner_opf` in inner_network_opf.jl is deliberately left untouched: it is the
# validated reference (checked against FilterDDP's own constraint callback to
# 6.7e-11), and this file is checked AGAINST it rather than assumed equivalent.
# Run this file directly to see that cross-check.

using JuMP
using Ipopt
using Printf

mutable struct PersistentInner
    model::Model
    data::Dict
    t::Int
    g::NamedTuple
    ps::VariableRef
    v::Vector{VariableRef}
    ell::Vector{VariableRef}
    P::Vector{VariableRef}
    Q::Vector{VariableRef}
    qnorm::Vector{VariableRef}
    bal_p::Dict{Int,ConstraintRef}
    rhs_const::Dict{Int,Float64}     # pD - pL contribution, fixed for this t
    allvars::Vector{VariableRef}
    allcons::Vector{Any}
    solved_once::Bool
end

"""
    persistent_inner(data, t; tol, max_iter)

Build the single-period network OPF once, with `P_B = 0`. Later dispatches are
applied by `solve_at!`. Hard voltage limits only -- the soft/penalty path stays
in `inner_opf`, where the coefficient logic lives.
"""
function persistent_inner(data::Dict, t::Int; tol::Float64=1e-10,
                          max_iter::Int=3000, silent::Bool=true)
    g = network_slice(data, t)
    nL = length(g.lines); nN = length(g.buses); nD = length(g.ders)

    model = Model(Ipopt.Optimizer)
    silent && set_silent(model)
    set_optimizer_attribute(model, "tol", tol)
    set_optimizer_attribute(model, "constr_viol_tol", tol)
    set_optimizer_attribute(model, "acceptable_tol", tol * 1e2)
    set_optimizer_attribute(model, "max_iter", max_iter)

    @variable(model, ps >= 0.0)
    @variable(model, qs)
    @variable(model, P[1:nL])
    @variable(model, Q[1:nL])
    @variable(model, v[1:nN])
    @variable(model, ell[1:nL] >= 0.0)
    @variable(model, -1.0 <= qnorm[1:nD] <= 1.0)
    @variable(model, soc_slack[1:nL] >= 0.0)
    for (k, j) in enumerate(g.buses)
        set_lower_bound(v[k], data[:Vminpu][j]^2)
        set_upper_bound(v[k], data[:Vmaxpu][j]^2)
    end
    set_start_value(ps, 1.0); set_start_value(qs, 0.3)
    for k in 1:nN; set_start_value(v[k], 1.0); end
    for e in 1:nL; set_start_value(ell[e], 1e-3); set_start_value(soc_slack[e], 1e-3); end

    @constraint(model, ps - sum(P[g.linepos[e]] for e in data[:L1set]) == 0.0)
    bal_p = Dict{Int,ConstraintRef}()
    rhs_const = Dict{Int,Float64}()
    for j in g.nonroot
        inc = (data[:parent][j], j)
        pL = j in data[:NLset] ? data[:p_L_pu][j, t] : 0.0
        pD = j in g.ders ? data[:p_D_pu][j, t] : 0.0
        r = data[:rdict_pu][inc]
        # built at P_B = 0; the constant that will carry P_B is the RHS
        bal_p[j] = @constraint(model,
            sum(P[g.linepos[(j, k)]] for k in data[:children][j]; init = 0.0)
            - P[g.linepos[inc]] + r * ell[g.linepos[inc]] - pD + pL == 0.0)
        rhs_const[j] = pD - pL
    end
    @constraint(model, qs - sum(Q[g.linepos[e]] for e in data[:L1set]) == 0.0)
    for j in g.nonroot
        inc = (data[:parent][j], j)
        qL = j in data[:NLset] ? data[:q_L_pu][j, t] : 0.0
        z = data[:xdict_pu][inc]
        qD = j in g.ders ? g.qmax[j] * qnorm[g.derpos[j]] : 0.0
        @constraint(model,
            sum(Q[g.linepos[(j, k)]] for k in data[:children][j]; init = 0.0)
            - Q[g.linepos[inc]] + z * ell[g.linepos[inc]] - qD + qL == 0.0)
    end
    for (e, (i, j)) in enumerate(g.lines)
        r, z = data[:rdict_pu][(i, j)], data[:xdict_pu][(i, j)]
        @constraint(model, v[g.buspos[j]] - v[g.buspos[i]] +
                           2 * (r * P[e] + z * Q[e]) - (r^2 + z^2) * ell[e] == 0.0)
    end
    for (e, (i, _)) in enumerate(g.lines)
        @constraint(model, P[e]^2 + Q[e]^2 - v[g.buspos[i]] * ell[e] + soc_slack[e] == 0.0)
    end
    @constraint(model, v[g.buspos[g.root]] - 1.05^2 == 0.0)

    price = data[:LoadShapeCost][t]
    @objective(model, Min, price * data[:kVA_B] * data[:delta_t_h] * ps)

    return PersistentInner(model, data, t, g, ps, collect(v), collect(ell),
                           collect(P), collect(Q), collect(qnorm), bal_p, rhs_const,
                           all_variables(model),
                           Any[c for c in all_constraints(model;
                               include_variable_in_set_constraints = true)], false)
end

"""
    solve_at!(pi, pb; warm=true)

Re-solve at dispatch `pb`. Only the balance right-hand sides change; the model,
its sparsity and its Ipopt instance are reused. With `warm=true` the previous
solution seeds the next solve.
"""
function solve_at!(pin::PersistentInner, pb::Vector{Float64}; warm::Bool=true)
    g = pin.g
    # Capture the previous solution BEFORE touching the model: editing any
    # right-hand side invalidates it, and querying afterwards throws
    # OptimizeNotCalled. Order matters here, not style.
    # Primals alone are not enough: restoring the BOUND multipliers is what
    # stops Ipopt re-walking the central path. With primals only this workload
    # took 22-28 interior-point iterations; with the duals it takes ~4.
    prev  = (warm && pin.solved_once) ? [value(vr) for vr in pin.allvars] : nothing
    prevy = prev === nothing ? nothing :
            [(try dual(cr) catch; 0.0 end) for cr in pin.allcons]
    for (b, j) in enumerate(g.batteries)
        haskey(pin.bal_p, j) || continue      # root-bus battery has no balance row
        set_normalized_rhs(pin.bal_p[j], pb[b] + pin.rhs_const[j])
    end
    if prev !== nothing
        set_optimizer_attribute(pin.model, "warm_start_init_point", "yes")
        set_optimizer_attribute(pin.model, "warm_start_bound_push", 1e-9)
        set_optimizer_attribute(pin.model, "warm_start_bound_frac", 1e-9)
        set_optimizer_attribute(pin.model, "warm_start_slack_bound_push", 1e-9)
        set_optimizer_attribute(pin.model, "warm_start_slack_bound_frac", 1e-9)
        set_optimizer_attribute(pin.model, "warm_start_mult_bound_push", 1e-9)
        set_optimizer_attribute(pin.model, "mu_init", 1e-7)
        for (vr, val) in zip(pin.allvars, prev)
            set_start_value(vr, val)
        end
        for (cr, val) in zip(pin.allcons, prevy)
            try; set_dual_start_value(cr, val); catch; end
        end
    end
    t0 = time()
    optimize!(pin.model)
    st = time() - t0
    status = termination_status(pin.model)
    ok = status in (MOI_.LOCALLY_SOLVED, MOI_.OPTIMAL, MOI_.ALMOST_LOCALLY_SOLVED)
    iters = try MOI_.get(pin.model, MOI_.BarrierIterations()) catch; -1 end
    ok && (pin.solved_once = true)
    if !ok
        return (; feasible = false, status = string(status), iterations = iters,
                solve_time = st, substation_cost = NaN, lambda_bal = Float64[],
                ps = NaN, vmin = NaN)
    end
    lam = Float64[haskey(pin.bal_p, j) ? dual(pin.bal_p[j]) : 0.0 for j in g.batteries]
    psv = value(pin.ps)
    price = pin.data[:LoadShapeCost][pin.t]
    return (; feasible = true, status = string(status), iterations = iters,
            solve_time = st,
            substation_cost = price * pin.data[:kVA_B] * pin.data[:delta_t_h] * psv,
            lambda_bal = lam, ps = psv,
            vmin = sqrt(max(minimum(value.(pin.v)), 0.0)))
end

# ------------------------------------------------------------- self-check ---
function _selfcheck(sys::String)
    data = deserialize(joinpath(normpath(joinpath(@__DIR__, "..", "..", "..")),
                                "ddp", "results", "network_filterddp",
                                "network_data_$(sys)_T3.jls"))
    nB = length(data[:Bset])
    pbr = Float64[data[:P_B_R_pu][j] for j in data[:Bset]]
    pin = persistent_inner(data, 1)
    worst_c = 0.0; worst_g = 0.0
    tp = 0.0; tr = 0.0; itp = 0; itr = 0
    for s in 1:6
        rng = MersenneTwister(4000 + s)
        pb = s == 1 ? zeros(nB) : [(2rand(rng) - 1) * 0.5 * pbr[b] for b in 1:nB]
        a = solve_at!(pin, pb)
        b = inner_opf(data, 1, pb)
        tp += a.solve_time; tr += b.solve_time; itp += a.iterations; itr += b.iterations
        a.feasible == b.feasible || error("feasibility disagreement at sample $s")
        if a.feasible
            worst_c = max(worst_c, abs(a.substation_cost - b.substation_cost))
            worst_g = max(worst_g, maximum(abs, a.lambda_bal .- b.lambda_bal))
        end
    end
    @printf("%-16s max |Phi| diff vs inner_opf = %.3e USD   max |grad| diff = %.3e
",
            sys, worst_c, worst_g)
    @printf("                 persistent %.3f s / %3d iters    reference %.3f s / %3d iters    SPEEDUP %.1fx
",
            tp, itp, tr, itr, tr / max(tp, 1e-9))
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    using LinearAlgebra, Random, Serialization
    include(joinpath(@__DIR__, "inner_network_opf.jl"))
    for sys in (isempty(ARGS) ? ["ieee123C_1ph"] : ARGS)
        _selfcheck(sys)
    end
end
