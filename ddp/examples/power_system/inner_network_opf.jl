# Inner single-period network OPF with the battery dispatch held FIXED.
#
# Reduced-space MPOPF probe (opt-in diagnostic; nothing here is imported by the
# production FilterDDP path). The proposed decomposition is:
#
#   outer  : B^{t-1}, P_B^t, battery dynamics, battery bounds
#   inner  : given a fixed P_B^t, solve every algebraic network quantity
#            (P_Subs, Q_Subs, branch P/Q, v, ell, DER reactive, SOCP slacks)
#
# This file provides the inner solve and a separate feasibility-restoration
# diagnostic. The equations, profiles, bounds, objective terms and units are
# mirrored from the centralized transcription in `ieee123c_filterddp.jl`
# (`build_model`/`control_layout`); every solve is additionally checked by
# evaluating FilterDDP's own stage-constraint callback on the recovered control
# vector, which is an independent code path from the JuMP model built here.
#
# Modelling decision, stated rather than hidden: the primary inner problem is
# NETWORK-only. The energy-slack row `B^{t-1} - dt*P_B^t - Bmin - s_E = 0` with
# `s_E in [0, Bmax-Bmin]` is an outer-layer battery bound -- with both B^{t-1}
# and P_B^t fixed it contains no optimisation freedom at all, it just evaluates
# whether the dispatch respects the state-of-charge box. Including it would
# conflate "this dispatch breaks the network" with "this dispatch breaks the
# battery", which is exactly the distinction this study needs. It is available
# via `include_energy_row=true` and reported separately.

using JuMP
using Ipopt
using LinearAlgebra
using Printf
using Serialization
using SparseArrays

const MOI_ = JuMP.MOI

"""
Geometry/profile slice that both the inner OPF and the restoration model need.
"""
function network_slice(data::Dict, t::Int)
    buses, lines = data[:Nset], data[:Lset]
    batteries, ders = data[:Bset], data[:Dset]
    root = data[:substationBus]
    nonroot = data[:Nm1set]
    buspos = Dict(j => k for (k, j) in enumerate(buses))
    linepos = Dict(e => k for (k, e) in enumerate(lines))
    batpos = Dict(j => k for (k, j) in enumerate(batteries))
    derpos = Dict(j => k for (k, j) in enumerate(ders))
    qmax = Dict(j => sqrt(max(0.0, data[:S_D_R][j]^2 - data[:p_D_pu][j, t]^2)) for j in ders)
    return (; buses, lines, batteries, ders, root, nonroot,
            buspos, linepos, batpos, derpos, qmax)
end

"""
    inner_opf(data, t, pb_fixed; kwargs...)

Solve the single-period network OPF with `P_B^t = pb_fixed` (per-unit, indexed
like `data[:Bset]`, positive = discharge/injection). Returns a NamedTuple of the
solution and every diagnostic the survey records.
"""
function inner_opf(data::Dict, t::Int, pb_fixed::Vector{Float64};
                   include_energy_row::Bool=false,
                   x_entering::Union{Nothing,Vector{Float64}}=nothing,
                   voltage_penalty::Float64=0.0,
                   tol::Float64=1e-10, max_iter::Int=3000, silent::Bool=true)
    g = network_slice(data, t)
    dt = data[:delta_t_h]
    pbase = data[:kVA_B]
    price = data[:LoadShapeCost][t]
    nB = length(g.batteries)
    length(pb_fixed) == nB || error("pb_fixed has $(length(pb_fixed)) entries, expected $nB")
    include_energy_row && x_entering === nothing &&
        error("include_energy_row=true requires x_entering")

    model = Model(Ipopt.Optimizer)
    silent && set_silent(model)
    set_optimizer_attribute(model, "tol", tol)
    set_optimizer_attribute(model, "constr_viol_tol", tol)
    set_optimizer_attribute(model, "acceptable_tol", tol * 1e2)
    set_optimizer_attribute(model, "max_iter", max_iter)

    nL = length(g.lines)
    nN = length(g.buses)
    nD = length(g.ders)

    # ---- variables (bounds identical to build_model's ControlLimits) --------
    @variable(model, ps >= 0.0)                       # P_Subs: lower bound only
    @variable(model, qs)                              # Q_Subs: free
    @variable(model, P[1:nL])
    @variable(model, Q[1:nL])
    @variable(model, v[k = 1:nN])
    @variable(model, ell[1:nL] >= 0.0)
    @variable(model, -1.0 <= qnorm[1:nD] <= 1.0)
    @variable(model, soc_slack[1:nL] >= 0.0)

    # Voltage limits: hard by default. With voltage_penalty > 0 they become
    # soft, enforced by an L1 exact penalty on nonnegative violation slacks.
    # "Situational" in the sense that matters: the penalty term is identically
    # zero at any dispatch the network can actually serve, so wherever the hard
    # problem is feasible the penalised problem has the same solution (provided
    # the coefficient exceeds the largest voltage-constraint dual). Where the
    # hard problem is infeasible, Phi_t is still defined and the penalty
    # measures how far outside F_t the dispatch sits.
    soft_voltage = voltage_penalty > 0.0
    s_vlo = nothing; s_vhi = nothing
    if soft_voltage
        @variable(model, svlo[1:nN] >= 0.0)
        @variable(model, svhi[1:nN] >= 0.0)
        for (k, j) in enumerate(g.buses)
            @constraint(model, v[k] >= data[:Vminpu][j]^2 - svlo[k])
            @constraint(model, v[k] <= data[:Vmaxpu][j]^2 + svhi[k])
        end
        s_vlo = svlo; s_vhi = svhi
    else
        for (k, j) in enumerate(g.buses)
            set_lower_bound(v[k], data[:Vminpu][j]^2)
            set_upper_bound(v[k], data[:Vmaxpu][j]^2)
        end
    end
    # warm start near flat voltage keeps Ipopt off the v*ell degeneracy
    set_start_value(ps, 1.0)
    set_start_value(qs, 0.3)
    for k in 1:nN; set_start_value(v[k], 1.0); end
    for e in 1:nL; set_start_value(ell[e], 1e-3); set_start_value(soc_slack[e], 1e-3); end

    # ---- equality constraints (mirrors `equations(x,u)` exactly) ------------
    @constraint(model, bal_p_root,
        ps - sum(P[g.linepos[e]] for e in data[:L1set]) == 0.0)
    bal_p = Dict{Int,ConstraintRef}()
    for j in g.nonroot
        inc_line = (data[:parent][j], j)
        pL = j in data[:NLset] ? data[:p_L_pu][j, t] : 0.0
        pD = j in g.ders ? data[:p_D_pu][j, t] : 0.0
        pb = j in g.batteries ? pb_fixed[g.batpos[j]] : 0.0
        r = data[:rdict_pu][inc_line]
        bal_p[j] = @constraint(model,
            sum(P[g.linepos[(j, k)]] for k in data[:children][j]; init = 0.0)
            - P[g.linepos[inc_line]] + r * ell[g.linepos[inc_line]]
            - pb - pD + pL == 0.0)
    end
    @constraint(model, bal_q_root,
        qs - sum(Q[g.linepos[e]] for e in data[:L1set]) == 0.0)
    bal_q = Dict{Int,ConstraintRef}()
    for j in g.nonroot
        inc_line = (data[:parent][j], j)
        qL = j in data[:NLset] ? data[:q_L_pu][j, t] : 0.0
        z = data[:xdict_pu][inc_line]
        qD = j in g.ders ? g.qmax[j] * qnorm[g.derpos[j]] : 0.0
        bal_q[j] = @constraint(model,
            sum(Q[g.linepos[(j, k)]] for k in data[:children][j]; init = 0.0)
            - Q[g.linepos[inc_line]] + z * ell[g.linepos[inc_line]]
            - qD + qL == 0.0)
    end
    drop = Vector{ConstraintRef}(undef, nL)
    for (e, (i, j)) in enumerate(g.lines)
        r, z = data[:rdict_pu][(i, j)], data[:xdict_pu][(i, j)]
        drop[e] = @constraint(model,
            v[g.buspos[j]] - v[g.buspos[i]] + 2 * (r * P[e] + z * Q[e])
            - (r^2 + z^2) * ell[e] == 0.0)
    end
    cone = Vector{ConstraintRef}(undef, nL)
    for (e, (i, _)) in enumerate(g.lines)
        cone[e] = @constraint(model,
            P[e]^2 + Q[e]^2 - v[g.buspos[i]] * ell[e] + soc_slack[e] == 0.0)
    end
    @constraint(model, vsub, v[g.buspos[g.root]] - 1.05^2 == 0.0)

    energy_slack = nothing
    if include_energy_row
        width = [(data[:soc_max][j] - data[:soc_min][j]) * data[:B_R_pu][j] for j in g.batteries]
        @variable(model, es[b = 1:nB])
        for b in 1:nB
            set_lower_bound(es[b], 0.0); set_upper_bound(es[b], width[b])
        end
        for (b, j) in enumerate(g.batteries)
            emin = data[:soc_min][j] * data[:B_R_pu][j]
            @constraint(model, x_entering[b] - dt * pb_fixed[b] - emin - es[b] == 0.0)
        end
        energy_slack = es
    end

    # ---- objective: identical time-local terms to analytic_objective -------
    # l = c^t * S_base * dt * P_Subs + C_B * S_base^2 * dt * sum(P_B^2)
    # with P_B fixed the second term is a constant; it is reported, not optimised.
    battery_term = data[:C_B] * pbase^2 * dt * sum(abs2, pb_fixed)
    if soft_voltage
        @objective(model, Min, price * pbase * dt * ps +
                               voltage_penalty * (sum(s_vlo) + sum(s_vhi)))
    else
        @objective(model, Min, price * pbase * dt * ps)
    end

    t_start = time()
    rss_start = Sys.maxrss()
    gc_start = Base.gc_bytes()
    optimize!(model)
    solve_time = time() - t_start
    alloc_mib = (Base.gc_bytes() - gc_start) / 2^20
    rss_delta_mib = max(Sys.maxrss() - rss_start, 0) / 2^20

    status = termination_status(model)
    ok = status in (MOI_.LOCALLY_SOLVED, MOI_.OPTIMAL, MOI_.ALMOST_LOCALLY_SOLVED)
    iters = try MOI_.get(model, MOI_.BarrierIterations()) catch; -1 end

    if !ok
        return (; feasible = false, status = string(status), iterations = iters,
                solve_time, alloc_mib, rss_delta_mib,
                objective = NaN, substation_cost = NaN, battery_term,
                soft_voltage, voltage_penalty, penalty_cost = NaN,
                total_violation = NaN, max_violation = NaN, n_violated = -1,
                ps = NaN, qs = NaN, vmin = NaN, vmax = NaN,
                ell_max = NaN, imag_max = NaN, qnorm_absmax = NaN, qnorm_absmean = NaN,
                min_bound_margin = NaN, ps_margin = NaN, reverse_export = false,
                filterddp_residual = NaN, validation_pass = false,
                loss = NaN, net_load = NaN, energy_slack_ok = missing)
    end

    psv = value(ps); qsv = value(qs)
    vv = value.(v); ellv = value.(ell); qn = value.(qnorm)
    Pv = value.(P); Qv = value.(Q); ssv = value.(soc_slack)
    vmag = sqrt.(max.(vv, 0.0))

    # Duals of the real-power balance at each battery bus. With JuMP's sign
    # convention for this equality row, this dual IS d(Phi_t)/d(P_B) -- the
    # analytic gradient of the inner optimal value with respect to that
    # battery's dispatch, available for free from the solve. Verified against
    # central differences to ~1e-9 relative in probe_reduced_value_function.jl.
    # A battery sitting on the substation bus has NO balance row: the root row
    # is `ps - sum(P_out)` and carries no pb term, so such a battery is inert in
    # this transcription and its true derivative is exactly zero, not missing.
    # (ieee2522C_1ph has one such battery, bus 1, rated 0.00426 pu. See
    # ddp/notes/REDUCED_SPACE_INNER_OPF_FEASIBILITY.md.)
    lambda_bal = Float64[]
    for j in g.batteries
        push!(lambda_bal, haskey(bal_p, j) ? dual(bal_p[j]) : 0.0)
    end

    # bound margins over every bounded variable
    margins = Float64[psv - 0.0]
    for k in 1:nN
        push!(margins, vv[k] - data[:Vminpu][g.buses[k]]^2)
        push!(margins, data[:Vmaxpu][g.buses[k]]^2 - vv[k])
    end
    append!(margins, ellv)
    append!(margins, ssv)
    for dd in 1:nD; push!(margins, qn[dd] + 1.0); push!(margins, 1.0 - qn[dd]); end

    net_load = sum(data[:p_L_pu][j, t] for j in data[:NLset]) -
               sum(data[:p_D_pu][j, t] for j in g.ders)
    loss = psv + sum(pb_fixed) - net_load   # real loss implied by the balance

    energy_ok = missing
    if include_energy_row
        esv = value.(energy_slack)
        width = [(data[:soc_max][j] - data[:soc_min][j]) * data[:B_R_pu][j] for j in g.batteries]
        energy_ok = all(-1e-8 .<= esv .<= width .+ 1e-8)
    end

    # Split the penalty back out so Phi_t stays comparable across hard and soft
    # runs: substation_cost is always the priced substation energy alone.
    vio_lo = soft_voltage ? value.(s_vlo) : zeros(nN)
    vio_hi = soft_voltage ? value.(s_vhi) : zeros(nN)
    total_violation = sum(vio_lo) + sum(vio_hi)
    penalty_cost = voltage_penalty * total_violation
    substation_cost = price * pbase * dt * psv

    return (; feasible = true, status = string(status), iterations = iters,
            solve_time, alloc_mib, rss_delta_mib,
            objective = substation_cost + battery_term,
            substation_cost, battery_term,
            soft_voltage, voltage_penalty, penalty_cost,
            total_violation, max_violation = max(maximum(vio_lo; init = 0.0),
                                                 maximum(vio_hi; init = 0.0)),
            n_violated = count(>(1e-9), vio_lo) + count(>(1e-9), vio_hi),
            ps = psv, qs = qsv,
            vmin = minimum(vmag), vmax = maximum(vmag),
            ell_max = maximum(ellv), imag_max = sqrt(maximum(ellv)),
            qnorm_absmax = maximum(abs, qn), qnorm_absmean = sum(abs, qn) / nD,
            min_bound_margin = minimum(margins),
            ps_margin = psv,
            reverse_export = psv <= 1e-7,
            loss, net_load,
            energy_slack_ok = energy_ok,
            lambda_bal,
            n_active_qnorm = count(x -> abs(abs(x) - 1.0) < 1e-6, qn),
            n_active_vmin = count(k -> vv[k] - data[:Vminpu][g.buses[k]]^2 < 1e-8, 1:nN),
            n_active_vmax = count(k -> data[:Vmaxpu][g.buses[k]]^2 - vv[k] < 1e-8, 1:nN),
            Pv, Qv, vv, ellv, qn, ssv)
end

"""
    filterddp_residual(ocp, idx, data, t, res, pb_fixed, x_entering)

Independent validation: rebuild the full FilterDDP control vector from the JuMP
solution and evaluate FilterDDP's own analytic stage-constraint callback.
Different code path, same physics -- a small residual means the JuMP model in
this file really is the centralized model.
"""
function filterddp_residual(ocp, idx, data::Dict, t::Int, res, pb_fixed::Vector{Float64},
                            x_entering::Vector{Float64})
    g = network_slice(data, t)
    nu = last(idx.energy_slack)      # layout is contiguous; last slot ends the vector
    u = zeros(nu)
    u[idx.ps] = res.ps
    u[idx.qs] = res.qs
    u[idx.P] .= res.Pv
    u[idx.Q] .= res.Qv
    u[idx.v] .= res.vv
    u[idx.ell] .= res.ellv
    u[idx.pb] .= pb_fixed
    u[idx.qnorm] .= res.qn
    u[idx.soc_slack] .= res.ssv
    for (b, j) in enumerate(g.batteries)
        emin = data[:soc_min][j] * data[:B_R_pu][j]
        u[idx.energy_slack[b]] = x_entering[b] - data[:delta_t_h] * pb_fixed[b] - emin
    end
    c = ocp.stage_constraints[t].c(x_entering, u)
    # the energy row is exactly satisfied by construction above; the network
    # rows are the real test
    return maximum(abs, c)
end

"""
    restore_feasibility(data, t, pb_fixed)

DIAGNOSTIC ONLY -- never a feasibility claim for the original problem. Adds
nonnegative violation slacks to each physical constraint class and minimises a
dominantly-weighted total violation, so an infeasible dispatch reports *which*
physical limit it breaks and by how much.
"""
function restore_feasibility(data::Dict, t::Int, pb_fixed::Vector{Float64};
                             tol::Float64=1e-9, max_iter::Int=3000)
    g = network_slice(data, t)
    pbase = data[:kVA_B]; dt = data[:delta_t_h]; price = data[:LoadShapeCost][t]
    nL = length(g.lines); nN = length(g.buses); nD = length(g.ders)

    model = Model(Ipopt.Optimizer); set_silent(model)
    set_optimizer_attribute(model, "tol", tol)
    set_optimizer_attribute(model, "max_iter", max_iter)

    @variable(model, ps)                       # relaxed: sign free here
    @variable(model, qs)
    @variable(model, P[1:nL]); @variable(model, Q[1:nL])
    @variable(model, v[1:nN]); @variable(model, ell[1:nL] >= 0.0)
    @variable(model, -1.0 <= qnorm[1:nD] <= 1.0)
    @variable(model, soc_slack[1:nL] >= 0.0)
    # violation slacks, one per physical class
    @variable(model, s_export >= 0.0)              # P_Subs >= 0 violation
    @variable(model, s_vlo[1:nN] >= 0.0)           # undervoltage
    @variable(model, s_vhi[1:nN] >= 0.0)           # overvoltage
    for k in 1:nN
        j = g.buses[k]
        @constraint(model, v[k] >= data[:Vminpu][j]^2 - s_vlo[k])
        @constraint(model, v[k] <= data[:Vmaxpu][j]^2 + s_vhi[k])
    end
    @constraint(model, ps >= -s_export)

    @constraint(model, ps - sum(P[g.linepos[e]] for e in data[:L1set]) == 0.0)
    for j in g.nonroot
        inc = (data[:parent][j], j)
        pL = j in data[:NLset] ? data[:p_L_pu][j, t] : 0.0
        pD = j in g.ders ? data[:p_D_pu][j, t] : 0.0
        pb = j in g.batteries ? pb_fixed[g.batpos[j]] : 0.0
        @constraint(model,
            sum(P[g.linepos[(j, k)]] for k in data[:children][j]; init = 0.0)
            - P[g.linepos[inc]] + data[:rdict_pu][inc] * ell[g.linepos[inc]]
            - pb - pD + pL == 0.0)
    end
    @constraint(model, qs - sum(Q[g.linepos[e]] for e in data[:L1set]) == 0.0)
    for j in g.nonroot
        inc = (data[:parent][j], j)
        qL = j in data[:NLset] ? data[:q_L_pu][j, t] : 0.0
        qD = j in g.ders ? g.qmax[j] * qnorm[g.derpos[j]] : 0.0
        @constraint(model,
            sum(Q[g.linepos[(j, k)]] for k in data[:children][j]; init = 0.0)
            - Q[g.linepos[inc]] + data[:xdict_pu][inc] * ell[g.linepos[inc]]
            - qD + qL == 0.0)
    end
    for (e, (i, j)) in enumerate(g.lines)
        r, z = data[:rdict_pu][(i, j)], data[:xdict_pu][(i, j)]
        @constraint(model, v[g.buspos[j]] - v[g.buspos[i]] + 2 * (r * P[e] + z * Q[e])
                           - (r^2 + z^2) * ell[e] == 0.0)
    end
    for (e, (i, _)) in enumerate(g.lines)
        @constraint(model, P[e]^2 + Q[e]^2 - v[g.buspos[i]] * ell[e] + soc_slack[e] == 0.0)
    end
    @constraint(model, v[g.buspos[g.root]] - 1.05^2 == 0.0)

    for k in 1:nN; set_start_value(v[k], 1.0); end
    for e in 1:nL; set_start_value(ell[e], 1e-3); set_start_value(soc_slack[e], 1e-3); end

    BIG = 1e6   # dominant weight: violation is minimised before any cost
    @objective(model, Min,
        BIG * (s_export + sum(s_vlo) + sum(s_vhi)) + price * pbase * dt * ps)
    optimize!(model)
    st = termination_status(model)
    ok = st in (MOI_.LOCALLY_SOLVED, MOI_.OPTIMAL, MOI_.ALMOST_LOCALLY_SOLVED)
    ok || return (; restored = false, status = string(st))
    return (; restored = true, status = string(st),
            export_violation = value(s_export),
            undervoltage_max = maximum(value.(s_vlo)),
            overvoltage_max = maximum(value.(s_vhi)),
            total_violation = value(s_export) + sum(value.(s_vlo)) + sum(value.(s_vhi)),
            ps = value(ps))
end
