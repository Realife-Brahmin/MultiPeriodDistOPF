# Starting point from a power flow (FILTERDDP_LOADFLOW_START=1 in the FilterDDP
# driver, IPOPT_LOADFLOW_START=1 in centralized_ipopt_matched.jl).
#
# The default start is flat: every voltage 1.0, every flow 0, every ell 1e-3.
# From there the solver first has to find the power flow itself, through the
# barrier and the filter (on large10k: about 70 of 100 iterations at steps of
# 0.03 or less; the substation import alone has to grow from 1 to several
# hundred p.u.). That part needs no optimisation: with the batteries idle and
# the DERs at zero reactive power, the network state of each period is one
# backward/forward sweep on the tree.
#
# The sweep solves the branch-flow equations of the driver exactly:
#   backward  P_recv = subtree load, ell = (P_recv^2 + Q_recv^2) / v_j,
#             P = P_recv + r ell,  Q = Q_recv + x ell
#   forward   v_j = v_i - 2 (r P + x Q) + (r^2 + x^2) ell
# and is repeated until the voltages stop changing.

# Network state of period t with idle batteries and zero DER reactive power:
# v in data[:Nset] order; P, Q, ell in data[:Lset] order; and the sweep count.
function loadflow_sweep(data, t::Int; sweeps::Int=50)
    buses, lines = data[:Nset], data[:Lset]
    root, children, parent = data[:substationBus], data[:children], data[:parent]
    buspos = Dict(j => k for (k, j) in enumerate(buses))
    linepos = Dict(e => k for (k, e) in enumerate(lines))
    order = Int[]                                   # parents before children
    stack = [root]
    while !isempty(stack)
        k = pop!(stack); push!(order, k)
        append!(stack, children[k])
    end
    NL, D = Set(data[:NLset]), Set(data[:Dset])
    v = fill(1.05^2, length(buses))
    P = zeros(length(lines)); Q = zeros(length(lines)); ell = zeros(length(lines))
    change = Inf; n = 0
    while n < sweeps && change > 1e-13
        n += 1
        for k in Iterators.reverse(order)
            k == root && continue
            e = (parent[k], k); le = linepos[e]
            pr = (k in NL ? data[:p_L_pu][k,t] : 0.0) - (k in D ? data[:p_D_pu][k,t] : 0.0)
            qr = (k in NL ? data[:q_L_pu][k,t] : 0.0)
            for c in children[k]
                lc = linepos[(k, c)]
                pr += P[lc]; qr += Q[lc]
            end
            ell[le] = (pr^2 + qr^2) / v[buspos[k]]
            P[le] = pr + data[:rdict_pu][e] * ell[le]
            Q[le] = qr + data[:xdict_pu][e] * ell[le]
        end
        change = 0.0
        for k in order
            k == root && continue
            e = (parent[k], k); le = linepos[e]
            r, x = data[:rdict_pu][e], data[:xdict_pu][e]
            vn = v[buspos[parent[k]]] - 2(r * P[le] + x * Q[le]) + (r^2 + x^2) * ell[le]
            change = max(change, abs(vn - v[buspos[k]]))
            v[buspos[k]] = vn
        end
    end
    change <= 1e-8 || error("power-flow sweep did not converge at period $t (change $change)")
    return v, P, Q, ell, n
end

# Fill the FilterDDP control vector of period t. `soc_slack` is the solver's
# own interior push for a variable bounded below by zero (kappa_1 = 0.01), and
# ell is raised to keep the SOC row satisfied with that slack.
function loadflow_start!(u::AbstractVector, data, idx, t::Int; soc_slack::Float64=0.01, sweeps::Int=50)
    v, P, Q, _, n = loadflow_sweep(data, t; sweeps=sweeps)
    lines = data[:Lset]
    buspos = Dict(j => k for (k, j) in enumerate(data[:Nset]))
    linepos = Dict(e => k for (k, e) in enumerate(lines))
    for (le, e) in enumerate(lines)
        u[idx.P[le]] = P[le]; u[idx.Q[le]] = Q[le]
        u[idx.soc_slack[le]] = soc_slack
        u[idx.ell[le]] = (P[le]^2 + Q[le]^2 + soc_slack) / v[buspos[e[1]]]
    end
    u[idx.v] .= v
    u[idx.ps] = sum(P[linepos[e]] for e in data[:L1set])
    u[idx.qs] = sum(Q[linepos[e]] for e in data[:L1set])
    return n
end
