# Voltage limits that cannot bind, identified from the network data alone (no
# OPF solution). Constraint screening; see ddp/notes/CONSTRAINT_SCREENING.md.
#
# In receiving-end form the voltage drop over line (i, j) is
#     v_i - v_j = 2 (r P_recv + x Q_recv) + |z|^2 ell,
# and the flow arriving at j is the net load of j's subtree plus the losses
# inside it. With ell >= 0 and r, x >= 0 this gives, for any feasible point,
#     P_recv >= Pmin_j,  Q_recv >= Qmin_j,
# the lossless subtree loads with every battery discharging at rating and
# every DER at its reactive limit. Hence
#     v_i - v_j >= d_j := 2 (r Pmin_j + x Qmin_j).                       (*)
#
#   substation  v_root is fixed by its own equality row, so its limits are
#               redundant. Its value 1.05^2 EQUALS the upper limit: the
#               variable has no interior and every full Newton step lands on
#               the bound.
#   upper bound summing (*) from the substation, v_j <= vub_j := 1.05^2 -
#               sum_path d. If vub_j is within the upper limit, that limit
#               cannot bind.
#   downhill    if d_j >= 0, voltage cannot rise from i to j. Then j's upper
#               limit is implied by i's, and i's lower limit by j's. (A bus
#               with no battery or DER below it always has d_j >= 0.)
#
#   P_Subs      the import is at least the sum of Pmin over the substation's
#               lines; if that is positive, P_Subs >= 0 cannot bind.
#
# All of these are exact for the relaxed branch-flow model: a limit dropped
# here is satisfied by every point that meets the remaining constraints. (So
# is ell >= 0, which the SOC row P^2 + Q^2 + s = v ell implies whenever s >= 0
# and v > 0; it needs no data and is handled in the driver.)

struct VoltageScreen
    buses::Vector{Int}          # data[:Nset] order
    keep_lo::BitMatrix          # [bus position, period]: lower limit still needed
    keep_up::BitMatrix
    vub::Matrix{Float64}        # lossless upper bound on v (squared p.u.)
    downhill::BitMatrix         # d_j >= 0 on the incoming line
    passive::BitVector          # no battery or DER in the subtree
    psubs_min::Vector{Float64}  # lossless lower bound on P_Subs, per period
end

function voltage_screen(data, T::Int)
    buses, root = data[:Nset], data[:substationBus]
    children, parent = data[:children], data[:parent]
    pos = Dict(j => n for (n, j) in enumerate(buses))
    order = Int[]                                   # parents before children
    stack = [root]
    while !isempty(stack)
        k = pop!(stack); push!(order, k)
        append!(stack, children[k])
    end
    B, D, NL = Set(data[:Bset]), Set(data[:Dset]), Set(data[:NLset])
    nN = length(buses)
    passive = falses(nN)
    for k in reverse(order)
        passive[pos[k]] = !(k in B) && !(k in D) && all(passive[pos[c]] for c in children[k])
    end
    vmin(j) = data[:Vminpu][j]^2; vmax(j) = data[:Vmaxpu][j]^2
    vub = zeros(nN, T); downhill = falses(nN, T)
    keep_lo = trues(nN, T); keep_up = trues(nN, T)
    pmin = zeros(nN); qmin = zeros(nN)
    psubs_min = zeros(T)
    for t in 1:T
        for k in reverse(order)
            n = pos[k]
            pmin[n] = (k in NL ? data[:p_L_pu][k,t] : 0.0) - (k in D ? data[:p_D_pu][k,t] : 0.0) -
                      (k in B ? data[:P_B_R_pu][k] : 0.0)
            qmin[n] = (k in NL ? data[:q_L_pu][k,t] : 0.0) -
                      (k in D ? sqrt(max(0.0, data[:S_D_R][k]^2 - data[:p_D_pu][k,t]^2)) : 0.0)
            for c in children[k]
                pmin[n] += pmin[pos[c]]; qmin[n] += qmin[pos[c]]
            end
        end
        nr = pos[root]
        psubs_min[t] = sum(pmin[pos[c]] for c in children[root]; init=0.0)
        vub[nr, t] = 1.05^2
        keep_lo[nr, t] = false; keep_up[nr, t] = false          # fixed by its equality row
        for k in order                                           # upper limits, root to leaves
            k == root && continue
            n, i = pos[k], parent[k]
            e = (i, k)
            (data[:rdict_pu][e] >= 0 && data[:xdict_pu][e] >= 0) || error("negative line impedance on $e")
            d = 2(data[:rdict_pu][e] * pmin[n] + data[:xdict_pu][e] * qmin[n])
            vub[n, t] = vub[pos[i], t] - d
            downhill[n, t] = d >= 0
            if vub[n, t] <= vmax(k)
                keep_up[n, t] = false
            elseif downhill[n, t] && vmax(k) >= vmax(i)
                keep_up[n, t] = false
            end
        end
        for k in reverse(order)                                  # lower limits, leaves to root
            k == root && continue
            n = pos[k]
            any(downhill[pos[c], t] && vmin(c) >= vmin(k) for c in children[k]) && (keep_lo[n, t] = false)
        end
    end
    return VoltageScreen(buses, keep_lo, keep_up, vub, downhill, passive, psubs_min)
end
