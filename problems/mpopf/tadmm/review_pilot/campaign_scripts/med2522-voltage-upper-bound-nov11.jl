using Serialization, Printf

function voltage_upper_bound(path)
    data = deserialize(path)
    buses = data[:Nset]
    root = data[:substationBus]
    children = data[:children]
    parent = data[:parent]
    r, x = data[:rdict_pu], data[:xdict_pu]
    @assert all(e -> r[e] >= 0 && x[e] >= 0, data[:Lset])
    loads, batteries, ders = Set(data[:NLset]), Set(data[:Bset]), Set(data[:Dset])
    order = [root]
    i = 1
    while i <= length(order)
        append!(order, children[order[i]])
        i += 1
    end
    @assert Set(order) == Set(buses) && length(order) == length(buses)
    pos = Dict(bus => i for (i,bus) in enumerate(order))

    # Relax storage energy constraints and grant all batteries maximum discharge
    # and all DERs maximum reactive injection. Series losses on every branch
    # can only lower voltage relative to the resulting lossless voltage.
    best = (margin=Inf, bus=0, period=0, upper=NaN, lower=NaN)
    for t in 1:data[:T]
        p = zeros(length(order))
        q = zeros(length(order))
        v = fill(1.05^2, length(order))
        for k in length(order):-1:2
            j = order[k]
            p[k] += (j in loads ? data[:p_L_pu][j,t] : 0.0) -
                    (j in batteries ? data[:P_B_R_pu][j] : 0.0) -
                    (j in ders ? data[:p_D_pu][j,t] : 0.0)
            q[k] += (j in loads ? data[:q_L_pu][j,t] : 0.0) -
                    (j in ders ? sqrt(max(0.0, data[:S_D_R][j]^2-data[:p_D_pu][j,t]^2)) : 0.0)
            ip = pos[parent[j]]
            p[ip] += p[k]
            q[ip] += q[k]
        end
        for k in 2:length(order)
            j = order[k]
            e = (parent[j],j)
            v[k] = v[pos[parent[j]]] - 2*(r[e]*p[k]+x[e]*q[k])
            lower = data[:Vminpu][j]^2
            margin = v[k]-lower
            if margin < best.margin
                best = (;margin,bus=j,period=t,upper=v[k],lower)
            end
        end
    end
    @printf("bus=%d period=%d upper_v_squared=%.12f required_v_squared=%.12f margin=%.12f\n",
        best.bus,best.period,best.upper,best.lower,best.margin)
    println(best.margin < -1e-8 ? "Voltage lower bound infeasibility certified." :
        "This relaxation does not certify infeasibility.")
end

voltage_upper_bound(only(ARGS))
