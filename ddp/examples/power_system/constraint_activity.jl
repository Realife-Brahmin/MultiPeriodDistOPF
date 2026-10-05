# Which inequality constraints are active at the MPOPF optimum, and which
# could have been known inactive beforehand? (Constraint screening, agenda of
# 2026-10-07.) Observation only: nothing is removed from any solve here.
#
# Solves the matched centralized model with Ipopt (accurate, independent of
# FilterDDP) by including centralized_ipopt_matched.jl, then
#
#   1. counts every inequality by type and how many are active at the optimum.
#      Active: slack at most 1e-6 of the feasible range, or multiplier larger
#      than slack (the interior-point identification). Near: not active, slack
#      within 1% of the range. Inactive: the rest.
#
#   2. evaluates the exact screening rules of voltage_screening.jl (they need
#      no OPF solution): the substation voltage is fixed, a lossless upper
#      bound on every voltage, and "voltage cannot rise along a line whose
#      subtree cannot export". It counts what the rules keep and checks their
#      premises against the solution.
#
#   3. writes the masks (which voltage limits are needed at the optimum, and
#      which each rule keeps) for the FilterDDP driver's screening switch.
#
#   REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1 \
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/constraint_activity.jl <system> <T> [ipopt_log]

include(joinpath(@__DIR__, "centralized_ipopt_matched.jl"))   # builds and solves `model`

const ACTIVE_TOL = 1e-6
const NEAR_TOL = 1e-2

lo_mult(x) = abs(dual(LowerBoundRef(x)))
up_mult(x) = abs(dual(UpperBoundRef(x)))
# 0 inactive, 1 near, 2 active
classify(slack, scale, mult) = (slack <= ACTIVE_TOL * scale || mult > slack) ? 2 : (slack <= NEAR_TOL * scale ? 1 : 0)

struct Tally; total::Int; active::Int; near::Int; elements::Int; never::Int; end
function tally(class::Function, elements, Tset)
    total = active = near = never = 0
    for e in elements
        ever = false
        for t in Tset
            c = class(e, t)
            total += 1
            c == 2 && (active += 1; ever = true)
            c == 1 && (near += 1)
        end
        ever || (never += 1)
    end
    return Tally(total, active, near, length(elements), never)
end

V, PBv, Bv, QDv = value.(v), value.(P_B), value.(B), value.(q_D)
ELLv, Pv, Qv, PSv = value.(ell), value.(P), value.(Q), value.(P_Subs)
vmin(j) = data[:Vminpu][j]^2; vmax(j) = data[:Vmaxpu][j]^2; vrange(j) = vmax(j) - vmin(j)
qmax(j, t) = sqrt(max(0.0, data[:S_D_R][j]^2 - data[:p_D_pu][j,t]^2))
pbr(j) = data[:P_B_R_pu][j]
bmin(j) = data[:soc_min][j] * data[:B_R_pu][j]; bmax(j) = data[:soc_max][j] * data[:B_R_pu][j]
soc_con = Dict{Tuple{Any,Int},Any}()
let cons = all_constraints(model, QuadExpr, MOI.LessThan{Float64}), k = 0
    length(cons) == length(Lset) * T || error("unexpected number of SOC rows")
    for t in Tset, e in Lset                    # creation order in centralized_ipopt_matched.jl
        soc_con[(e, t)] = cons[k += 1]
    end
end
soc_slack(e, t) = V[e[1],t] * ELLv[e,t] - Pv[e,t]^2 - Qv[e,t]^2

nonroot_buses = [j for j in Nset if j != root]
v_lo(j, t) = classify(V[j,t] - vmin(j), vrange(j), lo_mult(v[j,t]))
v_up(j, t) = classify(vmax(j) - V[j,t], vrange(j), up_mult(v[j,t]))
rows = Pair{String,Tally}[
    "voltage lower limit"            => tally(v_lo, nonroot_buses, Tset),
    "voltage upper limit"            => tally(v_up, nonroot_buses, Tset),
    "battery power, charge limit"    => tally((j,t) -> classify(PBv[j,t] + pbr(j), 2pbr(j), lo_mult(P_B[j,t])), Bset, Tset),
    "battery power, discharge limit" => tally((j,t) -> classify(pbr(j) - PBv[j,t], 2pbr(j), up_mult(P_B[j,t])), Bset, Tset),
    "battery energy, lower limit"    => tally((j,t) -> classify(Bv[j,t] - bmin(j), bmax(j) - bmin(j), lo_mult(B[j,t])), Bset, Tset),
    "battery energy, upper limit"    => tally((j,t) -> classify(bmax(j) - Bv[j,t], bmax(j) - bmin(j), up_mult(B[j,t])), Bset, Tset),
    "DER reactive, lower limit"      => tally((j,t) -> classify(QDv[j,t] + qmax(j,t), 2qmax(j,t), lo_mult(q_D[j,t])), Dset, Tset),
    "DER reactive, upper limit"      => tally((j,t) -> classify(qmax(j,t) - QDv[j,t], 2qmax(j,t), up_mult(q_D[j,t])), Dset, Tset),
    # ell is tiny in p.u., so the multiplier test misfires on it: call ell >= 0
    # active when the line carries no power (below 1 W). It is implied by the
    # SOC row and v > 0 in any case.
    "line current ell >= 0"          => tally((e,t) -> Pv[e,t]^2 + Qv[e,t]^2 <= 1e-12 ? 2 : 0, Lset, Tset),
    "SOC relaxation P^2+Q^2 <= v*ell" => tally((e,t) -> classify(soc_slack(e,t), 0.0, abs(dual(soc_con[(e,t)]))), Lset, Tset),
    "substation P_Subs >= 0"         => tally((e,t) -> classify(PSv[t], 0.0, lo_mult(P_Subs[t])), [1], Tset),
]

@printf("CONSTRAINT_ACTIVITY system=%s T=%d status=%s buses=%d lines=%d batteries=%d ders=%d\n",
        system, T, string(termination_status(model)), length(Nset), length(Lset), length(Bset), length(Dset))
@printf("  %-34s %9s %9s %7s %9s %9s | %8s %14s\n", "inequality", "total", "active", "share", "near(1%)", "inactive", "elements", "never active")
let tot = 0, act = 0, nr = 0
    for (name, r) in rows
        @printf("  %-34s %9d %9d %6.1f%% %9d %9d | %8d %8d (%5.1f%%)\n", name, r.total, r.active,
                100r.active / max(r.total, 1), r.near, r.total - r.active - r.near, r.elements, r.never, 100r.never / max(r.elements, 1))
        tot += r.total; act += r.active; nr += r.near
    end
    @printf("  %-34s %9d %9d %6.1f%% %9d %9d\n", "ALL INEQUALITIES", tot, act, 100act / tot, nr, tot - act - nr)
end
@printf("VOLTAGE_RANGE lowest %.4f pu (limit %.3f), highest away from the substation %.4f pu (limit %.3f), substation fixed at 1.05 pu\n",
        sqrt(minimum(V[j,t] for j in nonroot_buses, t in Tset)), minimum(data[:Vminpu][j] for j in Nset),
        sqrt(maximum(V[j,t] for j in nonroot_buses, t in Tset)), maximum(data[:Vmaxpu][j] for j in Nset))

let loose = [(e, t) for e in Lset, t in Tset if classify(soc_slack(e,t), 0.0, abs(dual(soc_con[(e,t)]))) < 2]
    zero_z = count(et -> data[:rdict_pu][et[1]] == 0, loose)
    zero_flow = count(et -> Pv[et[1],et[2]]^2 + Qv[et[1],et[2]]^2 <= 1e-12, loose)
    @printf("SOC_LOOSE %d rows on %d lines: zero resistance %d, zero flow %d, largest resistance %.2e pu, largest slack %.2e\n", length(loose),
            length(unique(first.(loose))), zero_z, zero_flow, isempty(loose) ? 0.0 : maximum(data[:rdict_pu][et[1]] for et in loose),
            isempty(loose) ? 0.0 : maximum(soc_slack(et[1], et[2]) for et in loose))
end

# ---- rules that need no OPF solution (voltage_screening.jl) ----------------
vs = voltage_screen(data, T)                 # voltage_screening.jl, included by the Ipopt script
children, parent = data[:children], data[:parent]
posN = Dict(j => n for (n, j) in enumerate(vs.buses))
nN = length(nonroot_buses)
rowsN = [posN[j] for j in nonroot_buses]                       # drop the substation row
keep_lo, keep_up = vs.keep_lo[rowsN, :], vs.keep_up[rowsN, :]
needed_lo = falses(nN, T); needed_up = falses(nN, T)           # active or near at the optimum
active_lo = falses(nN, T); active_up = falses(nN, T)
for (n, j) in enumerate(nonroot_buses), t in Tset
    a, b = v_lo(j,t), v_up(j,t)
    needed_lo[n,t] = a > 0; needed_up[n,t] = b > 0
    active_lo[n,t] = a == 2; active_up[n,t] = b == 2
end
# the rules' premises, checked on the solution
above_bound = count(V[j,t] > vs.vub[posN[j],t] + 1e-7 for j in nonroot_buses, t in Tset)
downhill_n = count(vs.downhill[rowsN, :])
uphill_seen = count(vs.downhill[posN[j],t] && V[j,t] - V[parent[j],t] > 1e-7 for j in nonroot_buses, t in Tset)
passive_n = count(vs.passive[rowsN])
pct(a, b) = 100a / max(b, 1)
@printf("SCREEN_RULES passive buses %d of %d (%.1f%%) | downhill lines %d of %d bus-periods (%.1f%%) | premises broken by the solution: above the upper bound %d, rise on a downhill line %d\n",
        passive_n, nN, pct(passive_n, nN), downhill_n, nN * T, pct(downhill_n, nN * T), above_bound, uphill_seen)
@printf("SCREEN_RULES lower limits kept %d of %d (%.1f%%), upper limits kept %d of %d (%.1f%%)\n",
        count(keep_lo), nN * T, pct(count(keep_lo), nN * T), count(keep_up), nN * T, pct(count(keep_up), nN * T))
bus_lo = count(any(keep_lo, dims=2)); bus_up = count(any(keep_up, dims=2))
@printf("SCREEN_RULES per bus (a limit is kept if any period needs it): lower %d of %d (%.1f%%), upper %d of %d (%.1f%%)\n",
        bus_lo, nN, pct(bus_lo, nN), bus_up, nN, pct(bus_up, nN))
@printf("SCREEN_RULES limits at the optimum that the rules drop as implied by a kept one: active %d lower, %d upper; near %d lower, %d upper\n",
        count(active_lo .& .!keep_lo), count(active_up .& .!keep_up),
        count(needed_lo .& .!active_lo .& .!keep_lo), count(needed_up .& .!active_up .& .!keep_up))
@printf("ORACLE needed at the optimum (active or near): %d lower, %d upper of %d each | buses: %d lower, %d upper of %d | kept by the rules but not needed: %d lower, %d upper\n",
        count(needed_lo), count(needed_up), nN * T, count(any(needed_lo, dims=2)), count(any(needed_up, dims=2)), nN,
        count(keep_lo .& .!needed_lo), count(keep_up .& .!needed_up))

# Lossless estimate of the lowest voltage: every battery charging at rating and
# every DER absorbing at its reactive limit. NOT a bound (losses lower v further).
order = Int[]
let stack = [root]
    while !isempty(stack)
        k = pop!(stack); push!(order, k)
        append!(stack, children[k])
    end
end
loadp(j, t) = j in data[:NLset] ? data[:p_L_pu][j,t] : 0.0
loadq(j, t) = j in data[:NLset] ? data[:q_L_pu][j,t] : 0.0
vlo_est = Dict{Tuple{Int,Int},Float64}()
for t in Tset
    pnet = Dict{Int,Float64}(); qnet = Dict{Int,Float64}()
    for k in reverse(order)
        pnet[k] = loadp(k,t) - (k in Dset ? data[:p_D_pu][k,t] : 0.0) + (k in Bset ? pbr(k) : 0.0) + sum(pnet[c] for c in children[k]; init=0.0)
        qnet[k] = loadq(k,t) + (k in Dset ? qmax(k,t) : 0.0) + sum(qnet[c] for c in children[k]; init=0.0)
    end
    vlo_est[(root,t)] = 1.05^2
    for k in order
        k == root && continue
        e = (parent[k], k)
        vlo_est[(k,t)] = vlo_est[(parent[k],t)] - 2(data[:rdict_pu][e] * pnet[k] + data[:xdict_pu][e] * qnet[k])
    end
end
@printf("LOWER_ESTIMATE lossless worst-case lowest voltage %.4f pu (realised %.4f, limit %.3f) | bus-periods with estimate below the limit: %d of %d\n",
        sqrt(max(0.0, minimum(values(vlo_est)))), sqrt(minimum(V[j,t] for j in nonroot_buses, t in Tset)),
        minimum(data[:Vminpu][j] for j in Nset), count(vlo_est[(j,t)] < vmin(j) for j in nonroot_buses, t in Tset), nN * T)

# where the needed voltage limits sit
depth = Dict(root => 0); for k in order; k == root || (depth[k] = depth[parent[k]] + 1); end
for (name, need) in (("lower", needed_lo), ("upper", needed_up))
    js = [j for (n, j) in enumerate(nonroot_buses) if any(need[n, :])]
    isempty(js) && continue
    @printf("VOLTAGE_NEEDED %s: %d buses (leaves %d, passive %d, with a device %d), depth %d..%d of %d, periods %s\n", name, length(js),
            count(j -> isempty(children[j]), js), count(j -> vs.passive[posN[j]], js), count(j -> (j in Bset || j in Dset), js),
            minimum(depth[j] for j in js), maximum(depth[j] for j in js), maximum(values(depth)),
            string(findall(t -> any(need[:, t]), 1:T)))
end

# battery limits implied by one another (presolve-style bound implication)
if !isempty(Bset)
    reach_lo(j, t) = data[:B0_pu][j] - t * dt * pbr(j); reach_hi(j, t) = data[:B0_pu][j] + t * dt * pbr(j)
    e_lo = count(reach_lo(j,t) >= bmin(j) for j in Bset, t in Tset); e_hi = count(reach_hi(j,t) <= bmax(j) for j in Bset, t in Tset)
    p_imp = count(dt * pbr(j) >= bmax(j) - bmin(j) for j in Bset)
    @printf("BATTERY_IMPLIED energy limits unreachable from B0 at rated power: %d lower, %d upper of %d each | batteries whose power limit is implied by the energy window: %d of %d\n",
            e_lo, e_hi, length(Bset) * T, p_imp, length(Bset))
end

maskfile = joinpath(dirname(ipopt_log), "voltage_masks_$(system)_T$(T).jls")
serialize(maskfile, (; system, T, buses=nonroot_buses, needed_lo, needed_up, active_lo, active_up, keep_lo, keep_up, objective=obj))
println("VOLTAGE_MASKS ", maskfile)
