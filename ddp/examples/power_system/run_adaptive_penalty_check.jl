# Does the tuning-free adaptive-rho scheme reproduce the hard-constrained
# feasibility verdict, without any per-system constant?
#
# Ground truth is the hard-limit inner OPF: it either solves (the dispatch is
# network-feasible) or reports LOCALLY_INFEASIBLE. The adaptive scheme never
# sees that answer; it only raises the penalty coefficient until the violation
# either vanishes or stops moving. This checks the two must-holds:
#
#   * verdict agreement  -- adaptive must call feasible exactly what the hard
#     solve calls feasible, with no false "feasible" (which would be a silent
#     voltage violation) and no false "infeasible" (which would discard a
#     legitimate dispatch);
#   * value agreement    -- where both say feasible, Phi_t must match.
#
# It also records the price: how many inner solves the adaptive scheme costs.
#
# Run:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/run_adaptive_penalty_check.jl [system] [T...]

using Printf
using Random
using Serialization
using Statistics

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "inner_network_opf.jl"))

system = length(ARGS) >= 1 ? ARGS[1] : "ieee2522C_1ph"
horizons = length(ARGS) >= 2 ? parse.(Int, ARGS[2:end]) : [3, 6, 12]

outdir = joinpath(REPO, "ddp", "results", "reduced_space"); mkpath(outdir)
csvpath = joinpath(outdir, "adaptive_penalty_check_$(system).csv")
io = open(csvpath, "w")
println(io, "system,horizon,time_index,demand_label,pattern,hard_status,hard_feasible," *
            "adaptive_verdict,adaptive_feasible,verdict_agrees,hard_phi_usd,adaptive_phi_usd," *
            "phi_gap_usd,final_rho,rounds,settled_violation,total_inner_solves,total_time_s")

agree = 0; total = 0; false_feas = 0; false_infeas = 0; inconclusive = 0
rounds_all = Int[]; gaps = Float64[]

for T in horizons
    datafile = joinpath(REPO, "ddp", "results", "network_filterddp",
                        "network_data_$(system)_T$(T).jls")
    isfile(datafile) || (@printf("SKIP T=%d\n", T); continue)
    data = deserialize(datafile)
    Bset = data[:Bset]; nB = length(Bset)
    pbr = Float64[data[:P_B_R_pu][j] for j in Bset]
    net = [sum(data[:p_L_pu][j, t] for j in data[:NLset]) -
           sum(data[:p_D_pu][j, t] for j in data[:Dset]) for t in 1:T]
    ord = sortperm(net)
    samples = [(ord[1], "low"), (ord[cld(length(ord), 2)], "medium"), (ord[end], "high")]

    root = data[:substationBus]
    depth = Dict{Int,Int}(root => 0); q = [root]
    while !isempty(q)
        j = popfirst!(q)
        for k in get(data[:children], j, Int[]); depth[k] = depth[j] + 1; push!(q, k); end
    end
    bd = [depth[j] for j in Bset]; md = median(bd)
    deep = findall(>=(md), bd); shallow = findall(<(md), bd)

    pats = Tuple{String,Vector{Float64}}[
        ("zero", zeros(nB)), ("all_max_charge", -copy(pbr)),
        ("all_max_discharge", copy(pbr))]
    v = zeros(nB); v[deep] .= -pbr[deep]; push!(pats, ("deep_max_charge", copy(v)))
    v = zeros(nB); v[shallow] .= -pbr[shallow]; push!(pats, ("shallow_max_charge", copy(v)))
    v = zeros(nB); v[deep] .= -pbr[deep]; v[shallow] .= pbr[shallow]
    push!(pats, ("opposing_deep_chg_shallow_dis", copy(v)))
    for s in 1:6
        rng = MersenneTwister(2000 + s)
        push!(pats, ("random$(s)", [(2rand(rng) - 1) * pbr[b] for b in 1:nB]))
    end

    for (t, label) in samples
        for (name, pb) in pats
            hard = inner_opf(data, t, pb)
            t0 = time()
            ad = inner_opf_adaptive(data, t, pb)
            dt_ad = time() - t0
            hf = hard.feasible
            af = ad.network_feasible
            ok = (af === missing) ? false : (af == hf)
            gap = (hf && af === true) ? abs(hard.substation_cost - ad.res.substation_cost) : NaN
            global total += 1
            ok && (global agree += 1)
            if af === missing
                global inconclusive += 1
            elseif af && !hf
                global false_feas += 1
            elseif !af && hf
                global false_infeas += 1
            end
            push!(rounds_all, ad.rounds)
            isnan(gap) || push!(gaps, gap)
            @printf(io, "%s,%d,%d,%s,%s,%s,%s,%s,%s,%s,%.6f,%.6f,%.3e,%.1e,%d,%.6e,%d,%.3f\n",
                system, T, t, label, name, hard.status, hf, ad.verdict,
                af === missing ? "missing" : string(af), ok,
                hf ? hard.substation_cost : NaN,
                ad.res.feasible ? ad.res.substation_cost : NaN,
                gap, ad.rho, ad.rounds,
                ad.res.feasible ? ad.res.total_violation : NaN,
                ad.rounds, dt_ad)
        end
        @printf("  T=%d t=%2d (%s) done\n", T, t, label); flush(stdout)
    end
end
close(io)

@printf("\nADAPTIVE verdict agreement: %d/%d\n", agree, total)
@printf("ADAPTIVE false 'feasible' (silent voltage violation): %d\n", false_feas)
@printf("ADAPTIVE false 'infeasible' (discarded a legal dispatch): %d\n", false_infeas)
@printf("ADAPTIVE inconclusive: %d\n", inconclusive)
@printf("ADAPTIVE inner solves per decision: min=%d max=%d mean=%.2f\n",
        minimum(rounds_all), maximum(rounds_all), mean(rounds_all))
isempty(gaps) || @printf("ADAPTIVE max |Phi gap| where both feasible: %.3e USD\n", maximum(gaps))
println("ADAPTIVE wrote=$csvpath")
