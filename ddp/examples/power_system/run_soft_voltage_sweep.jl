# Soft-voltage (L1 exact penalty) variant of the inner network OPF.
#
# Replaces the hard voltage box with nonnegative violation slacks carrying a
# linear penalty, so `Phi_t` is defined on the WHOLE battery-power box instead
# of only on `F_t`. Two things are being tested, and they are different:
#
#   1. EXACTNESS. Wherever the hard problem is feasible, the penalised problem
#      must return the same solution (the penalty term is then identically
#      zero). If it does not, the coefficient is too small and the "penalty"
#      is silently buying cost reductions with voltage violations.
#   2. COVERAGE. Wherever the hard problem is infeasible, the penalised problem
#      should converge cleanly and report how far outside `F_t` the dispatch is.
#
# A penalised solve that converges is NOT evidence that the dispatch is
# network-feasible. `total_violation > 0` means the original problem is
# infeasible at that point, full stop.
#
# Run:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/run_soft_voltage_sweep.jl [system] [T...]

using Printf
using Random
using Serialization
using Statistics

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "inner_network_opf.jl"))

system = length(ARGS) >= 1 ? ARGS[1] : "ieee2522C_1ph"
horizons = length(ARGS) >= 2 ? parse.(Int, ARGS[2:end]) : [3, 6, 12]
# coefficient sweep: USD per pu^2 of voltage-squared violation
penalties = [1e2, 1e4, 1e6]

outdir = joinpath(REPO, "ddp", "results", "reduced_space"); mkpath(outdir)
csvpath = joinpath(outdir, "soft_voltage_sweep_$(system).csv")
io = open(csvpath, "w")
println(io, "system,horizon,time_index,demand_label,pattern,penalty,hard_status,hard_feasible," *
            "soft_status,soft_feasible,hard_substation_cost,soft_substation_cost," *
            "phi_gap_usd,total_violation,max_violation,n_violated,penalty_cost," *
            "hard_iters,soft_iters,hard_time_s,soft_time_s,exactness_ok")

for T in horizons
    datafile = joinpath(REPO, "ddp", "results", "network_filterddp",
                        "network_data_$(system)_T$(T).jls")
    if !isfile(datafile)
        @printf("SKIP T=%d: no exported network data (%s)\n", T, basename(datafile)); continue
    end
    data = deserialize(datafile)
    Bset = data[:Bset]; nB = length(Bset)
    pbr = Float64[data[:P_B_R_pu][j] for j in Bset]
    net = [sum(data[:p_L_pu][j, t] for j in data[:NLset]) -
           sum(data[:p_D_pu][j, t] for j in data[:Dset]) for t in 1:T]
    ord = sortperm(net)
    samples = T >= 3 ? [(ord[1], "low"), (ord[cld(length(ord), 2)], "medium"), (ord[end], "high")] :
                       [(ord[end], "high")]

    root = data[:substationBus]
    depth = Dict{Int,Int}(root => 0); q = [root]
    while !isempty(q)
        j = popfirst!(q)
        for k in get(data[:children], j, Int[]); depth[k] = depth[j] + 1; push!(q, k); end
    end
    bd = [depth[j] for j in Bset]; md = median(bd)
    deep = findall(>=(md), bd); shallow = findall(<(md), bd)

    pats = Tuple{String,Vector{Float64}}[
        ("zero", zeros(nB)),
        ("all_max_charge", -copy(pbr)),
        ("all_max_discharge", copy(pbr)),
    ]
    v = zeros(nB); v[deep] .= -pbr[deep]; push!(pats, ("deep_max_charge", copy(v)))
    v = zeros(nB); v[shallow] .= -pbr[shallow]; push!(pats, ("shallow_max_charge", copy(v)))
    v = zeros(nB); v[deep] .= -pbr[deep]; v[shallow] .= pbr[shallow]
    push!(pats, ("opposing_deep_chg_shallow_dis", copy(v)))
    for s in 1:6
        rng = MersenneTwister(2000 + s)
        push!(pats, ("random$(s)", [(2rand(rng) - 1) * pbr[b] for b in 1:nB]))
    end

    @printf("T=%d  samples=%s  patterns=%d\n", T, [s[1] for s in samples], length(pats))
    for (t, label) in samples
        for (name, pb) in pats
            hard = inner_opf(data, t, pb)
            for rho in penalties
                soft = inner_opf(data, t, pb; voltage_penalty = rho)
                gap = (hard.feasible && soft.feasible) ?
                      abs(hard.substation_cost - soft.substation_cost) : NaN
                # exactness: where the hard problem is feasible the penalised
                # solve must reproduce it and report zero violation
                exact = hard.feasible ? (soft.feasible && gap < 1e-4 &&
                                         soft.total_violation < 1e-7) : true
                @printf(io, "%s,%d,%d,%s,%s,%.1e,%s,%s,%s,%s,%.6f,%.6f,%.3e,%.6e,%.6e,%d,%.6f,%d,%d,%.3f,%.3f,%s\n",
                    system, T, t, label, name, rho, hard.status, hard.feasible,
                    soft.status, soft.feasible,
                    hard.feasible ? hard.substation_cost : NaN,
                    soft.feasible ? soft.substation_cost : NaN,
                    gap, soft.feasible ? soft.total_violation : NaN,
                    soft.feasible ? soft.max_violation : NaN,
                    soft.feasible ? soft.n_violated : -1,
                    soft.feasible ? soft.penalty_cost : NaN,
                    hard.iterations, soft.iterations, hard.solve_time, soft.solve_time,
                    exact)
            end
        end
        @printf("  t=%2d (%s) done\n", t, label); flush(stdout)
    end
end
close(io)
println("SOFT_VOLTAGE wrote=$csvpath")
