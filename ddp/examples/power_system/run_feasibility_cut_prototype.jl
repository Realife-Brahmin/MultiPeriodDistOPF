# Can a Benders feasibility cut represent F_t without ever writing F_t down?
#
# The adaptive penalty says WHETHER a dispatch is servable and by HOW MUCH it
# misses. It does not say which way to move, so an outer layer can only reject.
# A feasibility cut supplies the direction. At an infeasible dispatch P_B0 solve
#
#     v(P_B) = min  (voltage violation)     s.t. network equations, P_B fixed
#
# and linearise:   v0 + g'(P_B - P_B0) <= 0,   g = d v / d P_B = balance duals.
#
# The inner problem is the BFM SOCP -- convex -- so v is a convex function of
# P_B and its linearisation is a GLOBAL under-estimator. The cut therefore
# cannot remove a dispatch the network can actually serve. That is the claim
# under test, and it is the whole reason the cut route is worth having: it is
# derived from duals, carries no tuned constant and no per-system bus list.
#
# Two things are checked, and only one of them is the headline:
#
#   1. VALIDITY (must hold). No cut generated at an infeasible dispatch may
#      exclude any dispatch the hard solve calls feasible. A single violation
#      here kills the approach -- it would mean silently discarding legal
#      operating points.
#   2. STRENGTH (informative). How many OTHER infeasible dispatches does one cut
#      already exclude? A cut that only removes its own point is valid but
#      useless; one that removes several is genuinely learning the shape of F_t.
#
# Plus a finite-difference check on g itself, because `pure_feasibility` changes
# the objective and hence what the balance duals mean -- the sign convention is
# re-verified here rather than assumed from the Phi_t probe.
#
# Run:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/run_feasibility_cut_prototype.jl [system] [T...]

using Printf
using Random
using Serialization
using Statistics

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "inner_network_opf.jl"))

system = length(ARGS) >= 1 ? ARGS[1] : "ieee2522C_1ph"
horizons = length(ARGS) >= 2 ? parse.(Int, ARGS[2:end]) : [3, 6, 12]

outdir = joinpath(REPO, "ddp", "results", "reduced_space"); mkpath(outdir)
cutcsv = joinpath(outdir, "feasibility_cuts_$(system).csv")
evalcsv = joinpath(outdir, "feasibility_cut_eval_$(system).csv")
cio = open(cutcsv, "w")
println(cio, "system,horizon,time_index,demand_label,source_pattern,v0,grad_norm," *
             "n_nonzero_grad,solve_time_s,fd_rel_err")
eio = open(evalcsv, "w")
println(eio, "system,horizon,time_index,source_pattern,target_pattern,target_hard_feasible," *
             "cut_value,excluded,valid")

n_cuts = 0; n_eval = 0
n_invalid = 0                 # cut excluded a hard-FEASIBLE dispatch: fatal
n_excl_infeas = 0; n_infeas_pairs = 0
fd_errs = Float64[]
cut_times = Float64[]

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

    # same 12-pattern grid as the adaptive check, so the verdicts line up
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
        # ground truth for every pattern at this (T, t)
        hardfeas = Dict{String,Bool}()
        for (name, pb) in pats
            hardfeas[name] = inner_opf(data, t, pb).feasible
        end
        infeas = [(n, pb) for (n, pb) in pats if !hardfeas[n]]
        @printf("T=%d t=%2d (%s): %d infeasible of %d\n",
                T, t, label, length(infeas), length(pats)); flush(stdout)
        isempty(infeas) && continue

        for (sname, spb) in infeas
            t0 = time()
            fs = inner_opf(data, t, spb; voltage_penalty = 1.0, pure_feasibility = true)
            dtc = time() - t0
            fs.feasible || (@printf("  cut solve failed at %s\n", sname); continue)
            v0 = fs.total_violation
            grad = copy(fs.lambda_bal)
            push!(cut_times, dtc)

            # finite-difference check of the gradient, in a random direction
            rng = MersenneTwister(777)
            d = [(2rand(rng) - 1) * pbr[b] for b in 1:nB]
            h = 1e-4
            fp = inner_opf(data, t, spb .+ h .* d; voltage_penalty = 1.0, pure_feasibility = true)
            fm = inner_opf(data, t, spb .- h .* d; voltage_penalty = 1.0, pure_feasibility = true)
            fd_rel = NaN
            if fp.feasible && fm.feasible
                fd = (fp.total_violation - fm.total_violation) / (2h)
                an = sum(grad .* d)
                fd_rel = abs(fd - an) / max(abs(fd), 1e-12)
                push!(fd_errs, fd_rel)
            end

            @printf(cio, "%s,%d,%d,%s,%s,%.6e,%.6e,%d,%.3f,%.3e\n",
                    system, T, t, label, sname, v0, sqrt(sum(abs2, grad)),
                    count(x -> abs(x) > 1e-9, grad), dtc, fd_rel)
            global n_cuts += 1

            # evaluate the cut on every pattern at this (T, t)
            for (tname, tpb) in pats
                cv = v0 + sum(grad .* (tpb .- spb))
                excluded = cv > 1e-8
                tf = hardfeas[tname]
                valid = !(excluded && tf)
                global n_eval += 1
                valid || (global n_invalid += 1)
                if !tf && tname != sname
                    global n_infeas_pairs += 1
                    excluded && (global n_excl_infeas += 1)
                end
                @printf(eio, "%s,%d,%d,%s,%s,%s,%.6e,%s,%s\n",
                        system, T, t, sname, tname, tf, cv, excluded, valid)
            end
        end
    end
end
close(cio); close(eio)

@printf("\nCUTS generated: %d\n", n_cuts)
@printf("CUT evaluations: %d\n", n_eval)
@printf("CUT VALIDITY violations (excluded a hard-feasible dispatch): %d\n", n_invalid)
@printf("CUT strength: excluded %d of %d other infeasible dispatches (%.0f%%)\n",
        n_excl_infeas, n_infeas_pairs,
        n_infeas_pairs == 0 ? 0.0 : 100 * n_excl_infeas / n_infeas_pairs)
isempty(fd_errs) || @printf("CUT gradient FD check: max rel err %.3e over %d dirs\n",
                            maximum(fd_errs), length(fd_errs))
isempty(cut_times) || @printf("CUT solve cost: mean %.2f s\n", mean(cut_times))
println("CUT wrote=$cutcsv")
println("CUT wrote=$evalcsv")
