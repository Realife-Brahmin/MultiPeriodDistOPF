# Time-to-near-optimality, recovered from logs of runs that went to full
# convergence (or were stopped). The scoring rule (user, 2026-09-19): a run counts
# the time until its iterate is FEASIBLE (primal infeasibility <= P, default
# 1e-6) AND within 0.5% of the centralized Ipopt objective. The strict-tolerance
# tail does not count. Future runs stop there on their own
# (FILTERDDP_NEAR_OPT_REFERENCE in DDP4OPF's solve.jl); this recovers the same
# number for runs made before that existed.
#
# Timing is exact: FILTERDDP_ITER_TIMING gives every backward/forward pass, and
# iterate k's residuals are known at the END of iteration k's backward pass
# (that is where the trace row is printed), so
#   t_k = setup + sum(total_s of every earlier pass) + backward_s of pass k,
# setup = (reported solve time) - sum(all total_s), about 3 s.
#
#   race mode:  julia near_opt_from_logs.jl race  <race_dir>
#               reference = Ipopt's objective on the identical problem
#   sweep mode: julia near_opt_from_logs.jl sweep <sweep_dir>
#               reference = the exact arm's final objective (same problem as the
#               diag arm; the horizon sweep predates the soft terminal SOC and has
#               no matched Ipopt run)

using Printf

const GAP = 0.005
const PRIMALS = (1e-4, 1e-5, 1e-6)

function trace_and_times(path)
    txt = read(path, String)
    trace = Dict{Int,Tuple{Float64,Float64}}()      # k => (objective, primal_inf)
    order = Int[]
    for m in eachmatch(r"(?m)^\s*(\d+)\s+(-?[\d.]+e[+-]\d+)\s+([\d.]+e[+-]\d+)\s+[\d.]+e[+-]\d+", txt)
        k = parse(Int, m.captures[1])
        haskey(trace, k) && continue
        trace[k] = (parse(Float64, m.captures[2]), parse(Float64, m.captures[3]))
        push!(order, k)
    end
    rows = [(parse(Int, m.captures[1]), parse(Float64, m.captures[2]),
             parse(Float64, m.captures[3]), m.captures[4])
            for m in eachmatch(r"FILTERDDP_ITER_TIMING iteration=(\d+) \S+ backward_s=([\d.]+) \S+ total_s=([\d.]+) outcome=(\w+)", txt)]
    sc = match(r"solve complete: ([\d.]+) s", txt)
    total = isempty(rows) ? 0.0 : sum(r[3] for r in rows)
    setup = sc === nothing ? 0.0 : max(0.0, parse(Float64, sc.captures[1]) - total)
    tk = Dict{Int,Float64}()
    cum = setup
    for (k, bwd, tot, outcome) in rows
        outcome != "barrier_update" && !haskey(tk, k) && (tk[k] = cum + bwd)
        cum += tot
    end
    wall = sc === nothing ? cum : parse(Float64, sc.captures[1])
    return trace, order, tk, wall, sc !== nothing
end

function first_near_opt(trace, order, tk, ref, P)
    for k in order
        obj, pr = trace[k]
        if pr <= P && abs(obj - ref) <= GAP * abs(ref) && haskey(tk, k)
            return k, tk[k], abs(obj - ref) / max(abs(ref), eps(Float64)), obj
        end
    end
    return -1, NaN, NaN, NaN
end

fmt(x) = isnan(x) ? "--" : @sprintf("%.1f", x)

mode, dir = ARGS[1], ARGS[2]
logdir = joinpath(dir, "logs")
out = joinpath(dir, "near_opt.csv")
cases = Tuple{String,Int,String,String}[]      # system, T, diag log, reference source
sysorder = Dict("ieee123C_1ph" => 1, "ieee2522C_1ph" => 2, "large10kC_1ph" => 3)
if mode == "race"
    for f in readdir(logdir)
        m = match(r"^fddp_diag_(\w+?)_T(\d+)\.log$", f)
        m === nothing && continue
        push!(cases, (m.captures[1], parse(Int, m.captures[2]), joinpath(logdir, f),
                      joinpath(logdir, "ipopt_$(m.captures[1])_T$(m.captures[2]).log")))
    end
else
    for f in readdir(logdir)
        m = match(r"^(\w+?)_T(\d+)_diag\.log$", f)
        m === nothing && continue
        push!(cases, (m.captures[1], parse(Int, m.captures[2]), joinpath(logdir, f),
                      joinpath(logdir, "$(m.captures[1])_T$(m.captures[2])_exact.log")))
    end
end
sort!(cases, by = c -> (get(sysorder, c[1], 9), c[2]))

refname = mode == "race" ? "ipopt" : "exact"
println("=" ^ 112)
println("TIME TO NEAR-OPTIMALITY: first iterate with primal_inf <= P and objective within 0.5% of the $(refname) reference")
println("=" ^ 112)
@printf("%-14s %4s | %10s | %22s %22s %22s | %9s\n", "system", "T",
        "$(refname) s", "P=1e-4: it / s", "P=1e-5: it / s", "P=1e-6: it / s", "ratio@1e-6")
println("-" ^ 112)
open(out, "w") do io
    println(io, "mode,system,T,reference_objective,reference_time_s,full_run_wall_s,full_run_finished,P,near_opt_iteration,near_opt_time_s,near_opt_objective,near_opt_rel_gap,ratio_to_reference")
    for (sys, T, dlog, rlog) in cases
        isfile(rlog) || continue
        rtxt = read(rlog, String)
        if mode == "race"
            mo = match(r"CENTRAL_IPOPT .* objective=([-\d.eE+]+) solve_time_s=([\d.]+)", rtxt)
            mo === nothing && continue
            ref, reft = parse(Float64, mo.captures[1]), parse(Float64, mo.captures[2])
        else
            mo = match(r"FilterDDP objective=([-\d.eE+]+)", rtxt)
            ms = match(r"solve complete: ([\d.]+) s", rtxt)
            (mo === nothing || ms === nothing) && continue
            ref, reft = parse(Float64, mo.captures[1]), parse(Float64, ms.captures[1])
            # the exact arm's own time to near-optimality is the fair comparison
            tr, od, tk, _, _ = trace_and_times(rlog)
            reft = first_near_opt(tr, od, tk, ref, 1e-6)[2]
        end
        trace, order, tk, wall, finished = trace_and_times(dlog)
        cells = String[]
        k6 = -1; t6 = NaN
        for P in PRIMALS
            k, t, g, obj = first_near_opt(trace, order, tk, ref, P)
            push!(cells, k < 0 ? "not reached" : @sprintf("%d / %.1f", k, t))
            P == 1e-6 && (k6 = k; t6 = t)
            @printf(io, "%s,%s,%d,%.12g,%.3f,%.3f,%d,%.0e,%d,%.3f,%.12g,%.3e,%.4f\n",
                    mode, sys, T, ref, reft, wall, finished, P, k, t, obj, g, t / reft)
        end
        @printf("%-14s %4d | %10s | %22s %22s %22s | %9s\n", sys, T, fmt(reft),
                cells..., isnan(t6) ? "--" : @sprintf("%.2fx", t6 / reft))
    end
end
println("-" ^ 112)
mode == "race" ?
    println("ipopt s = Ipopt's full solve on the identical problem. ratio = FilterDDP time to near-optimality / Ipopt time.") :
    println("exact s = the EXACT arm's own time to near-optimality (P = 1e-6). ratio = diag / exact, both to near-optimality.")
println("wrote $out")
