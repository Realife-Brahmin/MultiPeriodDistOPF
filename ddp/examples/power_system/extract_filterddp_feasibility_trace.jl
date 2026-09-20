# Convert one FilterDDP stdout log into a per-iterate feasibility trajectory.
# Usage: julia extract_filterddp_feasibility_trace.jl LOG CSV [IPOPT_OBJECTIVE]

using Printf

length(ARGS) >= 2 || error("usage: extract_filterddp_feasibility_trace.jl LOG CSV [IPOPT_OBJECTIVE]")
logpath, outpath = ARGS[1], ARGS[2]
reference = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : NaN
lines = readlines(logpath)

timings = NamedTuple[]
for line in lines
    m = match(r"FILTERDDP_ITER_TIMING iteration=(\d+) barrier_iteration=(\d+) backward_s=([\d.]+) forward_s=([\d.]+) total_s=([\d.]+) outcome=(\w+)", line)
    m === nothing && continue
    push!(timings, (iteration=parse(Int, m.captures[1]), barrier_iteration=parse(Int, m.captures[2]),
        backward=parse(Float64, m.captures[3]), total=parse(Float64, m.captures[5]),
        outcome=m.captures[6]))
end
timed_total = sum(r.total for r in timings; init=0.0)
solve_match = match(r"solve complete: ([\d.]+) s", join(lines, '\n'))
setup = solve_match === nothing ? 0.0 : max(0.0, parse(Float64, solve_match.captures[1]) - timed_total)
elapsed = Dict{Int,Float64}()
barrier_iteration = Dict{Int,Int}()
let running = setup
    for r in timings
        if r.outcome != "barrier_update" && !haskey(elapsed, r.iteration)
            elapsed[r.iteration] = running + r.backward
            barrier_iteration[r.iteration] = r.barrier_iteration
        end
        running += r.total
    end
end

feasibility = Dict{Int,Dict{String,String}}()
for line in lines
    startswith(line, "FILTERDDP_FEASIBILITY ") || continue
    d = Dict(m.captures[1] => m.captures[2] for m in eachmatch(r"(\w+)=(\S+)", line))
    feasibility[parse(Int, d["iteration"])] = d
end

row_pattern = r"^\s*(\d+)\s+([+-]?[\d.]+e[+-]\d+)\s+([+-]?[\d.]+e[+-]\d+)\s+([+-]?[\d.]+e[+-]\d+)\s+([+-]?[\d.]+e[+-]\d+)\s+([+-]?[\d.]+)\s+(\S+)\s+([+-]?[\d.]+e[+-]\d+)\s+(\d+)\s*$"
rows = Dict{Int,NamedTuple}()
for line in lines
    m = match(row_pattern, line)
    m === nothing && continue
    k = parse(Int, m.captures[1])
    obj = parse(Float64, m.captures[2])
    rows[k] = (iteration=k, objective=obj, primal_inf=parse(Float64, m.captures[3]),
        dual_inf=parse(Float64, m.captures[4]), complementarity_inf=parse(Float64, m.captures[5]),
        log10_barrier=parse(Float64, m.captures[6]), step_size=parse(Float64, m.captures[8]),
        backtracks=parse(Int, m.captures[9]))
end

keys_sorted = sort!(collect(keys(rows)))
mkpath(dirname(outpath))
open(outpath, "w") do io
    println(io, "iteration,barrier_iteration,elapsed_algorithm_s,objective,objective_rel_gap,primal_inf,dual_inf,complementarity_inf,barrier_mu,step_size,line_search_backtracks,equality_count,equality_rms,equality_max,equality_worst_stage,equality_worst_index,dynamics_count,dynamics_rms,dynamics_max,dynamics_worst_stage,dynamics_worst_index,bound_count,bound_rms,bound_max,bound_worst_stage,bound_worst_index,bound_worst_kind")
    for k in keys_sorted
        r = rows[k]
        f = get(feasibility, k, Dict{String,String}())
        value(name, default="") = get(f, name, default)
        gap = isnan(reference) ? NaN : abs(r.objective-reference) / max(abs(reference), eps(Float64))
        @printf(io, "%d,%d,%.9f,%.12g,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
            k, get(barrier_iteration, k, -1), get(elapsed, k, NaN), r.objective, gap,
            r.primal_inf, r.dual_inf, r.complementarity_inf, 10.0^r.log10_barrier,
            r.step_size, r.backtracks,
            value("equality_count"), value("equality_rms"), value("equality_max"),
            value("equality_worst_stage"), value("equality_worst_index"),
            value("dynamics_count"), value("dynamics_rms"), value("dynamics_max"),
            value("dynamics_worst_stage"), value("dynamics_worst_index"),
            value("bound_count"), value("bound_rms"), value("bound_max"),
            value("bound_worst_stage"), value("bound_worst_index"), value("bound_worst_kind"))
    end
end
println("wrote $outpath ($(length(keys_sorted)) unique iterations)")
