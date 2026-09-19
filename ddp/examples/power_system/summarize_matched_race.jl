# Ipopt vs FilterDDP (diagonal Hessian) on identical problems, one row per
# (system, T). Pure parsing of run_matched_ipopt_race.sh logs, plus the paper's
# old-family centralized sweep as a clearly labelled reference column.
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/summarize_matched_race.jl <race_dir>
#
# Writes <race_dir>/matched_race.csv and prints the table.

using Printf

race = ARGS[1]
const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
logdir = joinpath(race, "logs")

kv(s) = Dict(m.captures[1] => m.captures[2] for m in eachmatch(r"(\w+)=(\S+)", s))

function ipopt_row(sys, T)
    p = joinpath(logdir, "ipopt_$(sys)_T$(T).log")
    isfile(p) || return nothing
    txt = read(p, String)
    m = match(r"CENTRAL_IPOPT [^\n]*", txt)
    m === nothing && return (status = occursin("PIPELINE_STEP_OK", txt) ? "CRASHED" : "RUNNING",
                             iters = -1, solve = NaN, obj = NaN, gamma = NaN)
    d = kv(m.match)
    return (status = d["status"], iters = parse(Int, d["iterations"]),
            solve = parse(Float64, d["solve_time_s"]), obj = parse(Float64, d["objective"]),
            gamma = parse(Float64, get(d, "gamma", "NaN")))
end

function fddp_row(sys, T)
    p = joinpath(logdir, "fddp_diag_$(sys)_T$(T).log")
    isfile(p) || return nothing
    txt = read(p, String)
    m = match(r"solve complete: ([\d.]+) s, iterations=(\d+), status=(\d+)", txt)
    o = match(r"FilterDDP objective=([-\d.eE+]+)", txt)
    g = match(r"TERMINAL_SOC soft=\d gamma=([-\d.eE+]+)", txt)
    fb = match(r"factor_backed=(\d)", txt)
    base = (gamma = g === nothing ? NaN : parse(Float64, g.captures[1]),
            fb = fb === nothing ? 0 : parse(Int, fb.captures[1]))
    m === nothing && return (status = occursin("PIPELINE_STEP_OK", txt) ? "CRASHED" : "RUNNING",
                             iters = -1, wall = NaN, obj = NaN, base...)
    return (status = m.captures[3] == "0" ? "OK" : "FAIL" * m.captures[3],
            iters = parse(Int, m.captures[2]), wall = parse(Float64, m.captures[1]),
            obj = o === nothing ? NaN : parse(Float64, o.captures[1]), base...)
end

old = Dict{Tuple{String,Int},Tuple{Float64,Int}}()
oldcsv = joinpath(REPO, "ddp", "results", "centralized_ipopt", "centralized_ipopt_timing.csv")
if isfile(oldcsv)
    lines = readlines(oldcsv)
    hdr = [strip(h, '"') for h in split(lines[1], ',')]
    ix(n) = findfirst(==(n), hdr)
    for l in lines[2:end]
        f = [strip(x, '"') for x in split(l, ',')]
        old[(f[ix("system")], parse(Int, f[ix("T")]))] =
            (parse(Float64, f[ix("ipopt_reported_s")]), parse(Int, f[ix("iterations")]))
    end
end

keys_ = Tuple{String,Int}[]
for f in (isdir(logdir) ? readdir(logdir) : String[])
    m = match(r"^(?:ipopt|fddp_diag)_(\w+?)_T(\d+)\.log$", f)
    m === nothing || push!(keys_, (m.captures[1], parse(Int, m.captures[2])))
end
sysorder = Dict("ieee123C_1ph" => 1, "ieee2522C_1ph" => 2, "large10kC_1ph" => 3)
sort!(unique!(keys_), by = k -> (get(sysorder, k[1], 9), k[2]))

f1(x) = isnan(x) ? "--" : @sprintf("%.1f", x)

open(joinpath(race, "matched_race.csv"), "w") do io
    println(io, "system,T,gamma,ipopt_status,ipopt_iterations,ipopt_solve_s,ipopt_objective,diag_status,diag_iterations,diag_wall_s,diag_objective,factor_backed,obj_rel_diff,diag_over_ipopt,oldfamily_ipopt_s,oldfamily_iterations")
    println("=" ^ 118)
    println("IDENTICAL PROBLEMS: centralized Ipopt vs FilterDDP (diagonal Hessian)")
    println("periodic profile, C_B = 1e-3, soft terminal SOC with per-system gamma (terminal_soc_penalty.jl)")
    println("=" ^ 118)
    @printf("%-14s %4s | %-15s %5s %9s | %-7s %5s %9s | %9s %9s | %15s\n",
            "system", "T", "ipopt", "iters", "ipopt s", "diag", "iters", "diag s",
            "obj rel", "diag/ip", "old-family ref")
    println("-" ^ 118)
    for (sys, T) in keys_
        I = ipopt_row(sys, T); D = fddp_row(sys, T)
        rel = (I !== nothing && D !== nothing && !isnan(I.obj) && !isnan(D.obj)) ?
              abs(I.obj - D.obj) / abs(I.obj) : NaN
        ratio = (I !== nothing && D !== nothing && D.status == "OK" && !isnan(I.solve)) ?
                D.wall / I.solve : NaN
        o = get(old, (sys, T), (NaN, -1))
        gam = I !== nothing ? I.gamma : (D !== nothing ? D.gamma : NaN)
        @printf(io, "%s,%d,%.6e,%s,%d,%.3f,%.12g,%s,%d,%.3f,%.12g,%d,%.3e,%.4f,%.3f,%d\n",
                sys, T, gam,
                I === nothing ? "" : I.status, I === nothing ? -1 : I.iters,
                I === nothing ? NaN : I.solve, I === nothing ? NaN : I.obj,
                D === nothing ? "" : D.status, D === nothing ? -1 : D.iters,
                D === nothing ? NaN : D.wall, D === nothing ? NaN : D.obj,
                D === nothing ? 0 : D.fb, rel, ratio, o[1], o[2])
        @printf("%-14s %4d | %-15s %5s %9s | %-7s %5s %9s | %9s %9s | %15s\n",
                sys, T,
                I === nothing ? "--" : first(I.status, 15),
                I === nothing || I.iters < 0 ? "--" : string(I.iters),
                I === nothing ? "--" : f1(I.solve),
                D === nothing ? "--" : D.status,
                D === nothing || D.iters < 0 ? "--" : string(D.iters),
                D === nothing ? "--" : f1(D.wall),
                isnan(rel) ? "--" : @sprintf("%.1e", rel),
                isnan(ratio) ? "--" : @sprintf("%.2fx", ratio),
                isnan(o[1]) ? "--" : @sprintf("%.1f s/%d it", o[1], o[2]))
    end
end
println("-" ^ 118)
println("ipopt s = Ipopt solve time. diag s = FilterDDP wall. diag/ip < 1.00x means FilterDDP is faster.")
println("obj rel = |Ipopt - FilterDDP| / |Ipopt|: small means the two solved the same problem.")
println("old-family ref = the paper's centralized sweep (C_B ~ 8.8e-8, old sampling, free terminal SOC):")
println("  a DIFFERENT problem, shown only to locate the knee it exhibited.")
