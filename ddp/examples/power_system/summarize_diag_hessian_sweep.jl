# Parse run_diag_hessian_horizon_sweep.sh logs (plus the T = 3 agenda-pipeline
# logs, when present) into one CSV and one exact-vs-diag table. Pure parsing.
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/summarize_diag_hessian_sweep.jl <sweep_dir> [t3_logdir]
#
# Writes <sweep_dir>/diag_hessian_horizon.csv and prints the table to stdout.

using Printf

sweep = length(ARGS) >= 1 ? ARGS[1] : error("pass the sweep directory")
t3dir = length(ARGS) >= 2 ? ARGS[2] : ""

kv(line) = Dict(m.captures[1] => m.captures[2]
                for m in eachmatch(r"(\w+)=(-?[\d.]+(?:[eE][-+]?\d+)?)", line))

function parse_log(path, system, T, arm)
    isfile(path) || return nothing
    txt = read(path, String)
    lines = split(txt, '\n')
    m   = match(r"solve complete: ([\d.]+) s, iterations=(\d+), status=(\d+)", txt)
    obj = match(r"FilterDDP objective=([-\d.eE+]+) max_equality_residual=([-\d.eE+]+)", txt)
    env = match(r"PIPELINE_ENV .*factor_backed=(\d)", txt)
    tl = [kv(l) for l in lines if startswith(l, "FILTERDDP_TIMING")]
    nl = [kv(l) for l in lines if startswith(l, "FILTERDDP_NNZ")]
    s(rows, k) = isempty(rows) ? NaN : sum(parse(Float64, r[k]) for r in rows)
    mean(rows, k) = isempty(rows) ? NaN : s(rows, k) / length(rows)
    finished = occursin("PIPELINE_STEP_OK", txt)
    status = m === nothing ? (finished ? "CRASHED" : "RUNNING") :
             (m.captures[3] == "0" ? "OK" : "FAIL" * m.captures[3])
    return (system = system, T = T, arm = arm,
            factor_backed = env === nothing ? 0 : parse(Int, env.captures[1]),
            status = status,
            wall = m === nothing ? NaN : parse(Float64, m.captures[1]),
            iters = m === nothing ? -1 : parse(Int, m.captures[2]),
            objective = obj === nothing ? NaN : parse(Float64, obj.captures[1]),
            max_eq = obj === nothing ? NaN : parse(Float64, obj.captures[2]),
            stages = length(tl),
            factor_s = s(tl, "factor_s"), solve_s = s(tl, "solve_s"),
            assembly_s = s(tl, "kkt_assembly_s"), deriv_s = s(tl, "derivative_s"),
            update_s = s(tl, "update_s"),
            per_stage_factor_ms = isempty(tl) ? NaN : 1000 * s(tl, "factor_s") / length(tl),
            nnz_K = mean(nl, "nnz_K"), nnz_LU = mean(nl, "nnz_LU"))
end

runs = NamedTuple[]
# T = 3 arms from the agenda pipeline (A_<sys>_T3_<exact|diag1e-8>.log)
if !isempty(t3dir) && isdir(t3dir)
    for sys in ("ieee123C_1ph", "ieee2522C_1ph", "large10kC_1ph"),
        (arm, tag) in (("exact", "exact"), ("diag", "diag1e-8"))
        r = parse_log(joinpath(t3dir, "A_$(sys)_T3_$(tag).log"), sys, 3, arm)
        r === nothing || push!(runs, r)
    end
end
# this sweep: <sys>_T<T>_<arm>.log
logdir = joinpath(sweep, "logs")
for f in (isdir(logdir) ? readdir(logdir) : String[])
    mm = match(r"^(\w+?)_T(\d+)_(exact|diag)\.log$", f)
    mm === nothing && continue
    r = parse_log(joinpath(logdir, f), mm.captures[1], parse(Int, mm.captures[2]), mm.captures[3])
    r === nothing || push!(runs, r)
end

sysorder = Dict("ieee123C_1ph" => 1, "ieee2522C_1ph" => 2, "large10kC_1ph" => 3)
sort!(runs, by = r -> (get(sysorder, r.system, 9), r.T, r.arm))

open(joinpath(sweep, "diag_hessian_horizon.csv"), "w") do io
    println(io, "system,T,arm,factor_backed,status,wall_s,iterations,objective,max_equality_residual,stage_solves,per_stage_factor_ms,factor_s,solve_s,assembly_s,derivative_s,update_s,mean_nnz_K,mean_nnz_LU")
    for r in runs
        @printf(io, "%s,%d,%s,%d,%s,%.3f,%d,%.12g,%.3e,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.1f,%.1f\n",
                r.system, r.T, r.arm, r.factor_backed, r.status, r.wall, r.iters,
                r.objective, r.max_eq, r.stages, r.per_stage_factor_ms, r.factor_s,
                r.solve_s, r.assembly_s, r.deriv_s, r.update_s, r.nnz_K, r.nnz_LU)
    end
end

println("=" ^ 118)
println("EXACT vs DIAGONAL STAGE HESSIAN ACROSS HORIZONS  (periodic profile, C_B = 1e-3, floor 1e-8)")
println("=" ^ 118)
@printf("%-14s %4s %3s | %10s %5s %7s | %10s %5s %7s | %7s %7s %9s %11s\n",
        "system", "T", "FB", "exact wall", "iters", "status", "diag wall", "iters", "status",
        "d wall", "d iters", "d factor", "obj rel")
println("-" ^ 118)
for key in unique([(r.system, r.T) for r in runs])
    e = findfirst(r -> (r.system, r.T, r.arm) == (key..., "exact"), runs)
    d = findfirst(r -> (r.system, r.T, r.arm) == (key..., "diag"), runs)
    E = e === nothing ? nothing : runs[e]
    D = d === nothing ? nothing : runs[d]
    col(x) = x === nothing ? ("--", "--", "--") :
        (isnan(x.wall) ? "--" : @sprintf("%.1f", x.wall), x.iters < 0 ? "--" : string(x.iters), x.status)
    ce, cd = col(E), col(D)
    both = E !== nothing && D !== nothing && E.status == "OK" && D.status == "OK"
    pct(a, b) = both ? @sprintf("%+.1f%%", 100 * (b / a - 1)) : "--"
    fb = E !== nothing ? E.factor_backed : (D !== nothing ? D.factor_backed : 0)
    @printf("%-14s %4d %3d | %10s %5s %7s | %10s %5s %7s | %7s %7s %9s %11s\n",
            key[1], key[2], fb, ce..., cd...,
            both ? pct(E.wall, D.wall) : "--",
            both ? pct(E.iters, D.iters) : "--",
            both ? pct(E.per_stage_factor_ms, D.per_stage_factor_ms) : "--",
            both ? @sprintf("%.2e", abs(D.objective - E.objective) / abs(E.objective)) : "--")
end
println("-" ^ 118)
println("FB = FILTERDDP_FACTOR_BACKED_POLICY for both arms. d = diag relative to exact.")
println("T = 3 rows come from the agenda pipeline logs; large10k T = 3 is the pipeline run")
println("(a second run gave -57% wall; see FILTERDDP_DIAGONAL_HESSIAN.md).")
