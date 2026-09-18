# Turns run_agenda_pipeline.sh's raw logs and CSVs into one readable summary.
# Pure parsing -- no solver runs, no judgement calls. Safe to re-run any time.
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/summarize_agenda_pipeline.jl <outdir>

using Printf

out = length(ARGS) >= 1 ? ARGS[1] : error("pass the pipeline output directory")
logdir = joinpath(out, "logs")
isdir(logdir) || error("no logs directory under $out")

kv(line) = Dict(m.captures[1] => m.captures[2]
                for m in eachmatch(r"(\w+)=(-?[\d.]+(?:[eE][-+]?\d+)?)", line))

"""Pull every number the summary needs out of one solver log."""
function parse_run(path)
    isfile(path) || return nothing
    txt = read(path, String)
    m = match(r"solve complete: ([\d.]+) s, iterations=(\d+), status=(\d+)", txt)
    m === nothing && return nothing
    obj = match(r"FilterDDP objective=([-\d.eE+]+)", txt)
    tl = [kv(l) for l in eachline(IOBuffer(txt)) if startswith(l, "FILTERDDP_TIMING")]
    nl = [kv(l) for l in eachline(IOBuffer(txt)) if startswith(l, "FILTERDDP_NNZ")]
    s(rows, k) = isempty(rows) ? 0.0 : sum(parse(Float64, r[k]) for r in rows)
    mean(rows, k) = isempty(rows) ? 0.0 : s(rows, k) / length(rows)
    return (wall = parse(Float64, m.captures[1]),
            iters = parse(Int, m.captures[2]),
            status = parse(Int, m.captures[3]),
            objective = obj === nothing ? NaN : parse(Float64, obj.captures[1]),
            stages = length(tl),
            factor_s = s(tl, "factor_s"), solve_s = s(tl, "solve_s"),
            assembly_s = s(tl, "kkt_assembly_s"), deriv_s = s(tl, "derivative_s"),
            update_s = s(tl, "update_s"),
            per_stage_factor_ms = isempty(tl) ? NaN : 1000 * s(tl, "factor_s") / length(tl),
            nnz_K = mean(nl, "nnz_K"), nnz_LU = mean(nl, "nnz_LU"))
end

println("=" ^ 100)
println("AGENDA PIPELINE SUMMARY -- ", out)
println("=" ^ 100)

# ------------------------------------------------------------- A: Hessian ----
println("\nA. DIAGONAL STAGE HESSIAN  (converge? time benefit?)\n")
@printf("%-16s %-9s %8s %6s %7s %12s %14s %11s %11s %7s\n",
        "system", "arm", "wall_s", "iters", "status", "objective",
        "per-stage f/ms", "factor_s", "solve_s", "fill")
arms = ["exact", "diag1e-8", "diag1e-4", "diag1e-1"]
systems = ["ieee123C_1ph", "ieee2522C_1ph", "large10kC_1ph"]
base = Dict{String,Any}()
for sys in systems, arm in arms
    r = parse_run(joinpath(logdir, "A_$(sys)_T3_$(arm).log"))
    r === nothing && continue
    arm == "exact" && (base[sys] = r)
    fill = r.nnz_K > 0 ? r.nnz_LU / r.nnz_K : NaN
    @printf("%-16s %-9s %8.2f %6d %7s %12.6f %14.2f %11.2f %11.2f %7.3f\n",
            sys, arm, r.wall, r.iters, r.status == 0 ? "OK" : "FAIL$(r.status)",
            r.objective, r.per_stage_factor_ms, r.factor_s, r.solve_s, fill)
end
println("\n   deltas against each system's exact arm:")
for sys in systems
    haskey(base, sys) || continue
    b = base[sys]
    r = parse_run(joinpath(logdir, "A_$(sys)_T3_diag1e-8.log"))
    r === nothing && continue
    @printf("   %-16s wall %+6.1f%%   iters %+6.1f%%   per-stage factor %+6.1f%%   nnz_LU %+6.1f%%   obj rel diff %.2e\n",
            sys, 100*(r.wall/b.wall - 1), 100*(r.iters/b.iters - 1),
            100*(r.per_stage_factor_ms/b.per_stage_factor_ms - 1),
            100*(r.nnz_LU/b.nnz_LU - 1),
            abs(r.objective - b.objective) / max(abs(b.objective), eps()))
end

# ------------------------------------------------- B: stale-factor reuse -----
println("\n\nB. STALE FACTOR AS PRECONDITIONER  (does A^-1 vary slowly?)")
println("   break-even: refinement must converge in <= 1 + factor/solve steps (~3)\n")
for sys in systems
    csv = joinpath(out, "B_stale_factor_$(sys).csv")
    isfile(csv) || continue
    rows = [split(l, ',') for l in eachline(csv)][2:end]
    isempty(rows) && continue
    println("   $sys:")
    @printf("   %7s %5s %14s %14s %12s %9s %8s\n",
            "anchor", "lag", "dK/K", "dSigma/Sigma", "rho", "refine", "gmres")
    for r in rows
        lag = parse(Int, r[3])
        lag <= 3 || continue
        @printf("   %7s %5d %14.3e %14.3e %12.3e %9s %8s\n",
                r[1], lag, parse(Float64, r[4]), parse(Float64, r[5]),
                parse(Float64, r[6]), r[7], r[9])
    end
    best = minimum(parse(Int, r[7]) for r in rows if parse(Int, r[3]) >= 1 && parse(Int, r[7]) > 0)
    @printf("   -> best refinement step count at any lag >= 1: %d  (%s)\n\n",
            best, best <= 3 ? "PAYS" : "does not pay")
end

# ------------------------------------------------------- C: linear solvers ---
println("\nC. LINEAR SOLVERS  (UMFPACK vs MUMPS, threads, pivot tolerance)\n")
for sys in systems
    header_done = false
    for nt in (1, 2, 4, 8, 16)
        csv = joinpath(out, "C_solvers_$(sys)_threads$(nt).csv")
        isfile(csv) || continue
        if !header_done
            println("   $sys:")
            @printf("   %-22s %7s %11s %11s %11s %12s %11s\n",
                    "solver", "threads", "factor_s", "solve_s", "total_s", "nnz_LU", "residual")
            header_done = true
        end
        for l in Iterators.drop(eachline(csv), 1)
            r = split(l, ',')
            length(r) < 9 && continue
            @printf("   %-22s %7d %11.5f %11.5f %11.5f %12s %11.2e\n",
                    r[1], nt, parse(Float64, r[5]), parse(Float64, r[6]),
                    parse(Float64, r[7]), r[8] == "-1" ? "n/a" : r[8],
                    parse(Float64, r[9]))
        end
    end
    header_done && println()
end

println("=" ^ 100)
