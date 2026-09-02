using LinearAlgebra
using MUMPS
using Printf
using Serialization
using SparseArrays

input = length(ARGS) >= 1 ? ARGS[1] : error("pass a captured KKT .jls file")
repeats = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 7
payload = deserialize(input)
K, rhs = payload.K, payload.rhs
MUMPS.MPI.Initialized() || MUMPS.MPI.Init()

middle(xs) = sort(xs)[cld(length(xs), 2)]

function bench_umfpack(K, rhs, repeats)
    F = lu(K); X = F \ rhs # compile/warm
    factor_times = Float64[]; solve_times = Float64[]
    for _ in 1:repeats
        GC.gc()
        local Ft
        push!(factor_times, @elapsed Ft = lu(K))
        local Xt
        push!(solve_times, @elapsed Xt = Ft \ rhs)
        X = Xt
    end
    return middle(factor_times), middle(solve_times), norm(K*X-rhs)/norm(rhs)
end

function bench_mumps(K, rhs, repeats; sym=MUMPS.mumps_unsymmetric)
    make_factor() = begin
        F = MUMPS.Mumps{eltype(K)}(sym)
        MUMPS.associate_matrix!(F, K)
        MUMPS.suppress_display!(F)
        MUMPS.mumps_factorize!(F)
        F
    end
    F = make_factor(); X = F \ rhs # compile/warm
    finalize(F)
    factor_times = Float64[]; solve_times = Float64[]
    for _ in 1:repeats
        GC.gc()
        local Ft
        push!(factor_times, @elapsed Ft = make_factor())
        local Xt
        push!(solve_times, @elapsed Xt = Ft \ rhs)
        X = Xt
        finalize(Ft)
    end
    return middle(factor_times), middle(solve_times), norm(K*X-rhs)/norm(rhs)
end

results = NamedTuple[]
uf, us, ur = bench_umfpack(K, rhs, repeats)
push!(results, (solver="UMFPACK", factor_s=uf, solve_s=us, residual=ur))
mf, ms, mr = bench_mumps(K, rhs, repeats)
push!(results, (solver="MUMPS_unsymmetric", factor_s=mf, solve_s=ms, residual=mr))
try
    sf, ss, sr = bench_mumps(K, rhs, repeats; sym=MUMPS.mumps_symmetric)
    push!(results, (solver="MUMPS_symmetric", factor_s=sf, solve_s=ss, residual=sr))
catch err
    @warn "symmetric MUMPS benchmark failed" exception=(err, catch_backtrace())
end

@printf("KKT_BENCH size=%s nnz=%d rhs=%s repeats=%d\n",
    string(size(K)), nnz(K), string(size(rhs)), repeats)
for r in results
    @printf("KKT_BENCH solver=%s factor_s=%.9f solve_s=%.9f total_s=%.9f relative_residual=%.3e\n",
        r.solver, r.factor_s, r.solve_s, r.factor_s+r.solve_s, r.residual)
end

output = replace(input, r"\.jls$" => "_benchmark.csv")
open(output, "w") do io
    println(io, "solver,factor_s,solve_s,total_s,relative_residual")
    for r in results
        @printf(io, "%s,%.9f,%.9f,%.9f,%.6e\n",
            r.solver, r.factor_s, r.solve_s, r.factor_s+r.solve_s, r.residual)
    end
end
println("wrote $output")
MUMPS.MPI.Finalized() || MUMPS.MPI.Finalize()
