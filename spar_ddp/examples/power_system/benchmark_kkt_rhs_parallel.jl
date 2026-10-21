using LinearAlgebra
using Printf
using Serialization
using SparseArrays

input = length(ARGS) >= 1 ? ARGS[1] : error("pass a captured KKT .jls file")
workers = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
repeats = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 2
factor_mode = length(ARGS) >= 4 ? ARGS[4] : "independent"
factor_mode in ("independent", "shared") || error("factor mode must be independent or shared")
workers <= Threads.nthreads() || error("workers=$workers exceeds Julia threads=$(Threads.nthreads())")

payload = deserialize(input)
K, rhs = payload.K, payload.rhs
payload = nothing
GC.gc()

function column_ranges(ncols, workers)
    ranges = UnitRange{Int}[]
    for worker in 1:workers
        first_col = fld((worker - 1) * ncols, workers) + 1
        last_col = fld(worker * ncols, workers)
        push!(ranges, first_col:last_col)
    end
    ranges
end

function one_run(K, rhs, workers, factor_mode)
    ranges = column_ranges(size(rhs, 2), workers)
    factors = Vector{Any}(undef, workers)
    factor_s = if factor_mode == "shared"
        local shared_factor
        elapsed = @elapsed shared_factor = lu(K)
        fill!(factors, shared_factor)
        elapsed
    else
        @elapsed Threads.@sync for worker in 1:workers
            Threads.@spawn factors[worker] = lu(K)
        end
    end

    blocks = Vector{Matrix{eltype(rhs)}}(undef, workers)
    for worker in 1:workers
        blocks[worker] = Matrix(@view rhs[:, ranges[worker]])
    end
    solve_s = @elapsed Threads.@sync for worker in 1:workers
        Threads.@spawn ldiv!(factors[worker], blocks[worker])
    end

    solution_norm = sqrt(sum(sum(abs2, block) for block in blocks))
    sample_relative_residuals = zeros(Float64, workers)
    Threads.@sync for worker in 1:workers
        Threads.@spawn begin
            col = first(ranges[worker])
            x = @view blocks[worker][:, 1]
            b = @view rhs[:, col]
            sample_relative_residuals[worker] = norm(K * x - b) / norm(b)
        end
    end
    max_sample_residual = maximum(sample_relative_residuals)
    return factor_s, solve_s, solution_norm, max_sample_residual
end

# Warm compilation on one column without retaining the temporary factor.
warm_factor = lu(K)
warm_block = Matrix(@view rhs[:, 1:1])
ldiv!(warm_factor, warm_block)
warm_factor = nothing
warm_block = nothing
GC.gc()

factor_times = Float64[]
solve_times = Float64[]
solution_norms = Float64[]
residuals = Float64[]
for _ in 1:repeats
    GC.gc()
    factor_s, solve_s, solution_norm, residual = one_run(K, rhs, workers, factor_mode)
    push!(factor_times, factor_s)
    push!(solve_times, solve_s)
    push!(solution_norms, solution_norm)
    push!(residuals, residual)
    GC.gc()
end

median_value(xs) = sort(xs)[cld(length(xs), 2)]
factor_s = median_value(factor_times)
solve_s = median_value(solve_times)
solution_norm = median_value(solution_norms)
max_sample_residual = maximum(residuals)
maxrss_MiB = Sys.maxrss() / 2.0^20

@printf("KKT_PARALLEL workers=%d factor_mode=%s julia_threads=%d repeats=%d factor_s=%.9f solve_s=%.9f total_s=%.9f solution_norm=%.12e max_sample_relative_residual=%.3e maxrss_MiB=%.3f\n",
    workers, factor_mode, Threads.nthreads(), repeats, factor_s, solve_s, factor_s + solve_s,
    solution_norm, max_sample_residual, maxrss_MiB)
