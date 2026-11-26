using LinearAlgebra
using Printf
using Serialization
using SparseArrays

input = length(ARGS) >= 1 ? ARGS[1] : error("pass a captured KKT .jls file")
repeats = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3
payload = deserialize(input)
K, rhs = payload.K, payload.rhs
ncols = size(rhs, 2)

requested = length(ARGS) >= 3 ? parse.(Int, split(ARGS[3], ',')) :
    [1, 4, 16, 32, 64, 128, 256, 512, ncols]
blocks = unique(sort(clamp.(requested, 1, ncols)))
last(blocks) == ncols || push!(blocks, ncols)

median_value(xs) = sort(xs)[cld(length(xs), 2)]
factor_s = @elapsed F = lu(K)

function solve_in_blocks!(F, rhs, block_width)
    solve_s = 0.0
    copy_s = 0.0
    allocated = 0
    solution_norm_sq = 0.0
    for first_col in 1:block_width:size(rhs, 2)
        last_col = min(first_col + block_width - 1, size(rhs, 2))
        width = last_col - first_col + 1
        work = Matrix{eltype(rhs)}(undef, size(rhs, 1), width)
        copy_s += @elapsed copyto!(work, @view rhs[:, first_col:last_col])
        stats = @timed ldiv!(F, work)
        solve_s += stats.time
        allocated += stats.bytes
        solution_norm_sq += sum(abs2, work)
    end
    return solve_s, copy_s, allocated, sqrt(solution_norm_sq)
end

# Compile the vector and matrix paths before collecting measurements.
solve_in_blocks!(F, @view(rhs[:, 1:1]), 1)
solve_in_blocks!(F, @view(rhs[:, 1:min(2, ncols)]), min(2, ncols))

rows = NamedTuple[]
for block_width in blocks
    solve_times = Float64[]
    copy_times = Float64[]
    allocations = Int[]
    norms = Float64[]
    for _ in 1:repeats
        GC.gc()
        solve_s, copy_s, bytes, solution_norm = solve_in_blocks!(F, rhs, block_width)
        push!(solve_times, solve_s)
        push!(copy_times, copy_s)
        push!(allocations, bytes)
        push!(norms, solution_norm)
    end
    push!(rows, (
        block_width=block_width,
        blocks=cld(ncols, block_width),
        solve_s=median_value(solve_times),
        copy_s=median_value(copy_times),
        solve_alloc_MiB=median_value(allocations) / 2.0^20,
        solution_norm=median_value(norms),
    ))
end

# Independent residual check for two representative columns.
sample_cols = unique([1, ncols])
sample_rhs = Matrix(@view rhs[:, sample_cols])
sample_solution = copy(sample_rhs)
ldiv!(F, sample_solution)
relative_residual = norm(K * sample_solution - sample_rhs) / norm(sample_rhs)

@printf("KKT_RHS_BLOCK_BENCH size=%s nnz=%d rhs=%s repeats=%d factor_s=%.9f sample_relative_residual=%.3e\n",
    string(size(K)), nnz(K), string(size(rhs)), repeats, factor_s, relative_residual)
for row in rows
    @printf("KKT_RHS_BLOCK_BENCH block_width=%d blocks=%d solve_s=%.9f copy_s=%.9f solve_alloc_MiB=%.3f solution_norm=%.12e\n",
        row.block_width, row.blocks, row.solve_s, row.copy_s,
        row.solve_alloc_MiB, row.solution_norm)
end

output = replace(input, r"\.jls$" => "_rhs_blocks.csv")
open(output, "w") do io
    println(io, "block_width,blocks,factor_s,solve_s,copy_s,solve_alloc_MiB,solution_norm,sample_relative_residual")
    for row in rows
        @printf(io, "%d,%d,%.9f,%.9f,%.9f,%.6f,%.12e,%.6e\n",
            row.block_width, row.blocks, factor_s, row.solve_s, row.copy_s,
            row.solve_alloc_MiB, row.solution_norm, relative_residual)
    end
end
println("wrote $output")
