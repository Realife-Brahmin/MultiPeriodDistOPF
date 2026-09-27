# Binned sparsity patterns for the ordering figures: the captured K, the matrix
# after each configuration's ordering, and the resulting factor.
#
# UMFPACK configs: reordered = K[p,q] from the factorization itself; factor =
#   the actual L and U (both in pivot order).
# MUMPS configs: the pivot order comes from MUMPS's own analysis (SYM_PERM);
#   reordered = K[o,o]; factor = the exact Cholesky pattern of that ordering
#   on an SPD surrogate of |K| (MUMPS keeps its factors internal). This is the
#   pattern before delayed pivots, i.e. what the ordering itself implies.
#
# Patterns are binned onto an nb x nb grid of nonzero counts (large10k has
# ~1e5 rows), with per-matrix statistics, for plot_kkt_ordering_patterns.py.
#
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/kkt_ordering_patterns.jl \
#         <capture.jls> <tag> <outdir> <config>[,<config>...] [bins=600]

using LinearAlgebra, SparseArrays, Serialization, Printf, MUMPS

const U = SparseArrays.UMFPACK
const LSS = SparseArrays.LibSuiteSparse

input, tag, outdir = ARGS[1], ARGS[2], ARGS[3]
configs = split(ARGS[4], ',')
nb = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 600
mkpath(outdir)

cap = deserialize(input)
K = SparseMatrixCSC{Float64, Int64}(cap.K)
n = size(K, 1)
nb = min(nb, n)

# Orderings push the dense battery-curvature block and its coupling rows to the
# end, which is where fill differs; a second, finer grid covers the trailing
# `tail` rows and columns.
const tail = min(n, 4000)
const nbt = min(300, tail)
function write_grid(path, I, J, lo, span, bins, header)
    C = zeros(Int, bins, bins)
    for k in eachindex(I)
        (I[k] > lo && J[k] > lo) || continue
        C[cld((I[k] - lo) * bins, span), cld((J[k] - lo) * bins, span)] += 1
    end
    open(path, "w") do io
        println(io, header)
        for i in 1:bins
            println(io, join(C[i, :], ","))
        end
    end
end
function write_bins(name, S)
    I, J, _ = findnz(S)
    bw = isempty(I) ? 0 : maximum(abs.(I .- J))
    path = joinpath(outdir, "$(tag)__$(name).csv")
    write_grid(path, I, J, 0, n, nb, @sprintf("# n=%d nnz=%d bins=%d bandwidth=%d", n, nnz(S), nb, bw))
    lo = n - tail
    ntail = count(k -> I[k] > lo && J[k] > lo, eachindex(I))
    write_grid(joinpath(outdir, "$(tag)__$(name)__tail.csv"), I, J, lo, tail, nbt,
               @sprintf("# n=%d nnz=%d bins=%d bandwidth=%d tail=%d", n, ntail, nbt, bw, tail))
    @printf("  %-40s nnz=%10d  bandwidth=%6d  trailing-%d nnz=%d\n", name, nnz(S), bw, tail, ntail)
end

write_bins("K", K)

umf = Dict(
    "UMFPACK_default"      => (nothing, nothing),
    "UMFPACK_unsym_colamd" => (LSS.UMFPACK_STRATEGY_UNSYMMETRIC, LSS.UMFPACK_ORDERING_AMD),
    "UMFPACK_sym_amd"      => (LSS.UMFPACK_STRATEGY_SYMMETRIC, LSS.UMFPACK_ORDERING_AMD),
    "UMFPACK_unsym_metis"  => (LSS.UMFPACK_STRATEGY_UNSYMMETRIC, LSS.UMFPACK_ORDERING_METIS),
    "UMFPACK_sym_metis"    => (LSS.UMFPACK_STRATEGY_SYMMETRIC, LSS.UMFPACK_ORDERING_METIS),
    "UMFPACK_cholmod"      => (nothing, LSS.UMFPACK_ORDERING_CHOLMOD),
    "UMFPACK_best"         => (nothing, LSS.UMFPACK_ORDERING_BEST))
mumps_icntl = Dict(
    "MUMPS_sym_auto" => [7 => 7], "MUMPS_sym_amd" => [7 => 0], "MUMPS_sym_amf" => [7 => 2],
    "MUMPS_sym_qamd" => [7 => 6], "MUMPS_sym_pord" => [7 => 4], "MUMPS_sym_scotch" => [7 => 3],
    "MUMPS_sym_metis" => [7 => 5], "MUMPS_sym_metis_compressed" => [7 => 5, 12 => 2],
    "MUMPS_sym_ipopt" => [7 => 7, 8 => 77, 6 => 7], "MUMPS_sym_amd_plain" => [7 => 0, 12 => 1],
    "MUMPS_sym_amf_plain" => [7 => 2, 12 => 1], "MUMPS_sym_metis_plain" => [7 => 5, 12 => 1])

const Ssurr = let A = abs.(K); A + spdiagm(0 => vec(sum(A; dims=2)) .+ 1.0) end
MUMPS.MPI.Initialized() || MUMPS.MPI.Init()

for c in configs
    if haskey(umf, c)
        strat, ord = umf[c]
        ctrl = U.get_umfpack_control(Float64, Int64)
        strat === nothing || (ctrl[LSS.UMFPACK_STRATEGY + 1] = strat)
        ord === nothing || (ctrl[LSS.UMFPACK_ORDERING + 1] = ord)
        F = lu(K; control=ctrl)
        write_bins("$(c)__reordered", K[F.p, F.q])
        write_bins("$(c)__factor", F.L + F.U)
    elseif haskey(mumps_icntl, c)
        m = MUMPS.Mumps{Float64}(MUMPS.mumps_symmetric, MUMPS.get_icntl(), MUMPS.default_cntl64)
        MUMPS.suppress_display!(m)
        for (i, v) in mumps_icntl[c]; MUMPS.set_icntl!(m, i, v; displaylevel=0); end
        MUMPS.associate_matrix!(m, K)
        MUMPS.set_job!(m, MUMPS.ANALYZE); MUMPS.invoke_mumps!(m)
        m.infog[1] < 0 && error("MUMPS analysis failed: INFOG(1)=$(m.infog[1])")
        order = invperm(copy(unsafe_wrap(Array, m.sym_perm, n)))
        finalize(m)
        write_bins("$(c)__reordered", K[order, order])
        write_bins("$(c)__factor", sparse(ldlt(Symmetric(Ssurr[order, order]); perm=collect(1:n)).LD))
    else
        error("unknown config $c")
    end
end
MUMPS.MPI.Finalized() || MUMPS.MPI.Finalize()
