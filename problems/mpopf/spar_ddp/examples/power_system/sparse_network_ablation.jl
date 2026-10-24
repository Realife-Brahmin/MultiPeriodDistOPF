# Experimental ablation only. Default remains the radial tree backend.
# This backend uses UMFPACK on the network submatrix and constructs the same
# retained/halved battery system. It does NOT disable the temporal scan or
# form a dense inverse. Only the exact-curvature, half-block path is supported.
function sparse_network_kkt_factor(lay, K)
    _half_block() && _lean_exact() || error("umfpack network ablation needs HALF_BLOCK=1 and LEAN_EXACT=1")
    network = sort!(reduce(vcat, lay.own))
    return SparseNetworkKKTFactor(lu(K[network, network]), network, ones(1, 1))
end

function sparse_network_solve!(x, fac::SparseNetworkKKTFactor, b)
    x[fac.network] = fac.factor \ b[fac.network]
    return x
end

# Shared entry point used by the process/thread pass's policy recovery.
function _tp_network_solve!(x, lay, fac::SparseNetworkKKTFactor, ooff, pre, z, b)
    return sparse_network_solve!(x, fac, b)
end

function _tree_kkt_prepare(lay, stat::TreeKKTStatic, fac::SparseNetworkKKTFactor, K::SparseMatrixCSC)
    E, F = lay.E, fac.network
    S = Matrix(K[E, E])
    KEF = K[E, F]
    width = parse(Int, get(ENV, "FILTERDDP_SPARSE_SCHUR_BATCH", "16"))
    width > 0 || error("FILTERDDP_SPARSE_SCHUR_BATCH must be positive")
    # At most length(F) x width dense response entries at once. Never N^-1.
    for first in 1:width:length(E)
        cols = first:min(first + width - 1, length(E))
        response = fac.factor \ Matrix(K[F, E[cols]])
        @views S[:, cols] .-= KEF * response
    end
    # The same independent energy auxiliaries as on the tree path.
    if !isempty(lay.iso)
        D = K[lay.iso, lay.iso]
        isdiag(D) || error("energy auxiliary block is not diagonal")
        all(!iszero, diag(D)) || error("singular energy auxiliary block")
        S .-= Matrix(K[E, lay.iso] * Diagonal(1 ./ diag(D)) * K[lay.iso, E])
    end
    all(isfinite, S) || error("nonfinite sparse-network Schur complement")
    return S, zeros(length(E), 1)
end

function _tree_kkt_finish0(lay, stat::TreeKKTStatic, fac::SparseNetworkKKTFactor,
                           K::SparseMatrixCSC, S, Up, nB::Int)
    nE = length(lay.E)
    nE == 2nB || error("expected power/energy interface of size 2nB")
    energy = nB+1:nE
    NN = S[energy, energy]
    # Energy rows have no network coupling, so this is exact structural zero.
    isdiag(NN) || error("energy Schur block must be diagonal")
    all(!iszero, diag(NN)) || error("singular energy Schur block")
    # Preserve the existing grouped ordering, messages, and half-block solves.
    nP = Int[]; Ninv = Vector{Float64}[]
    PN = SparseMatrixCSC{Float64,Int}[]; NP = SparseMatrixCSC{Float64,Int}[]
    prow = Int[]
    for rg in stat.granges
        ids = stat.perm[rg]
        np = count(<=(nB), ids)
        all(k -> (ids[k] <= nB) == (k <= np), eachindex(ids)) || error("unexpected interface ordering")
        powers, energies = ids[1:np], ids[np+1:end]
        outside = setdiff(energy, energies)
        iszero(S[powers, outside]) && iszero(S[outside, powers]) || error("cross-group power/energy coupling")
        push!(nP, np)
        push!(Ninv, 1 ./ diag(S[energies, energies]))
        push!(PN, sparse(S[powers, energies])); push!(NP, sparse(S[energies, powers]))
        append!(prow, rg[1:np])
    end
    bat = stat.perm[prow]
    # Same second elimination as production: 2nB -> nB, in grouped power order.
    Sh0 = S[bat, bat] - S[bat, energy] * (Diagonal(1 ./ diag(NN)) * S[energy, bat])
    return TreeKKTBlock0(Sh0, bat, nP, Ninv, PN, NP, prow, Up)
end

function _tree_kkt_finish(lay, stat::TreeKKTStatic, fac::SparseNetworkKKTFactor,
                          K::SparseMatrixCSC, Bd, Up, dense_check::Bool, Vc)
    error("umfpack network ablation supports only the exact-curvature half-block path; no fallback allowed")
end

# Capture is opt-in and belongs to a separate diagnostic run, never a timed arm.
function capture_network_ablation(lay, K)
    directory = ENV["TREE_ABLATION_CAPTURE_DIR"]
    mkpath(directory)
    file = joinpath(directory, "network_$(getpid()).jls")
    SPARDDP.Serialization.serialize(file, (K=copy(K), lay=lay,
        julia_version=string(VERSION), backend=get(ENV,"FILTERDDP_NETWORK_BACKEND","tree")))
    println("TREE_ABLATION_CAPTURE file=", file)
end
