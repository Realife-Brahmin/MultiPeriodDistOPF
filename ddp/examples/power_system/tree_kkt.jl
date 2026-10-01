# Block-tree elimination of a FilterDDP stage KKT on a radial network: the
# hand-written counterpart of MUMPS's Schur complement onto the battery rows
# (ddp/notes/PARALLEL_IN_TIME.md, Section 4).
#
# The stage KKT K = [H cu'; cu 0] of ieee123c_filterddp.jl groups naturally by
# bus: a non-root bus j owns its incoming line's P, Q, ell, SOC slack, its
# voltage v_j, its DER control, and its own P-balance, Q-balance, voltage-drop
# and SOC rows. Those groups couple only to the parent bus and the children
# (checked against K's pattern when the layout is built), so eliminating the
# groups from the leaves to the substation is block Gaussian elimination on
# the feeder tree: every elimination is a small dense LU (~10 unknowns) and
# passes a 3x3 update to the parent. What is kept is E = {battery powers,
# energy rows}; the energy slacks couple only to E and are eliminated exactly
# on their own. The Schur complement onto E is
#     S = K_EE - K_E,Fc * Ninv[Fc, Fc] * K_Fc,E - K_E,iso * inv(K_iso) * K_iso,E,
# where N is the network part, Fc the network rows coupled to E (the
# P-balance rows of battery buses), and Ninv[Fc, Fc] is obtained by one sweep
# up the battery paths and one sweep down the union of those paths.

using LinearAlgebra, SparseArrays

struct TreeKKTLayout
    n::Int
    root::Int
    post::Vector{Int}                  # nodes, children before parent, root last
    parent::Vector{Int}                # 0 for the root
    children::Vector{Vector{Int}}
    own::Vector{Vector{Int}}           # K indices eliminated at each node
    iface::Vector{Vector{Int}}         # K indices of the parent coupled to own
    iface_pos::Vector{Vector{Int}}     # their positions in the parent's own
    E::Vector{Int}                     # kept: battery powers, energy rows
    iso::Vector{Int}                   # energy slacks (coupled only to E)
    Fc::Vector{Int}                    # network indices coupled to E
    own_Fc::Vector{Vector{Int}}        # per node: positions in Fc of its own Fc indices
    own_Fc_pos::Vector{Vector{Int}}    # ... and their positions in own
    cols::Vector{Vector{Int}}          # per node: Fc positions in its subtree, block order
    # nzval positions for the dense local blocks (pattern is fixed)
    oo::Vector{Vector{NTuple{3,Int}}}  # (row, col, nz) within own x own
    oi::Vector{Vector{NTuple{3,Int}}}  # own x iface
end

# K index -> owner: node id (> 0), -1 for E, -2 for iso.
function tree_kkt_layout(data, idx, nu::Int, K::SparseMatrixCSC)
    buses, lines, nonroot = data[:Nset], data[:Lset], data[:Nm1set]
    rootbus = data[:substationBus]
    N, L, nB = length(buses), length(lines), length(data[:Bset])
    buspos = Dict(j => k for (k, j) in enumerate(buses))
    linepos = Dict(e => k for (k, e) in enumerate(lines))
    derpos = Dict(j => k for (k, j) in enumerate(data[:Dset]))
    n = size(K, 1)
    root = buspos[rootbus]
    own = [Int[] for _ in 1:N]
    # Rows follow `equations` in ieee123c_filterddp.jl: P balance (root, then
    # nonroot), Q balance (root, nonroot), voltage drop and SOC per line, root
    # voltage, energy per battery.
    append!(own[root], [idx.ps, idx.qs, idx.v[root], nu + 1, nu + N + 1, nu + 2N + 2L + 1])
    for (k, j) in enumerate(nonroot)
        b = buspos[j]
        e = linepos[(data[:parent][j], j)]
        append!(own[b], [idx.P[e], idx.Q[e], idx.ell[e], idx.soc_slack[e], idx.v[b],
                         nu + 1 + k, nu + N + 1 + k, nu + 2N + e, nu + 2N + L + e])
    end
    for (d, j) in enumerate(data[:Dset])
        push!(own[buspos[j]], idx.qnorm[d])
    end
    E = vcat(collect(idx.pb), [nu + 2N + 2L + 1 + b for b in 1:nB])
    iso = collect(idx.energy_slack)

    owner = zeros(Int, n)
    for j in 1:N, i in own[j]
        owner[i] == 0 || error("index $i assigned twice")
        owner[i] = j
    end
    for i in E; owner[i] == 0 || error("E index $i assigned twice"); owner[i] = -1; end
    for i in iso; owner[i] == 0 || error("iso index $i assigned twice"); owner[i] = -2; end
    all(!=(0), owner) || error("$(count(==(0), owner)) K indices not assigned to any block")

    parent = zeros(Int, N)
    children = [Int[] for _ in 1:N]
    for j in nonroot
        b, p = buspos[j], buspos[data[:parent][j]]
        parent[b] = p
        push!(children[p], b)
    end
    # Post-order without recursion.
    post = Int[]; stack = [(root, false)]
    while !isempty(stack)
        j, done = pop!(stack)
        if done
            push!(post, j)
        else
            push!(stack, (j, true))
            for c in children[j]; push!(stack, (c, false)); end
        end
    end

    # Interfaces, and a check that the pattern really is a block tree.
    rows = rowvals(K)
    iface = [Int[] for _ in 1:N]
    for j in 1:N, c in own[j], p in nzrange(K, c)
        i = rows[p]; o = owner[i]
        if o > 0 && o != j
            if o == parent[j]
                push!(iface[j], i)
            elseif parent[o] != j
                error("K couples node $j to node $o, which is neither parent nor child")
            end
        elseif o == -2
            error("energy slack $i couples to the network")
        end
    end
    for i in iso, p in nzrange(K, i)
        owner[rows[p]] in (-1, -2) || error("energy slack $i couples outside E")
    end
    iface = [sort!(unique(v)) for v in iface]
    ownpos = [Dict(i => t for (t, i) in enumerate(own[j])) for j in 1:N]
    iface_pos = [j == root ? Int[] : [ownpos[parent[j]][i] for i in iface[j]] for j in 1:N]

    # Network indices coupled to E.
    isE = falses(n); isE[E] .= true
    Fc = Int[]
    for c in E, p in nzrange(K, c)
        i = rows[p]
        owner[i] > 0 && push!(Fc, i)
    end
    Fc = sort!(unique(Fc))
    Fcpos = Dict(i => t for (t, i) in enumerate(Fc))
    own_Fc = [Int[] for _ in 1:N]; own_Fc_pos = [Int[] for _ in 1:N]
    for i in Fc
        j = owner[i]
        push!(own_Fc[j], Fcpos[i]); push!(own_Fc_pos[j], ownpos[j][i])
    end
    cols = [Int[] for _ in 1:N]
    for j in post
        cols[j] = vcat(own_Fc[j], (cols[c] for c in children[j])...)
    end

    # nzval positions of the local blocks.
    oo = [NTuple{3,Int}[] for _ in 1:N]
    oi = [NTuple{3,Int}[] for _ in 1:N]
    for j in 1:N
        ipos = Dict(i => t for (t, i) in enumerate(iface[j]))
        for (tc, c) in enumerate(own[j]), p in nzrange(K, c)
            i = rows[p]
            if owner[i] == j
                push!(oo[j], (ownpos[j][i], tc, p))
            elseif haskey(ipos, i)
                push!(oi[j], (tc, ipos[i], p))      # K[i, c], i in iface: stored as (own col, iface row)
            end
        end
    end
    return TreeKKTLayout(n, root, post, parent, children, own, iface, iface_pos,
                         E, iso, Fc, own_Fc, own_Fc_pos, cols, oo, oi)
end

struct TreeKKTFactor
    F::Vector{Any}                    # LU of each node's own block
    W::Vector{Matrix{Float64}}        # A_oo \ A_oi
    Aio::Vector{Matrix{Float64}}      # A_io (iface x own)
end

# Block LU from the leaves to the root. K is symmetric, so A_oi = A_io'.
function tree_kkt_factor(lay::TreeKKTLayout, K::SparseMatrixCSC)
    N = length(lay.own)
    nz = nonzeros(K)
    F = Vector{Any}(undef, N)
    W = Vector{Matrix{Float64}}(undef, N)
    Aio = Vector{Matrix{Float64}}(undef, N)
    C = Vector{Matrix{Float64}}(undef, N)
    for j in lay.post
        m = length(lay.own[j]); q = length(lay.iface[j])
        A = zeros(m, m)
        for (r, c, p) in lay.oo[j]; A[r, c] = nz[p]; end
        for c in lay.children[j]
            pos = lay.iface_pos[c]
            @views A[pos, pos] .+= C[c]
        end
        Fj = lu!(A)
        F[j] = Fj
        if j != lay.root
            Bio = zeros(q, m)
            for (c, r, p) in lay.oi[j]; Bio[r, c] = nz[p]; end
            Wj = Fj \ Matrix(Bio')
            W[j] = Wj; Aio[j] = Bio
            C[j] = -(Bio * Wj)
        end
    end
    return TreeKKTFactor(F, W, Aio)
end

# Ninv[Fc, Fc]: sweep up the battery paths, then down their union.
function tree_kkt_fc_block(lay::TreeKKTLayout, fac::TreeKKTFactor)
    N = length(lay.own); nF = length(lay.Fc)
    Z = Vector{Matrix{Float64}}(undef, N)
    msg = Vector{Matrix{Float64}}(undef, N)
    for j in lay.post
        k = length(lay.cols[j]); k == 0 && continue
        m = length(lay.own[j])
        R = zeros(m, k)
        for (t, pos) in enumerate(lay.own_Fc_pos[j]); R[pos, t] = 1.0; end
        off = length(lay.own_Fc[j])
        for c in lay.children[j]
            kc = length(lay.cols[c]); kc == 0 && continue
            @views R[lay.iface_pos[c], off+1:off+kc] .+= msg[c]
            off += kc
        end
        Zj = fac.F[j] \ R
        Z[j] = Zj
        j != lay.root && (msg[j] = -(fac.Aio[j] * Zj))
    end
    Y = zeros(nF, nF)
    yroot = zeros(length(lay.own[lay.root]), nF)
    yroot[:, lay.cols[lay.root]] .= Z[lay.root]
    stack = [(lay.root, yroot)]
    while !isempty(stack)
        j, yj = pop!(stack)
        for (t, pos) in zip(lay.own_Fc[j], lay.own_Fc_pos[j])
            @views Y[t, :] .= yj[pos, :]
        end
        for c in lay.children[j]
            isempty(lay.cols[c]) && continue
            yI = yj[lay.iface_pos[c], :]
            yc = -(fac.W[c] * yI)
            @views yc[:, lay.cols[c]] .+= Z[c]
            push!(stack, (c, yc))
        end
    end
    return Y
end

# Schur complement of K onto E, in the order of lay.E.
function tree_kkt_schur(lay::TreeKKTLayout, K::SparseMatrixCSC)
    fac = tree_kkt_factor(lay, K)
    Y = tree_kkt_fc_block(lay, fac)
    KEF = K[lay.E, lay.Fc]                      # sparse: one entry per battery
    S = Matrix(K[lay.E, lay.E])
    S .-= (KEF * Y) * KEF'
    # Each energy slack couples to its own energy row only.
    Epos = Dict(i => t for (t, i) in enumerate(lay.E))
    rows, vals = rowvals(K), nonzeros(K)
    for i in lay.iso
        d = K[i, i]
        nzE = [(Epos[rows[p]], vals[p]) for p in nzrange(K, i) if haskey(Epos, rows[p])]
        for (a, va) in nzE, (b, vb) in nzE
            S[a, b] -= va * vb / d
        end
    end
    return S
end
