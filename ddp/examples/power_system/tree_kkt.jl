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

# Dense LU with partial pivoting for the ~10x10 node blocks, in plain Julia:
# LAPACK/OpenBLAS calls on blocks this small cost far more in call and
# threading overhead than in arithmetic (~68 us per bus at large10k).
struct SmallLU
    A::Matrix{Float64}                # unit L below the diagonal, U on and above
    perm::Vector{Int}                 # row i of the factored matrix = row perm[i] of the original
end

function smalllu!(A::Matrix{Float64})
    n = size(A, 1); perm = collect(1:n)
    @inbounds for k in 1:n
        p = k; amax = abs(A[k, k])
        for i in k+1:n
            a = abs(A[i, k]); a > amax && (amax = a; p = i)
        end
        amax == 0 && throw(SingularException(k))
        if p != k
            for jj in 1:n; A[k, jj], A[p, jj] = A[p, jj], A[k, jj]; end
            perm[k], perm[p] = perm[p], perm[k]
        end
        inv_d = 1 / A[k, k]
        for i in k+1:n
            l = A[i, k] * inv_d
            A[i, k] = l
            l == 0 && continue
            for jj in k+1:n; A[i, jj] -= l * A[k, jj]; end
        end
    end
    return SmallLU(A, perm)
end

function Base.:\(F::SmallLU, B::AbstractMatrix)
    A = F.A; n = size(A, 1); w = size(B, 2)
    X = Matrix{Float64}(undef, n, w)
    @inbounds for c in 1:w
        for i in 1:n; X[i, c] = B[F.perm[i], c]; end
        for i in 2:n
            s = X[i, c]
            for k in 1:i-1; s -= A[i, k] * X[k, c]; end
            X[i, c] = s
        end
        for i in n:-1:1
            s = X[i, c]
            for k in i+1:n; s -= A[i, k] * X[k, c]; end
            X[i, c] = s / A[i, i]
        end
    end
    return X
end

# Small dense product without BLAS, for the same reason.
function smallmul(A::AbstractMatrix, B::AbstractMatrix)
    m, k = size(A); w = size(B, 2)
    C = zeros(m, w)
    @inbounds for c in 1:w, kk in 1:k
        b = B[kk, c]
        b == 0 && continue
        for i in 1:m; C[i, c] += A[i, kk] * b; end
    end
    return C
end

struct TreeKKTFactor
    F::Vector{SmallLU}                # LU of each node's own block
    W::Vector{Matrix{Float64}}        # A_oo \ A_oi
    Aio::Vector{Matrix{Float64}}      # A_io (iface x own)
    Mroot::Matrix{Float64}            # the root block with all contributions, unfactored
end

# Block LU from the leaves to the root. K is symmetric, so A_oi = A_io'.
function tree_kkt_factor(lay::TreeKKTLayout, K::SparseMatrixCSC)
    N = length(lay.own)
    nz = nonzeros(K)
    F = Vector{SmallLU}(undef, N)
    W = Vector{Matrix{Float64}}(undef, N)
    Aio = Vector{Matrix{Float64}}(undef, N)
    C = Vector{Matrix{Float64}}(undef, N)
    Mroot = zeros(0, 0)
    for j in lay.post
        m = length(lay.own[j]); q = length(lay.iface[j])
        A = zeros(m, m)
        for (r, c, p) in lay.oo[j]; A[r, c] = nz[p]; end
        for c in lay.children[j]
            pos = lay.iface_pos[c]; Cc = C[c]
            @inbounds for b in eachindex(pos), a in eachindex(pos)
                A[pos[a], pos[b]] += Cc[a, b]
            end
        end
        j == lay.root && (Mroot = copy(A))
        Fj = smalllu!(A)
        F[j] = Fj
        if j != lay.root
            Bio = zeros(q, m)
            for (c, r, p) in lay.oi[j]; Bio[r, c] = nz[p]; end
            Wj = Fj \ Bio'
            W[j] = Wj; Aio[j] = Bio
            Cj = smallmul(Bio, Wj); Cj .*= -1
            C[j] = Cj
        end
    end
    return TreeKKTFactor(F, W, Aio, Mroot)
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
        j != lay.root && (msg[j] = smallmul(fac.Aio[j], Zj) .* -1)
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
            yc = smallmul(fac.W[c], yI); yc .*= -1
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

# ---------------------------------------------------------------- structured --
# The substation (root) is the only place feeders meet, through the 3 root
# unknowns each feeder's top bus couples to. With N the network part, block
# elimination with the root last gives, exactly,
#     Ninv[Fc, Fc] = Bd + Psi * inv(M) * Psi',
# Bd = blockdiag over feeders of the feeder-only inverse, M the root block
# with all feeder contributions, and Psi the feeder rows' response to the
# root unknowns (-1 on root-owned rows). Under the diagonal Hessian K_EE
# couples a battery only to its own energy row, so the Schur complement onto
# E is block diagonal by feeder plus a rank-(root size) term:
#     S = D - U * inv(M) * U',   U = K_E,Fc * Psi,
# and S \ R follows from Woodbury with one small solve per feeder.
#
# Everything that depends only on the sparsity pattern (groups, positions in
# nonzeros(K)) is computed once in TreeKKTStatic.

struct TreeKKTStatic
    feeder::Vector{Int}                       # node -> feeder top (0 for the root)
    gkeys::Vector{Int}                        # group key: feeder top, or -1 for root-coupled E
    perm::Vector{Int}                         # E positions, grouped (permuted order)
    granges::Vector{UnitRange{Int}}           # each group's range in the permuted order
    kee::Vector{Vector{NTuple{3,Int}}}        # per group: (local row, local col, nz) of K_EE
    kef::Vector{NTuple{3,Int}}                # (E pos, Fc pos, nz) of K_E,Fc
    gkef::Vector{Vector{NTuple{3,Int}}}       # per group: (local row, feeder-local Fc col, nz)
    iso_d::Vector{Int}                        # nz of each energy slack's diagonal
    iso_e::Vector{Vector{NTuple{2,Int}}}      # per slack: (E pos, nz)
    giso::Vector{Vector{NTuple{5,Int}}}       # per group: (local a, local b, nz a, nz b, slack)
    colpos::Vector{Dict{Int,Int}}             # per feeder top: Fc pos -> local column
    ooff::Vector{Int}                         # node -> offset of its own block in a flat buffer
    nN::Int
    # Downward sweep: only the rows a node's children and its own Fc entries
    # read are computed (the bus's two balance rows and its voltage).
    need::Vector{Vector{Int}}                 # per node: positions in own
    iface_in_need::Vector{Vector{Int}}        # per node: its iface rows within need[parent]
    fc_in_need::Vector{Vector{Int}}           # per node: its own Fc entries within need
    cidx::Vector{Vector{Int}}                 # per node: feeder-local columns of cols[node]
    depth::Vector{Int}                        # depth below the feeder top (top = 1)
end

function _feeder_of(lay)
    N = length(lay.own)
    feeder = zeros(Int, N)
    stack = [(c, c) for c in lay.children[lay.root]]
    while !isempty(stack)
        j, f = pop!(stack)
        feeder[j] = f
        for c in lay.children[j]; push!(stack, (c, f)); end
    end
    return feeder
end

function tree_kkt_static(lay, K::SparseMatrixCSC)
    N = length(lay.own); nF = length(lay.Fc); nE = length(lay.E)
    feeder = _feeder_of(lay)
    owner_node = zeros(Int, nF)
    for j in 1:N, t in lay.own_Fc[j]; owner_node[t] = j; end
    Epos = Dict(i => t for (t, i) in enumerate(lay.E))
    Fcpos = Dict(i => t for (t, i) in enumerate(lay.Fc))
    rows = rowvals(K)
    E_group = zeros(Int, nE)
    kef = NTuple{3,Int}[]
    for (a, i) in enumerate(lay.E), p in nzrange(K, i)
        t = get(Fcpos, rows[p], 0); t == 0 && continue
        push!(kef, (a, t, p))
        f = feeder[owner_node[t]]; g = f == 0 ? -1 : f
        E_group[a] == 0 || E_group[a] == g || error("E index couples to two feeders")
        E_group[a] = g
    end
    for _ in 1:2, (a, i) in enumerate(lay.E), p in nzrange(K, i)   # follow K_EE
        b = get(Epos, rows[p], 0); b == 0 && continue
        if E_group[a] == 0 && E_group[b] != 0
            E_group[a] = E_group[b]
        elseif E_group[a] != 0 && E_group[b] != 0 && E_group[a] != E_group[b]
            error("K_EE couples two feeders (exact Hessian?): no block structure")
        end
    end
    all(!=(0), E_group) || error("E index with no feeder")
    gkeys = sort!(unique(E_group))
    gidx = Dict(g => t for (t, g) in enumerate(gkeys))
    members = [Int[] for _ in gkeys]
    for a in 1:nE; push!(members[gidx[E_group[a]]], a); end
    perm = vcat(members...)
    granges = UnitRange{Int}[]; off = 0
    loc = zeros(Int, nE)
    for m in members
        push!(granges, off+1:off+length(m)); off += length(m)
        for (t, a) in enumerate(m); loc[a] = t; end
    end
    kee = [NTuple{3,Int}[] for _ in gkeys]
    for (a, i) in enumerate(lay.E), p in nzrange(K, i)
        b = get(Epos, rows[p], 0); b == 0 && continue
        push!(kee[gidx[E_group[a]]], (loc[b], loc[a], p))            # K[E_b, E_a]
    end
    colpos = [Dict{Int,Int}() for _ in 1:N]
    for top in lay.children[lay.root], (s, t) in enumerate(lay.cols[top]); colpos[top][t] = s; end
    gkef = [NTuple{3,Int}[] for _ in gkeys]
    for (a, t, p) in kef
        g = E_group[a]; g == -1 && continue
        push!(gkef[gidx[g]], (loc[a], colpos[g][t], p))
    end
    iso_d = Int[]; iso_e = Vector{NTuple{2,Int}}[]
    giso = [NTuple{5,Int}[] for _ in gkeys]
    for (si, i) in enumerate(lay.iso)
        es = NTuple{2,Int}[]; dp = 0
        for p in nzrange(K, i)
            r = rows[p]
            r == i && (dp = p)
            haskey(Epos, r) && push!(es, (Epos[r], p))
        end
        dp == 0 && error("energy slack without a diagonal entry")
        push!(iso_d, dp); push!(iso_e, es)
        for (a, pa) in es, (b, pb) in es
            E_group[a] == E_group[b] || error("energy slack couples two groups")
            push!(giso[gidx[E_group[a]]], (loc[a], loc[b], pa, pb, si))
        end
    end
    ooff = zeros(Int, N); nN = 0
    for j in 1:N; ooff[j] = nN; nN += length(lay.own[j]); end
    need = [Int[] for _ in 1:N]
    for j in 1:N
        isempty(lay.cols[j]) && continue
        v = copy(lay.own_Fc_pos[j])
        for c in lay.children[j]; isempty(lay.cols[c]) || append!(v, lay.iface_pos[c]); end
        need[j] = sort!(unique(v))
    end
    iface_in_need = [Int[] for _ in 1:N]; fc_in_need = [Int[] for _ in 1:N]
    cidx = [Int[] for _ in 1:N]; depth = zeros(Int, N)
    for j in 1:N
        (isempty(lay.cols[j]) || j == lay.root) && continue
        f = feeder[j]
        cidx[j] = [colpos[f][t] for t in lay.cols[j]]
        fc_in_need[j] = [findfirst(==(pos), need[j]) for pos in lay.own_Fc_pos[j]]
        if j != f
            p = lay.parent[j]
            iface_in_need[j] = [findfirst(==(pos), need[p]) for pos in lay.iface_pos[j]]
        end
        d = 1; k = j
        while k != f; k = lay.parent[k]; d += 1; end
        depth[j] = d
    end
    return TreeKKTStatic(feeder, gkeys, perm, granges, kee, kef, gkef, iso_d, iso_e, giso,
                         colpos, ooff, nN, need, iface_in_need, fc_in_need, cidx, depth)
end

# The dense blocks here are small or mid-sized (a feeder's batteries), where
# OpenBLAS's threading costs more than it gains: a 498 x 498 LU takes 22 ms on
# ten threads and 3 ms on one. Run them on one BLAS thread.
const _BLAS_DEFAULT = Ref(BLAS.get_num_threads())
function _blas1(f)
    n = BLAS.get_num_threads()
    n == 1 && return f()
    _BLAS_DEFAULT[] = n
    BLAS.set_num_threads(1)
    try
        return f()
    finally
        BLAS.set_num_threads(n)
    end
end
# A large dense block inside a _blas1 region: restore the threads for it.
function _blas_full(f, n::Int)
    (n < 1000 || BLAS.get_num_threads() == _BLAS_DEFAULT[]) && return f()
    BLAS.set_num_threads(_BLAS_DEFAULT[])
    try
        return f()
    finally
        BLAS.set_num_threads(1)
    end
end

struct TreeKKTStructured
    perm::Vector{Int}
    granges::Vector{UnitRange{Int}}
    Dfac::Vector{LU{Float64, Matrix{Float64}, Vector{Int}}}
    Up::Matrix{Float64}             # U in the permuted (grouped) order, nE x m_root
    DUp::Matrix{Float64}            # D \ U, same order
    small::LU{Float64, Matrix{Float64}, Vector{Int}}   # M - U' (D \ U)
    M::Matrix{Float64}
    Dfull::Matrix{Float64}          # (verification only) dense D in E order
    # Exact battery curvature: D is no longer block diagonal, so the whole
    # Schur complement S = D + Vc - U inv(M) U' is factored densely instead.
    Sfac::Union{Nothing, LU{Float64, Matrix{Float64}, Vector{Int}}}
end

tree_kkt_structured(lay, stat::TreeKKTStatic, fac::TreeKKTFactor, K::SparseMatrixCSC;
                    dense_check::Bool=false, Vc=nothing) =
    _blas1(() -> _tree_kkt_structured(lay, stat, fac, K, dense_check, Vc))

function _tree_kkt_structured(lay, stat::TreeKKTStatic, fac::TreeKKTFactor, K::SparseMatrixCSC, dense_check::Bool, Vc)
    Bd, Up = _tree_kkt_prepare(lay, stat, fac, K)
    return _tree_kkt_finish(lay, stat, fac, K, Bd, Up, dense_check, Vc)
end

# The structured form in two steps, split by what they read.
#
# _tree_kkt_prepare reads the network factor and K_E,Fc only: the feeder
# inverse blocks Bd and U = K_E,Fc * Psi. Like tree_kkt_factor it never reads
# K_EE, which is the only place the next stage's value function enters a stage
# matrix (the battery-power diagonal under the diagonal Hessian, Vc under the
# exact one). So factor + prepare of every stage can be done before the
# backward sweep, all stages at once (checked by tree_kkt_phase_check.jl,
# which overwrites K_EE with NaN and compares).
#
# _tree_kkt_finish reads K_EE and Vc: the feeder blocks of the Schur
# complement and their factorization. This is what the sweep does in sequence.
function _tree_kkt_prepare(lay, stat::TreeKKTStatic, fac::TreeKKTFactor, K::SparseMatrixCSC)
    N = length(lay.own); nF = length(lay.Fc); nE = length(lay.E)
    root = lay.root; mr = length(lay.own[root])
    nz = nonzeros(K)

    # Upward sweep of the unit Fc columns, stopping below the root.
    Z = Vector{Matrix{Float64}}(undef, N)
    msg = Vector{Matrix{Float64}}(undef, N)
    @inbounds for j in lay.post
        j == root && continue
        k = length(lay.cols[j]); k == 0 && continue
        m = length(lay.own[j])
        R = zeros(m, k)
        for (t, pos) in enumerate(lay.own_Fc_pos[j]); R[pos, t] = 1.0; end
        off = length(lay.own_Fc[j])
        for c in lay.children[j]
            kc = length(lay.cols[c]); kc == 0 && continue
            pos = lay.iface_pos[c]; mc = msg[c]
            for cc in 1:kc, t in eachindex(pos); R[pos[t], off+cc] += mc[t, cc]; end
            off += kc
        end
        Zj = fac.F[j] \ R
        Z[j] = Zj
        mj = smallmul(fac.Aio[j], Zj); mj .*= -1
        msg[j] = mj
    end

    # Downward sweeps per feeder: Bd with zero root values, Psi from unit root
    # values. Each feeder works on its own k_f + q columns.
    Bd = Vector{Matrix{Float64}}(undef, N)
    Psi = zeros(nF, mr)
    @inbounds for top in lay.children[root]
        kf = length(lay.cols[top]); kf == 0 && continue
        colpos = stat.colpos[top]
        Bf = zeros(kf, kf)
        ipos = lay.iface_pos[top]; q = length(ipos)
        # One buffer per depth holds the needed rows of the node being
        # processed at that depth; a node's rows stay valid while its subtree
        # is processed, since only deeper buffers are written meanwhile.
        w = kf + q
        bufs = Matrix{Float64}[]
        ensure(d, nr) = (while length(bufs) < d; push!(bufs, zeros(max(nr, 4), w)); end;
                         size(bufs[d], 1) < nr && (bufs[d] = zeros(nr, w)); bufs[d])
        nd = stat.need[top]; y = ensure(1, length(nd))
        Ztop = Z[top]; Wtop = fac.W[top]
        for (r, pos) in enumerate(nd)
            for cc in 1:kf; y[r, cc] = Ztop[pos, cc]; end
            for cc in 1:q; y[r, kf+cc] = -Wtop[pos, cc]; end
        end
        stack = [top]
        while !isempty(stack)
            j = pop!(stack); d = stat.depth[j]
            if j != top
                # Rows of this node from its parent's (depth d-1), computed
                # now, when the node is processed, so siblings cannot clash.
                yp = bufs[d-1]; ndj = stat.need[j]; yj = ensure(d, length(ndj))
                Wj = fac.W[j]; iin = stat.iface_in_need[j]; Zj = Z[j]; ci = stat.cidx[j]
                for (r, pos) in enumerate(ndj)
                    for cc in 1:w
                        v = 0.0
                        for qq in eachindex(iin); v -= Wj[pos, qq] * yp[iin[qq], cc]; end
                        yj[r, cc] = v
                    end
                    for s in eachindex(ci); yj[r, ci[s]] += Zj[pos, s]; end
                end
            end
            yj = bufs[d]
            for (t, r) in zip(lay.own_Fc[j], stat.fc_in_need[j])
                row = colpos[t]
                for cc in 1:kf; Bf[row, cc] = yj[r, cc]; end
                for cc in 1:q; Psi[t, ipos[cc]] = -yj[r, kf+cc]; end
            end
            for c in lay.children[j]
                isempty(lay.cols[c]) || push!(stack, c)
            end
        end
        Bd[top] = Bf
    end
    @inbounds for (t, pos) in zip(lay.own_Fc[root], lay.own_Fc_pos[root])
        Psi[t, pos] = -1.0
    end

    # U = K_E,Fc * Psi (permuted order).
    invp = invperm(stat.perm)
    Up = zeros(nE, mr)
    @inbounds for (a, t, p) in stat.kef
        v = nz[p]; r = invp[a]
        for c in 1:mr; Up[r, c] += v * Psi[t, c]; end
    end
    return Bd, Up
end

function _tree_kkt_finish(lay, stat::TreeKKTStatic, fac::TreeKKTFactor, K::SparseMatrixCSC, Bd, Up,
                          dense_check::Bool, Vc)
    nE = length(lay.E); mr = length(lay.own[lay.root])
    nz = nonzeros(K)
    invp = invperm(stat.perm)
    # The diagonal blocks D_g.
    ng = length(stat.gkeys)
    Dfac = Vector{LU{Float64, Matrix{Float64}, Vector{Int}}}(undef, ng)
    DUp = similar(Up)
    Dfull = dense_check ? zeros(nE, nE) : zeros(0, 0)
    exact = !isnothing(Vc)
    Sdense = exact ? zeros(nE, nE) : zeros(0, 0)      # permuted (grouped) order
    @inbounds for gi in 1:ng
        rg = stat.granges[gi]; n = length(rg)
        Dg = zeros(n, n)
        for (b, a, p) in stat.kee[gi]; Dg[b, a] = nz[p]; end
        if stat.gkeys[gi] != -1
            Bf = Bd[stat.gkeys[gi]]
            ge = stat.gkef[gi]
            for (a, fa, pa) in ge, (b, fb, pb) in ge
                Dg[a, b] -= nz[pa] * Bf[fa, fb] * nz[pb]
            end
        end
        for (a, b, pa, pb, si) in stat.giso[gi]
            Dg[a, b] -= nz[pa] * nz[pb] / nz[stat.iso_d[si]]
        end
        dense_check && (Dfull[stat.perm[rg], stat.perm[rg]] .= Dg)
        if exact
            Sdense[rg, rg] .= Dg
            continue
        end
        F = lu!(Dg)
        Dfac[gi] = F
        DUp[rg, :] = F \ Up[rg, :]
    end
    if exact
        nB = size(Vc, 1)
        (size(Vc) == (nB, nB) && nB <= nE) || error("battery curvature has the wrong size")
        @inbounds for b in 1:nB, a in 1:nB
            Sdense[invp[a], invp[b]] += Vc[a, b]      # battery powers are the first n_B of E
        end
        Sdense .-= Up * (fac.Mroot \ Up')
        Sfac = _blas_full(() -> lu!(Sdense), nE)
        return TreeKKTStructured(stat.perm, stat.granges, Dfac, Up, DUp, lu(fac.Mroot), fac.Mroot, Dfull, Sfac)
    end
    small = lu!(fac.Mroot - Up' * DUp)
    return TreeKKTStructured(stat.perm, stat.granges, Dfac, Up, DUp, small, fac.Mroot, Dfull, nothing)
end

# S \ R for S = D - U inv(M) U' (Woodbury): X = D\R + D\U * ((M - U' D\U) \ (U' D\R)).
# Works in the grouped order so that each group is a contiguous row block.
function tree_kkt_schur_solve(st::TreeKKTStructured, R::AbstractMatrix)
    Rp = R[st.perm, :]
    if !isnothing(st.Sfac)
        ldiv!(st.Sfac, Rp)
        X = Matrix{Float64}(undef, size(Rp))
        X[st.perm, :] = Rp
        return X
    end
    _blas1() do
        for (g, rg) in enumerate(st.granges)
            ldiv!(st.Dfac[g], view(Rp, rg, :))
        end
    end
    T = st.small \ (st.Up' * Rp)
    mul!(Rp, st.DUp, T, 1.0, 1.0)
    X = Matrix{Float64}(undef, size(Rp))
    X[st.perm, :] = Rp
    return X
end

function tree_kkt_schur_dense(st::TreeKKTStructured)
    invp = invperm(st.perm)
    U = st.Up[invp, :]
    return st.Dfull - U * (st.M \ U')
end

# ------------------------------------------------------------- full solver --
# One object per stage that replaces the stage's sparse LU: full solves
# K \ b (feedforward column, forward-pass policy) and the battery rows of the
# feedback columns. With N the network part and E last,
#     y   = N \ b_F
#     x_E = S \ (b_E - K_EF y - K_E,iso (b_iso ./ d_iso))
#     x_F = y - N \ (K_FE x_E),   x_iso = (b_iso - K_iso,E x_E) ./ d_iso.

struct TreeKKTSolver
    lay::TreeKKTLayout
    stat::TreeKKTStatic
    fac::TreeKKTFactor
    st::TreeKKTStructured
    nz::Vector{Float64}                    # copy of the K values the couplings are read from
    pre::Vector{Int}                       # nodes, parents before children
    z::Vector{Float64}                     # flat work buffer for the network solves
    y1::Vector{Float64}
    y2::Vector{Float64}
    t2::Vector{Float64}
    # Seconds spent, by what the work needs from the next stage (see
    # _tree_kkt_prepare): [1] factor + prepare (nothing), [2] finish (the value
    # function), and accumulated over the solves [3] first network solve
    # (nothing), [4] battery block solve, [5] second network solve (this
    # stage's battery solution only).
    times::Vector{Float64}
end

function tree_kkt_solver(lay::TreeKKTLayout, stat::TreeKKTStatic, K::SparseMatrixCSC; Vc=nothing)
    t0 = time_ns()
    fac = tree_kkt_factor(lay, K)
    Bd, Up = _blas1(() -> _tree_kkt_prepare(lay, stat, fac, K))
    nzc = copy(nonzeros(K)); pre = reverse(lay.post)
    z = zeros(stat.nN); y1 = zeros(lay.n); y2 = zeros(lay.n); t2 = zeros(lay.n)
    t1 = time_ns()
    st = _blas1(() -> _tree_kkt_finish(lay, stat, fac, K, Bd, Up, false, Vc))
    times = [(t1 - t0) / 1e9, (time_ns() - t1) / 1e9, 0.0, 0.0, 0.0]
    return TreeKKTSolver(lay, stat, fac, st, nzc, pre, z, y1, y2, t2, times)
end
tree_kkt_solver(lay::TreeKKTLayout, K::SparseMatrixCSC) = tree_kkt_solver(lay, tree_kkt_static(lay, K), K)

# x[network rows] = N \ b[network rows], with no allocation: z holds every
# node's own block back to back, and must be zero on entry (it is on exit).
function tree_kkt_network_solve!(x::AbstractVector, s::TreeKKTSolver, b::AbstractVector)
    lay, fac, z, ooff = s.lay, s.fac, s.z, s.stat.ooff
    @inbounds for j in lay.post
        own = lay.own[j]; m = length(own); o = ooff[j]
        F = fac.F[j]; A = F.A; perm = F.perm
        # right-hand side = b + the children's messages (accumulated in z), pivoted
        for t in 1:m; x[own[t]] = b[own[t]] + z[o+t]; end
        for t in 1:m; z[o+t] = x[own[perm[t]]]; end
        for i in 2:m
            v = z[o+i]
            for k in 1:i-1; v -= A[i, k] * z[o+k]; end
            z[o+i] = v
        end
        for i in m:-1:1
            v = z[o+i]
            for k in i+1:m; v -= A[i, k] * z[o+k]; end
            z[o+i] = v / A[i, i]
        end
        if j != lay.root
            po = ooff[lay.parent[j]]; pos = lay.iface_pos[j]; Bio = fac.Aio[j]
            for q in eachindex(pos)
                v = 0.0
                for t in 1:m; v += Bio[q, t] * z[o+t]; end
                z[po+pos[q]] -= v                                # message to the parent
            end
        end
    end
    @inbounds for j in s.pre
        own = lay.own[j]; m = length(own); o = ooff[j]
        if j != lay.root
            po = ooff[lay.parent[j]]; pos = lay.iface_pos[j]; W = fac.W[j]
            for q in eachindex(pos)
                yI = z[po+pos[q]]
                yI == 0 && continue
                for t in 1:m; z[o+t] -= W[t, q] * yI; end
            end
        end
        for t in 1:m; x[own[t]] = z[o+t]; end
    end
    fill!(z, 0.0)
    return x
end

function tree_kkt_solve!(s::TreeKKTSolver, b::AbstractVector)
    lay, stat, nz = s.lay, s.stat, s.nz
    y, y2, t2 = s.y1, s.y2, s.t2
    t_a = time_ns()
    tree_kkt_network_solve!(y, s, b)
    t_b = time_ns()
    nE = length(lay.E)
    rE = Matrix{Float64}(undef, nE, 1)
    @inbounds for a in 1:nE; rE[a, 1] = b[lay.E[a]]; end
    @inbounds for (a, t, p) in stat.kef; rE[a, 1] -= nz[p] * y[lay.Fc[t]]; end
    @inbounds for (si, i) in enumerate(lay.iso)
        w = b[i] / nz[stat.iso_d[si]]
        for (a, p) in stat.iso_e[si]; rE[a, 1] -= nz[p] * w; end
    end
    xE = tree_kkt_schur_solve(s.st, rE)
    t_c = time_ns()
    fill!(t2, 0.0)
    @inbounds for (a, f, p) in stat.kef; t2[lay.Fc[f]] += nz[p] * xE[a, 1]; end
    tree_kkt_network_solve!(y2, s, t2)
    @inbounds for (si, i) in enumerate(lay.iso)
        v = b[i]
        for (a, p) in stat.iso_e[si]; v -= nz[p] * xE[a, 1]; end
        b[i] = v / nz[stat.iso_d[si]]
    end
    @inbounds for j in eachindex(lay.own), i in lay.own[j]; b[i] = y[i] - y2[i]; end
    @inbounds for a in 1:nE; b[lay.E[a]] = xE[a, 1]; end
    tm = s.times
    tm[3] += (t_b - t_a) / 1e9; tm[4] += (t_c - t_b) / 1e9; tm[5] += (time_ns() - t_c) / 1e9
    return b
end

function LinearAlgebra.ldiv!(s::TreeKKTSolver, B::AbstractVecOrMat)
    if B isa AbstractVector
        tree_kkt_solve!(s, B)
    else
        for c in axes(B, 2); tree_kkt_solve!(s, @view B[:, c]); end
    end
    return B
end

# Battery rows (in lay.E order) of K \ R for right-hand sides supported on E.
tree_kkt_battery_rows(s::TreeKKTSolver, RE::AbstractMatrix) = tree_kkt_schur_solve(s.st, RE)

# Driver hook for FILTERDDP_TREE_KKT=1. The layout depends only on the network
# and on K's sparsity pattern, so it is built once and rebuilt only if the
# pattern changes.
function install_tree_kkt_hook(data, idx, nu::Int)
    state = Ref{Any}(nothing)
    DDP4OPF.TREE_KKT_HOOK[] = function (K, Vc=nothing)
        st = state[]
        if isnothing(st) || st.colptr != K.colptr || st.rowval != K.rowval
            lay = tree_kkt_layout(data, idx, nu, K)
            st = (lay=lay, stat=tree_kkt_static(lay, K), colptr=copy(K.colptr), rowval=copy(K.rowval))
            state[] = st
            println("TREE_KKT layout: nodes=", length(lay.own), " feeders=",
                    length(lay.children[lay.root]), " nE=", length(lay.E),
                    " groups=", length(st.stat.gkeys),
                    " max_own=", maximum(length, lay.own), " max_iface=", maximum(length, lay.iface))
        end
        return tree_kkt_solver(st.lay, st.stat, K; Vc=Vc)
    end
    return nothing
end

function DDP4OPF.battery_block_rows(F::TreeKKTSolver, K, E, R)
    E == F.lay.E || error("battery rows differ from the tree layout's")
    return tree_kkt_battery_rows(F, R)
end

DDP4OPF.stage_phase_times(F::TreeKKTSolver) = F.times
