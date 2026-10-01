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
# with all feeder contributions (already factored in tree_kkt_factor), and
# Psi the feeder rows' response to the root unknowns (-1 on root-owned rows).
# Under the diagonal Hessian K_EE couples a battery only to its own energy
# row, so the Schur complement onto E is block diagonal by feeder plus a
# rank-(root size) term:
#     S = D - U * inv(M) * U',   U = K_E,Fc * Psi,
# and S \ R follows from Woodbury with one small solve per feeder.

function _feeder_of(lay::TreeKKTLayout)
    N = length(lay.own)
    feeder = zeros(Int, N)                      # 0 for the root
    stack = [(c, c) for c in lay.children[lay.root]]
    while !isempty(stack)
        j, f = pop!(stack)
        feeder[j] = f
        for c in lay.children[j]; push!(stack, (c, f)); end
    end
    return feeder
end

struct TreeKKTStructured
    groups::Vector{Vector{Int}}     # positions in E of each group (feeders, then root)
    Dfac::Vector{Any}               # LU of each group's diagonal block
    U::Matrix{Float64}              # nE x m_root
    M::Matrix{Float64}              # the root block M (small)
    Dfull::Matrix{Float64}          # (verification only) dense D; empty in production
end

function tree_kkt_structured(lay::TreeKKTLayout, fac::TreeKKTFactor, K::SparseMatrixCSC; dense_check::Bool=false)
    N = length(lay.own); nF = length(lay.Fc); nE = length(lay.E)
    root = lay.root; mr = length(lay.own[root])
    feeder = _feeder_of(lay)
    owner_node = zeros(Int, nF)
    for j in 1:N, t in lay.own_Fc[j]; owner_node[t] = j; end
    Fc_feeder = [feeder[owner_node[t]] for t in 1:nF]

    # Upward sweep of the unit Fc columns, stopping below the root.
    Z = Vector{Matrix{Float64}}(undef, N)
    msg = Vector{Matrix{Float64}}(undef, N)
    for j in lay.post
        j == root && continue
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
        Z[j] = fac.F[j] \ R
        msg[j] = smallmul(fac.Aio[j], Z[j]); msg[j] .*= -1
    end

    # Downward sweeps per feeder: Bd (feeder-only inverse at its Fc rows) with
    # zero root values, and Psi from unit root values (zero right-hand side).
    Bd = Dict{Int, Matrix{Float64}}()         # feeder top => k_f x k_f in cols[top] order
    Psi = zeros(nF, mr)
    for top in lay.children[root]
        kf = length(lay.cols[top]); kf == 0 && continue
        colpos = Dict(t => s for (s, t) in enumerate(lay.cols[top]))
        Bf = zeros(kf, kf)
        ipos = lay.iface_pos[top]               # root positions this feeder touches
        q = length(ipos)
        # y (own x (kf + q)): first kf columns = Fc unit responses, last q = root basis responses
        ytop = zeros(length(lay.own[top]), kf + q)
        @views ytop[:, 1:kf] .= Z[top]
        @views ytop[:, kf+1:end] .= -fac.W[top]          # y = -W * e_k
        stack = [(top, ytop)]
        while !isempty(stack)
            j, yj = pop!(stack)
            for (t, pos) in zip(lay.own_Fc[j], lay.own_Fc_pos[j])
                @views Bf[colpos[t], :] .= yj[pos, 1:kf]
                @views Psi[t, ipos] .= .-yj[pos, kf+1:end]  # Psi = -(response)
            end
            for c in lay.children[j]
                isempty(lay.cols[c]) && continue
                yc = smallmul(fac.W[c], @view yj[lay.iface_pos[c], :]); yc .*= -1
                cidx = [colpos[t] for t in lay.cols[c]]
                @views yc[:, cidx] .+= Z[c]
                push!(stack, (c, yc))
            end
        end
        Bd[top] = Bf
    end
    for (t, pos) in zip(lay.own_Fc[root], lay.own_Fc_pos[root])
        Psi[t, pos] = -1.0
    end

    # Group E by feeder through each E index's Fc coupling (energy rows follow
    # their battery power through K_EE).
    Epos = Dict(i => t for (t, i) in enumerate(lay.E))
    Fcpos = Dict(i => t for (t, i) in enumerate(lay.Fc))
    rows, vals = rowvals(K), nonzeros(K)
    E_group = zeros(Int, nE)                     # feeder top, or -1 for the root group
    for (a, i) in enumerate(lay.E), p in nzrange(K, i)
        t = get(Fcpos, rows[p], 0); t == 0 && continue
        g = Fc_feeder[t] == 0 ? -1 : Fc_feeder[t]
        E_group[a] == 0 || E_group[a] == g || error("E index couples to two feeders")
        E_group[a] = g
    end
    for _ in 1:2, (a, i) in enumerate(lay.E), p in nzrange(K, i)  # propagate along K_EE
        b = get(Epos, rows[p], 0); b == 0 && continue
        if E_group[a] == 0 && E_group[b] != 0
            E_group[a] = E_group[b]
        elseif E_group[a] != 0 && E_group[b] != 0 && E_group[a] != E_group[b]
            error("K_EE couples two feeders (exact Hessian?): no block structure")
        end
    end
    all(!=(0), E_group) || error("E index with no feeder")
    gkeys = sort!(unique(E_group))
    groups = [findall(==(g), E_group) for g in gkeys]

    # U = K_E,Fc * Psi and the diagonal blocks D_g.
    KEF = K[lay.E, lay.Fc]
    U = Matrix(KEF * Psi)
    Dfac = Vector{Any}(undef, length(groups))
    Dfull = dense_check ? zeros(nE, nE) : zeros(0, 0)
    for (gi, g) in enumerate(gkeys)
        ga = groups[gi]
        Dg = Matrix(K[lay.E[ga], lay.E[ga]])
        if g != -1
            fcs = lay.cols[g]                        # Fc positions of the feeder, Bd order
            Kg = Matrix(KEF[ga, fcs])
            Dg .-= Kg * Bd[g] * Kg'
        end
        # energy-slack terms: each slack couples to one energy row
        gaset = Dict(a => s for (s, a) in enumerate(ga))
        for i in lay.iso
            nzE = [(gaset[Epos[rows[p]]], vals[p]) for p in nzrange(K, i)
                   if haskey(Epos, rows[p]) && haskey(gaset, Epos[rows[p]])]
            isempty(nzE) && continue
            d = K[i, i]
            for (a, va) in nzE, (b, vb) in nzE; Dg[a, b] -= va * vb / d; end
        end
        dense_check && (Dfull[ga, ga] .= Dg)
        Dfac[gi] = lu(Dg)
    end
    return TreeKKTStructured(groups, Dfac, U, fac.Mroot, Dfull)
end

# S \ R for S = D - U inv(M) U' (Woodbury): X = D\R + D\U * ((M - U' D\U) \ (U' D\R)).
function tree_kkt_schur_solve(st::TreeKKTStructured, R::AbstractMatrix)
    DR = Matrix{Float64}(undef, size(R)); DU = similar(st.U)
    for (g, ga) in enumerate(st.groups)
        F = st.Dfac[g]::LU{Float64, Matrix{Float64}, Vector{Int}}
        DR[ga, :] = F \ R[ga, :]
        DU[ga, :] = F \ st.U[ga, :]
    end
    small = st.M - st.U' * DU
    return DR + DU * (small \ (st.U' * DR))
end

tree_kkt_schur_dense(st::TreeKKTStructured) =
    st.Dfull - st.U * (st.M \ st.U')

# ------------------------------------------------------------- full solver --
# One object per stage that replaces the stage's sparse LU: full solves
# K \ b (feedforward column, forward-pass policy) and the battery rows of the
# feedback columns. With N the network part and E last,
#     y   = N \ b_F
#     x_E = S \ (b_E - K_EF y - K_E,iso (b_iso ./ d_iso))
#     x_F = y - N \ (K_FE x_E),   x_iso = (b_iso - K_iso,E x_E) ./ d_iso.

struct TreeKKTSolver
    lay::TreeKKTLayout
    fac::TreeKKTFactor
    st::TreeKKTStructured
    KEF::SparseMatrixCSC{Float64,Int}      # E x Fc
    KEiso::SparseMatrixCSC{Float64,Int}    # E x iso
    diso::Vector{Float64}
    pre::Vector{Int}                       # nodes, parents before children
end

function tree_kkt_solver(lay::TreeKKTLayout, K::SparseMatrixCSC)
    fac = tree_kkt_factor(lay, K)
    st = tree_kkt_structured(lay, fac, K)
    return TreeKKTSolver(lay, fac, st, K[lay.E, lay.Fc], K[lay.E, lay.iso],
                         Float64[K[i, i] for i in lay.iso], reverse(lay.post))
end

# x[network rows] = N \ b[network rows]; other entries of x untouched.
function tree_kkt_network_solve!(x::AbstractVector, s::TreeKKTSolver, b::AbstractVector)
    lay, fac = s.lay, s.fac
    N = length(lay.own)
    z = Vector{Matrix{Float64}}(undef, N)
    msg = Vector{Matrix{Float64}}(undef, N)
    @inbounds for j in lay.post
        own = lay.own[j]; m = length(own)
        r = Matrix{Float64}(undef, m, 1)
        for t in 1:m; r[t, 1] = b[own[t]]; end
        for c in lay.children[j]
            pos = lay.iface_pos[c]; mc = msg[c]
            for t in eachindex(pos); r[pos[t], 1] += mc[t, 1]; end
        end
        z[j] = fac.F[j] \ r
        if j != lay.root
            mj = smallmul(fac.Aio[j], z[j]); mj .*= -1
            msg[j] = mj
        end
    end
    @inbounds for j in s.pre
        own = lay.own[j]; zj = z[j]
        if j != lay.root
            p = lay.parent[j]; pos = lay.iface_pos[j]; W = fac.W[j]; zp = z[p]
            for q in eachindex(pos)
                yI = zp[pos[q], 1]
                yI == 0 && continue
                for t in eachindex(own); zj[t, 1] -= W[t, q] * yI; end
            end
        end
        for t in eachindex(own); x[own[t]] = zj[t, 1]; end
    end
    return x
end

function tree_kkt_solve!(s::TreeKKTSolver, b::AbstractVector)
    lay = s.lay
    y = zeros(lay.n)
    tree_kkt_network_solve!(y, s, b)
    bE = b[lay.E]; biso = b[lay.iso]
    rE = bE - s.KEF * y[lay.Fc] - s.KEiso * (biso ./ s.diso)
    xE = vec(tree_kkt_schur_solve(s.st, reshape(rE, :, 1)))
    t = zeros(lay.n)
    t[lay.Fc] .= s.KEF' * xE
    y2 = zeros(lay.n)
    tree_kkt_network_solve!(y2, s, t)
    xiso = (biso .- s.KEiso' * xE) ./ s.diso
    for j in eachindex(lay.own), i in lay.own[j]; b[i] = y[i] - y2[i]; end
    b[lay.E] .= xE
    b[lay.iso] .= xiso
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
    DDP4OPF.TREE_KKT_HOOK[] = function (K)
        st = state[]
        if isnothing(st) || st.colptr != K.colptr || st.rowval != K.rowval
            lay = tree_kkt_layout(data, idx, nu, K)
            st = (lay=lay, colptr=copy(K.colptr), rowval=copy(K.rowval))
            state[] = st
            println("TREE_KKT layout: nodes=", length(lay.own), " feeders=",
                    length(lay.children[lay.root]), " nE=", length(lay.E),
                    " max_own=", maximum(length, lay.own), " max_iface=", maximum(length, lay.iface))
        end
        return tree_kkt_solver(st.lay, K)
    end
    return nothing
end

function DDP4OPF.battery_block_rows(F::TreeKKTSolver, K, E, R)
    E == F.lay.E || error("battery rows differ from the tree layout's")
    return tree_kkt_battery_rows(F, R)
end
