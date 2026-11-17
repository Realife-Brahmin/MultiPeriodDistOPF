# FILTERDDP_VALUE_BLOCKS=1: the battery value function as one small block per
# feeder, where the network makes that exact. Opt-in; off by default.
#
# The pairwise pass (paired_pass.jl, procs_pass.jl) carries each period as a
# function G_t(a, b) of the battery energies, three n_B x n_B blocks, and every
# combination and every stage's battery block costs n_B^3. On large10k (1020
# batteries) that is 60-80% of a measured solve.
#
# Why the blocks are exact. The tree solver already writes the halved battery
# block of a stage as (tree_kkt.jl, _tree_kkt_finish0)
#
#     Shat0 = blockdiag_g(H0_g) - U_P inv(M) U_P',
#
# one block H0_g per feeder below the substation, and a term through the root
# block M (the substation's own unknowns with every feeder's contribution W
# folded in). A feeder reaches the root only through the root unknowns of its
# interface, I: the substation voltage and the multipliers of its two balances.
# The other root unknowns, Z, are the substation's injections P_Subs, Q_Subs
# (and the multiplier of the fixed voltage). When the injections are free and
# their cost is linear, Z has no diagonal entry, and the root block is
#
#     M = [ 0   J ]        inv(M) = [ *   inv(J') ]
#         [ J'  W ]                 [ inv(J)  0   ]
#
# with J square and nonsingular: each interface unknown is pinned by one free
# unknown. The (I, I) block of inv(M) is zero whatever W is, so the root term
# vanishes and the stage's battery block is block diagonal by feeder. A sum, a
# Schur complement and an inverse of matrices that are block diagonal in the
# same partition are block diagonal in it, so every G_t, every combination and
# every value function is too. Physically: a slack bus that holds its voltage
# and trades any amount at a fixed price decouples the feeders hanging off it.
#
# It stops being exact when the substation's injection carries curvature: an
# active or barrier-carrying bound on P_Subs (FILTERDDP_SCREEN without psubs),
# or a quadratic cost on it. Nothing is assumed here: every stage of every
# pass computes a bound on the root term it leaves out (relative to the
# largest block entry) and the stage fails if that is above
# FILTERDDP_VALUE_BLOCKS_TOL (default 1e-10), with a BLOCK_VALUE_FAIL line.
# The paired pass then hands that pass to the solver's own sequential one
# (correct, slow); the process pass stops the solve with status 1.
#
# What changes with it:
#   - the stage's battery block is never assembled at n_B x n_B: the blocks
#     H0_g are kept (bv_finish0);
#   - a period's function is three n_g x n_g blocks per feeder (BVElem), in
#     the first part of the same storage slot the dense blocks would use;
#   - a combination is one small combination per feeder (the same kernel,
#     _pp_combine!, on each block);
#   - a stage's battery solve factors one small block per feeder and uses the
#     tree solver's structured solve (per-feeder factors and a Woodbury
#     correction through the root), as the diagonal-Hessian path does;
#   - with blocks this small a combination costs less than handing it to
#     another process, so the value functions are built by the ordinary
#     recursion from the last period to the first, in one process
#     (FILTERDDP_VALUE_BLOCKS_SWEEP_MAX, default 128: the largest block for
#     which that is done; above it the rounds are used, block by block).
#
# With one feeder (ieee123, med2522) there is one block of all the batteries:
# the switch is then inactive and the dense code runs unchanged.
#
# Included at the end of paired_pass.jl.

_bv_requested() = get(ENV, "FILTERDDP_VALUE_BLOCKS", "0") != "0"
_bv_tol() = parse(Float64, get(ENV, "FILTERDDP_VALUE_BLOCKS_TOL", "1e-10"))
_bv_sweep_max() = parse(Int, get(ENV, "FILTERDDP_VALUE_BLOCKS_SWEEP_MAX", "128"))
_bv_audit() = get(ENV, "FILTERDDP_VALUE_BLOCKS_AUDIT", "0") != "0"

const _BV_ACTIVE = Ref(false)             # decided at install, from the network
const _BV_N = Ref(0)                      # periods
const _BV_SHARED = Ref{Any}(nothing)      # the periods' shared segments (procs), or nothing
const _BV_LAYOUT = Ref{Any}(nothing)      # set by the first stage
const _BV_SLOTS = Any[]                   # one BVElem per period
const _BV_BUFS = Vector{Float64}[]        # their storage when it is not shared
const _BV_END = Ref{Any}(nothing)         # the zero value function after the last period
const _BV_LOCK = ReentrantLock()
const _BV_DROP = Ref(0.0)                 # largest relative bound on the root term left out, this process
const _BV_FAILED = Ref(false)

struct BVLayout
    idx::Vector{Vector{Int}}              # per block: its batteries (states), in the order of the block's rows
    gid::Vector{Int}                      # per block: its group in the tree solver
    gof::Vector{Int}                      # battery -> block
    lof::Vector{Int}                      # battery -> row in its block
    off::Vector{Int}                      # block -> offset of its numbers in a period's storage
    len::Int                              # numbers per period
    bat::Vector{Int}                      # the tree solver's battery of each grouped power row (for checking)
    nP::Vector{Int}
end

# A period's function, block by block. el[g] wraps the block's storage; once
# the period is closed with the value function after it, val[g] (el[g].Maa and
# el[g].ma) is the value function before the period.
struct BVElem
    el::Vector{PintElem}
    val::Vector{PintVal}
end
struct BVVal
    val::Vector{PintVal}
end

# The halved battery block of a stage without the next stage's curvature, as
# its feeder blocks. The fields that _pp_elem reads have TreeKKTBlock0's names.
struct BVBlock0
    H0::Vector{Matrix{Float64}}           # per group: D_PP - D_PN inv(D_NN) D_NP
    bat::Vector{Int}
    nP::Vector{Int}
    Ninv::Vector{Vector{Float64}}
    PN::Vector{SparseMatrixCSC{Float64, Int}}
    NP::Vector{SparseMatrixCSC{Float64, Int}}
    prow::Vector{Int}
    Up::Matrix{Float64}
    cols::Vector{Int}                     # the columns of Up that are not zero (root unknowns the feeders reach)
    drop::Float64                         # bound on the root term left out, relative to the largest block entry
end

function _bv_fail(msg::String)
    if !_BV_FAILED[]
        println("BLOCK_VALUE_FAIL ", msg, "; rerun with FILTERDDP_VALUE_BLOCKS=0")
        _BV_FAILED[] = true
    end
    return nothing
end

# Feeders (branches directly below the substation) that hold batteries, from
# the network data alone. More than one, and no battery at the substation
# itself: the blocks are worth having.
function bv_groups_from_data(data)
    sub = data[:substationBus]; par = data[:parent]
    tops = Int[]
    for b in data[:Bset]
        b == sub && return 0, 0
        k = b
        while par[k] != sub; k = par[k]; end
        push!(tops, k)
    end
    u = unique(tops)
    return length(u), maximum(count(==(k), tops) for k in u)
end

# Numbers per period with one block per feeder (the size of a shared segment).
function bv_block_length(data)
    sub = data[:substationBus]; par = data[:parent]
    cnt = Dict{Int, Int}()
    for b in data[:Bset]
        k = b
        while par[k] != sub; k = par[k]; end
        cnt[k] = get(cnt, k, 0) + 1
    end
    return sum(3m * m + 2m for m in values(cnt))
end

# N: the number of periods if known here (0: paired_pass sets it with its storage).
function bv_install!(data, N::Int, shared; quiet::Bool=false)
    _BV_LAYOUT[] = nothing; empty!(_BV_SLOTS); empty!(_BV_BUFS); _BV_END[] = nothing
    _BV_DROP[] = 0.0; _BV_FAILED[] = false
    _BV_N[] = N; _BV_SHARED[] = shared
    ng, mx = bv_groups_from_data(data)
    fit = _paired_mode() != "0" && _paired_inplace() && _paired_rounds() && _half_block() && SPARDDP._lean_exact_enabled()
    _BV_ACTIVE[] = _bv_requested() && ng > 1 && fit
    (_bv_requested() && !quiet) && println("BLOCK_VALUE requested: feeders with batteries=", ng, " largest=", mx,
        " batteries=", length(data[:Bset]), " active=", _BV_ACTIVE[],
        _BV_ACTIVE[] ? "" : (ng <= 1 ? " (one block: the dense pass runs unchanged)" : " (needs the pairwise in-place pass with the halved exact block)"))
    return _BV_ACTIVE[]
end

# ---------------------------------------------------------------- the stage --

# _tree_kkt_finish0 without the n_B x n_B assembly: the same group blocks, the
# same conditions, and a bound on the root term in place of the term.
function bv_finish0(lay, stat::TreeKKTStatic, fac, K::SparseMatrixCSC, Bd, Up, nB::Int)
    nE = length(lay.E); nz = nonzeros(K)
    ng = length(stat.gkeys)
    any(==(-1), stat.gkeys) && return _bv_fail("a battery row couples to the substation directly")
    Dgs = Vector{Matrix{Float64}}(undef, ng)
    @inbounds for gi in 1:ng
        rg = stat.granges[gi]; n = length(rg)
        Dg = zeros(n, n)
        for (b, a, p) in stat.kee[gi]; Dg[b, a] = nz[p]; end
        Bf = Bd[stat.gkeys[gi]]
        ge = stat.gkef[gi]
        for (a, fa, pa) in ge, (b, fb, pb) in ge
            Dg[a, b] -= nz[pa] * Bf[fa, fb] * nz[pb]
        end
        for (a, b, pa, pb, si) in stat.giso[gi]
            Dg[a, b] -= nz[pa] * nz[pb] / nz[stat.iso_d[si]]
        end
        Dgs[gi] = Dg
    end
    isP = falses(nE)
    for (a, t, p) in stat.kef; isP[a] = true; end
    nPv = zeros(Int, ng)
    for gi in 1:ng
        rg = stat.granges[gi]; n = length(rg); Dg = Dgs[gi]
        np = count(l -> isP[stat.perm[rg[l]]], 1:n)
        ok = np >= 1 && all(l -> isP[stat.perm[rg[l]]] == (l <= np), 1:n)
        for l2 in np+1:n, l1 in np+1:n
            ok &= l1 == l2 ? Dg[l1, l1] != 0 : Dg[l1, l2] == 0
        end
        ok &= all(iszero, @view Up[rg[np+1:n], :])
        ok || return _bv_fail("the battery block of a feeder cannot be halved")
        nPv[gi] = np
    end
    Ninv = Vector{Float64}[]; prow = Int[]
    PN = SparseMatrixCSC{Float64, Int}[]; NP = SparseMatrixCSC{Float64, Int}[]
    H0 = Vector{Matrix{Float64}}(undef, ng)
    hmax = 0.0
    for gi in 1:ng
        rg = stat.granges[gi]; n = length(rg); np = nPv[gi]; Dg = Dgs[gi]
        push!(Ninv, [1 / Dg[l, l] for l in np+1:n])
        push!(PN, sparse(Dg[1:np, np+1:n])); push!(NP, sparse(Dg[np+1:n, 1:np]))
        append!(prow, rg[1:np])
        PN[gi] == sparse(NP[gi]') || return _bv_fail("the energy rows of a feeder are not symmetric")
        blk = Dg[1:np, 1:np]
        blk .-= Matrix(PN[gi] * (Diagonal(Ninv[gi]) * NP[gi]))
        H0[gi] = blk
        @inbounds for v in blk; hmax = max(hmax, abs(v)); end
    end
    n = length(prow)
    bat = [stat.perm[prow[k]] for k in 1:n]
    (n == nB && all(<=(nB), bat)) || return _bv_fail("the power rows are not the batteries")
    # The root term U_P inv(M) U_P' is left out. Its entry (r, s) is
    # (U_P C)[r, :] . U_P[s, :] with C the block of inv(M) on the root unknowns
    # the feeders reach, so its largest entry is at most the product of the
    # largest row norms of U_P C and of U_P.
    mr = size(Up, 2)
    cols = [c for c in 1:mr if any(!iszero, @view Up[:, c])]
    drop = 0.0
    if !isempty(cols)
        C = inv(Matrix(fac.Mroot))[cols, cols]
        a = 0.0; b = 0.0; q = length(cols)
        tmp = zeros(q)
        @inbounds for r in prow
            fill!(tmp, 0.0); s2 = 0.0
            for (ci, c) in enumerate(cols)
                u = Up[r, c]; s2 += u * u
                u == 0 && continue
                for cj in 1:q; tmp[cj] += u * C[ci, cj]; end
            end
            b = max(b, s2)
            t2 = 0.0
            for cj in 1:q; t2 += tmp[cj] * tmp[cj]; end
            a = max(a, t2)
        end
        drop = sqrt(a) * sqrt(b) / max(hmax, floatmin())
    end
    isfinite(drop) && drop <= _bv_tol() ||
        return _bv_fail(@sprintf("the feeders are coupled through the substation (root term %.3e of the largest block entry, tolerance %.1e)", drop, _bv_tol()))
    drop > _BV_DROP[] && (_BV_DROP[] = drop)
    return BVBlock0(H0, bat, nPv, Ninv, PN, NP, prow, Up, cols, drop)
end

function _bv_wrap(p::Ptr{Float64}, L::BVLayout)
    el = PintElem[]; val = PintVal[]
    for (bi, st) in enumerate(L.idx)
        m = length(st); q = p + L.off[bi] * sizeof(Float64)
        blk = i -> unsafe_wrap(Array, q + (i - 1) * m * m * sizeof(Float64), (m, m))
        vec = i -> unsafe_wrap(Array, q + (3m * m + (i - 1) * m) * sizeof(Float64), (m,))
        e = PintElem(blk(1), blk(2), blk(3), vec(1), vec(2))
        push!(el, e); push!(val, PintVal(e.Maa, e.ma))
    end
    return BVElem(el, val)
end

function _bv_make_layout(bat::Vector{Int}, nP::Vector{Int}, nB::Int)
    idx = Vector{Int}[]; gid = Int[]; off = Int[]
    gof = zeros(Int, nB); lof = zeros(Int, nB)
    o = 0; len = 0
    for (gi, np) in enumerate(nP)
        st = bat[o+1:o+np]; o += np
        push!(idx, st); push!(gid, gi); push!(off, len)
        for (k, i) in enumerate(st); gof[i] = length(idx); lof[i] = k; end
        len += 3np * np + 2np
    end
    all(!=(0), gof) || error("block value pass: a battery is in no feeder block")
    return BVLayout(idx, gid, gof, lof, off, len, copy(bat), copy(nP))
end

# The storage of every period in this layout: the first numbers of the shared
# segments (procs), or arrays of its own.
function _bv_adopt!(L::BVLayout)
    N = _BV_N[]; N > 0 || error("block value pass: the number of periods is not set")
    empty!(_BV_SLOTS); empty!(_BV_BUFS)
    sh = _BV_SHARED[]
    for t in 1:N
        if isnothing(sh)
            buf = Vector{Float64}(undef, L.len); push!(_BV_BUFS, buf)
            push!(_BV_SLOTS, _bv_wrap(pointer(buf), L))
        else
            length(sh[t]) >= L.len || error("block value pass: a shared segment is too short")
            push!(_BV_SLOTS, _bv_wrap(pointer(sdata(sh[t])), L))
        end
    end
    _BV_END[] = BVVal([PintVal(zeros(length(st), length(st)), zeros(length(st))) for st in L.idx])
    _BV_LAYOUT[] = L
    return L
end

function bv_layout!(b0::BVBlock0)
    L = _BV_LAYOUT[]
    if isnothing(L)
        lock(_BV_LOCK) do
            isnothing(_BV_LAYOUT[]) && _bv_adopt!(_bv_make_layout(b0.bat, b0.nP, length(b0.bat)))
        end
        L = _BV_LAYOUT[]
    end
    (L.bat == b0.bat && L.nP == b0.nP) || error("block value pass: a stage has other feeder blocks than the first one")
    return L::BVLayout
end

# The stage's function into its slot: _pp_elem!'s sums in the same order,
# entry by entry, on the blocks. nothing if a term couples two blocks.
function bv_elem!(t::Int, b0::BVBlock0, f::Vector{Float64}, rhat::Vector{Float64},
                  Rhat::SparseMatrixCSC, D::SparseMatrixCSC, dvec::Vector{Float64})
    L = bv_layout!(b0); x = _BV_SLOTS[t]::BVElem
    gof, lof = L.gof, L.lof
    p = 1 ./ f
    @inbounds for (bi, st) in enumerate(L.idx)
        H = b0.H0[L.gid[bi]]; e = x.el[bi]; m = length(st)
        for l in 1:m
            pl = p[st[l]]
            for k in 1:m; e.Mbb[k, l] = p[st[k]] * H[k, l] * pl; end
        end
        copyto!(e.Mab, e.Mbb)                                  # Np + Rp, Rp[i, j] = Rhat[j, i] p[j]
    end
    rv = rowvals(Rhat); nzv = nonzeros(Rhat)
    @inbounds for i in 1:size(Rhat, 2), q in nzrange(Rhat, i)
        j = rv[q]
        gof[i] == gof[j] || return _bv_fail("the cost or the energy rows couple batteries of two feeders")
        x.el[gof[i]].Mab[lof[i], lof[j]] += nzv[q] * p[j]
    end
    @inbounds for e in x.el; copyto!(e.Maa, e.Mab); end        # ... + Rp' + D
    @inbounds for i in 1:size(Rhat, 2), q in nzrange(Rhat, i)
        j = rv[q]
        x.el[gof[i]].Maa[lof[j], lof[i]] += nzv[q] * p[j]
    end
    rd = rowvals(D); nzd = nonzeros(D)
    @inbounds for j in 1:size(D, 2), q in nzrange(D, j)
        r = rd[q]
        gof[r] == gof[j] || return _bv_fail("the state cost couples batteries of two feeders")
        x.el[gof[j]].Maa[lof[r], lof[j]] += nzd[q]
    end
    @inbounds for (bi, st) in enumerate(L.idx)
        e = x.el[bi]
        e.Mab .*= -1
        for (k, i) in enumerate(st)
            pr = p[i] * rhat[i]
            e.ma[k] = pr + dvec[i]
            e.mb[k] = -pr
        end
    end
    return x
end

# ----------------------------------------------------------- combinations --

# Scratch for blocks of each size, per task.
function _bv_works()
    return get!(() -> Dict{Int, PairWork}(), task_local_storage(), :bv_work)::Dict{Int, PairWork}
end
_bv_work(d::Dict{Int, PairWork}, m::Int) = get!(() -> PairWork(m), d, m)

function _pp_combine!(x::BVElem, y::BVElem, ::PairWork)
    d = _bv_works()
    @inbounds for g in eachindex(x.el)
        _pp_combine!(x.el[g], y.el[g], _bv_work(d, length(x.el[g].ma)))
    end
    return x
end
function _pp_combine!(x::BVElem, y::BVVal, ::PairWork)
    d = _bv_works()
    @inbounds for g in eachindex(x.el)
        _pp_combine!(x.el[g], y.val[g], _bv_work(d, length(x.el[g].ma)))
    end
    return BVVal(x.val)
end

# The value function before every period by the ordinary recursion, last
# period first, for the feeder blocks part, part+parts, ...: the feeders are
# independent, so the parts can run at the same time (threads, or processes
# on the shared segments). Afterwards slot t holds the value function before
# period t, as after the rounds.
function bv_sweep!(N::Int, part::Int=1, parts::Int=1)
    L = _BV_LAYOUT[]::BVLayout
    d = _bv_works()
    yend = (_BV_END[]::BVVal).val
    @inbounds for g in part:parts:length(L.idx)
        w = _bv_work(d, length(L.idx[g]))
        y = yend[g]
        for t in N:-1:1
            y = _pp_combine!((_BV_SLOTS[t]::BVElem).el[g], y, w)
        end
    end
    return nothing
end
bv_sweep_here() = maximum(length, (_BV_LAYOUT[]::BVLayout).idx) <= _bv_sweep_max()
bv_blocks() = length((_BV_LAYOUT[]::BVLayout).idx)

# The sweep in this process, in `parts` parts. timed (the one-worker count):
# the parts one after another, each timed without collections; returns the sum
# and the slowest part. Otherwise on the Julia threads.
function bv_sweep_parts!(N::Int, parts::Int, timed::Bool)
    parts = max(1, min(parts, bv_blocks()))
    if timed
        s = 0.0; mx = 0.0
        for p in 1:parts
            t0 = time_ns(); c0 = Base.gc_num().total_time
            bv_sweep!(N, p, parts)
            el = (time_ns() - t0) / 1e9 - (Base.gc_num().total_time - c0) / 1e9
            s += el; mx = max(mx, el)
        end
        return s, mx
    end
    t0 = time_ns()
    @sync for p in 1:parts
        Threads.@spawn bv_sweep!(N, p, parts)
    end
    el = (time_ns() - t0) / 1e9
    return el, el
end

# --------------------------------------------- the stage's battery solution --

# Step 3 of a stage up to its battery solution, given the value function after
# it as blocks: _pp_solve with one factor per feeder. The battery rows of the
# feedback are formed here too, block by block.
function _pp_solve(S, Vn::BVVal)
    lay, stat, fac = S.lay, S.stat, S.fac
    b0 = S.b0::BVBlock0
    L = _BV_LAYOUT[]::BVLayout
    nB = length(S.active); nx = nB
    f = nonzeros(S.fu_active)                               # one entry per battery, in order (checked by _pp_elem)
    V̂x = Vector{Float64}(undef, nx)
    @inbounds for (bi, st) in enumerate(L.idx)
        v = Vn.val[bi].v
        for (k, i) in enumerate(st); V̂x[i] = v[k]; end
    end
    Qû = S.lu_ + S.fu' * V̂x + S.mubar
    ng = length(stat.gkeys)
    length(L.idx) == ng || return _bv_fail("a feeder group holds no battery")
    # The battery cross term B = f .* V (+ lux) is block diagonal: kept per block
    # (Bblk) and, for the stored policy, as a sparse matrix.
    Bblk = Vector{Matrix{Float64}}(undef, length(L.idx))
    Dfac = Vector{LU{Float64, Matrix{Float64}, Vector{Int}}}(undef, ng)
    DUp = zeros(size(b0.Up))
    nnzB = sum(length(st)^2 for st in L.idx)
    BI = Vector{Int}(undef, nnzB); BJ = Vector{Int}(undef, nnzB); BVv = Vector{Float64}(undef, nnzB); q = 0
    @inbounds for (bi, st) in enumerate(L.idx)
        gi = L.gid[bi]; V = Vn.val[bi].V; m = length(st); H = b0.H0[gi]
        Sh = Matrix{Float64}(undef, m, m); Bb = Matrix{Float64}(undef, m, m)
        for l in 1:m, k in 1:m
            bkl = f[st[k]] * V[k, l]
            Bb[k, l] = bkl
            Sh[k, l] = H[k, l] + bkl * f[st[l]]
        end
        Bblk[bi] = Bb
        F = lu!(Sh; check=false)
        issuccess(F) || return nothing
        Dfac[gi] = F
        rP = stat.granges[gi][1:m]
        DUp[rP, b0.cols] = F \ b0.Up[rP, b0.cols]
    end
    if S.lux_in_B
        lx_ = sparse(S.lux[S.active, :]); rv = rowvals(lx_); nzv = nonzeros(lx_)
        @inbounds for j in 1:size(lx_, 2), p in nzrange(lx_, j)
            i = rv[p]
            L.gof[i] == L.gof[j] || return _bv_fail("the cost couples a battery power to the energy of another feeder")
            Bblk[L.gof[i]][L.lof[i], L.lof[j]] += nzv[p]
        end
    end
    @inbounds for (bi, st) in enumerate(L.idx)
        Bb = Bblk[bi]; m = length(st)
        for l in 1:m, k in 1:m
            q += 1; BI[q] = st[k]; BJ[q] = st[l]; BVv[q] = Bb[k, l]
        end
    end
    B_active = sparse(BI, BJ, BVv, nB, nx)
    Qû = Qû + S.cuϕ
    small = Matrix{Float64}(fac.Mroot)
    isempty(b0.cols) || (small[b0.cols, b0.cols] .-= b0.Up[:, b0.cols]' * DUp[:, b0.cols])
    small = lu!(small; check=false)
    issuccess(small) || return nothing
    st = TreeKKTStructured(stat.perm, stat.granges, Dfac, b0.Up, DUp, small, fac.Mroot, zeros(0, 0), nothing,
                           b0.nP, b0.Ninv, b0.PN, b0.NP, b0.prow)
    s = TreeKKTSolver(lay, stat, fac, st, S.nzc, S.pre, S.z, S.y1, S.y2, S.t2, zeros(5))
    b = S.b; nz = S.nzc; y = S.y1
    @inbounds for j in S.active; b[j] = -Qû[j]; end
    nE = length(lay.E)
    rE = Matrix{Float64}(undef, nE, 1)
    @inbounds for a in 1:nE; rE[a, 1] = b[lay.E[a]]; end
    @inbounds for (a, t, p) in stat.kef; rE[a, 1] -= nz[p] * y[lay.Fc[t]]; end
    @inbounds for (si, i) in enumerate(lay.iso)
        w = b[i] / nz[stat.iso_d[si]]
        for (a, p) in stat.iso_e[si]; rE[a, 1] -= nz[p] * w; end
    end
    xE = tree_kkt_schur_solve(s.st, rE)
    all(isfinite, xE) || return nothing
    XB = _bv_feedback_rows(S, L, Dfac, Bblk, nB, nx)
    isnothing(XB) && return nothing
    return (; s, xE, Qû, B_active, XB)
end

# Battery rows of the feedback, (H0_g + curvature) \ (-B + S_PN inv(S_NN) cxE)
# on each feeder's block: tree_kkt_feedback_rows block by block.
function _bv_feedback_rows(S, L::BVLayout, Dfac, Bblk::Vector{Matrix{Float64}}, nB::Int, nx::Int)
    stat = S.stat; b0 = S.b0::BVBlock0; cxE = S.cxE
    nE = length(stat.perm)
    size(cxE) == (nE - nB, nx) || return _bv_fail("the energy rows have another shape than expected")
    ecol = zeros(Int, nE - nB); ecx = zeros(nE - nB)
    rv = rowvals(cxE); nzv = nonzeros(cxE)
    for j in 1:nx, q in nzrange(cxE, j)
        ecol[rv[q]] == 0 || return _bv_fail("an energy row has more than one entry")
        ecol[rv[q]] = j; ecx[rv[q]] = nzv[q]
    end
    nnzB = sum(length(st)^2 for st in L.idx)
    XI = Vector{Int}(undef, nnzB); XJ = Vector{Int}(undef, nnzB); XV = Vector{Float64}(undef, nnzB); q = 0
    @inbounds for (bi, st) in enumerate(L.idx)
        gi = L.gid[bi]; m = length(st)
        X = -Bblk[bi]
        rg = stat.granges[gi]; PNg = b0.PN[gi]; ninv = b0.Ninv[gi]
        prv = rowvals(PNg); pnz = nonzeros(PNg)
        for l in 1:length(rg) - m
            er = stat.perm[rg[m + l]] - nB; j = ecol[er]; j == 0 && continue
            L.gof[j] == bi || return _bv_fail("an energy row couples batteries of two feeders")
            c = ecx[er]
            for q in nzrange(PNg, l)
                X[prv[q], L.lof[j]] += pnz[q] * ninv[l] * c
            end
        end
        ldiv!(Dfac[gi], X)
        for l in 1:m, k in 1:m
            q += 1; XI[q] = st[k]; XJ[q] = st[l]; XV[q] = X[k, l]
        end
    end
    return sparse(XI, XJ, XV, nB, nx)
end

# FILTERDDP_VALUE_BLOCKS_AUDIT=1 (checks only, not for timed runs): the dense
# block of the same stage beside the blocks. Returns the largest entry of the
# dense block outside the feeder blocks and the largest difference inside
# them, both relative to the largest block entry.
function bv_audit(lay, stat, fac, K, Bd, Up, nB::Int, b0::BVBlock0)
    d0 = _tree_kkt_finish0(lay, stat, fac, K, Bd, Up, nB)
    isnothing(d0) && return (NaN, NaN)
    hmax = 0.0
    for H in b0.H0, v in H; hmax = max(hmax, abs(v)); end
    out = 0.0; ins = 0.0; o = 0
    grp = zeros(Int, nB)
    for (gi, np) in enumerate(b0.nP); grp[o+1:o+np] .= gi; o += np; end
    offs = cumsum(vcat(0, b0.nP))
    @inbounds for l in 1:nB, k in 1:nB
        if grp[k] == grp[l]
            g = grp[k]
            ins = max(ins, abs(d0.Sh0[k, l] - b0.H0[g][k - offs[g], l - offs[g]]))
        else
            out = max(out, abs(d0.Sh0[k, l]))
        end
    end
    return (out / hmax, ins / hmax)
end
