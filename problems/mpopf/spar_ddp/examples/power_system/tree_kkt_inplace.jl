# FILTERDDP_INPLACE_FACTOR=1: the network factorization of a stage written into
# the storage of that stage's previous factorization.
#
# Julia's collector takes 17-22% of a FilterDDP solve on med2522 and large10k
# (FILTERDDP_GC_TIME, 2026-10-11). tree_kkt_factor allocates, for every node of
# the feeder tree, at every stage of every backward pass: the node's block, its
# pivot vector, the coupling block, the solved coupling block and the
# contribution to the parent, plus two temporaries. Everything but the
# temporaries is kept until the stage's next factorization, so the collector
# has several objects per bus and period to trace, all replaced once per pass.
#
# Here the storage of a stage's factor is kept and overwritten: same
# arithmetic in the same order (smalllu!, the solve with the adjoint coupling
# block and the product of tree_kkt_factor, written out in place), so the
# factor is the same to the last bit. A stage is recognised by its K, which
# the KKT pattern cache holds under (solver, stage). The old solver object of
# the stage shares the overwritten factor and is invalid from then on; the
# backward pass replaces it before anything can use it.
#
# Included at the end of tree_kkt.jl.

_inplace_factor() = get(ENV, "FILTERDDP_INPLACE_FACTOR", "0") != "0"

struct FactorSlot
    lay::Any
    fac::TreeKKTFactor
    C::Vector{Matrix{Float64}}                # contribution of each node to its parent
end

# One slot per (solver, stage): the key under which the KKT pattern cache holds
# the stage's K. A new solve gets the same keys back, so nothing accumulates.
const _FACTOR_SLOTS = Dict{Tuple{UInt, Int}, FactorSlot}()
const _FACTOR_SLOTS_LOCK = ReentrantLock()

function _factor_key(K::SparseMatrixCSC)
    for (key, entry) in SPARDDP._KKT_PATTERN_CACHE
        entry.K === K && return key
    end
    return nothing
end

# smalllu!, with the pivot vector reset instead of created.
function _smalllu_inplace!(A::Matrix{Float64}, perm::Vector{Int})
    n = size(A, 1)
    @inbounds for i in 1:n; perm[i] = i; end
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
    return nothing
end

# X = F \ B' for B (q x n), written into X (n x q): the loops of \(::SmallLU, ::AbstractMatrix).
function _smallsolve_adjoint!(X::Matrix{Float64}, A::Matrix{Float64}, perm::Vector{Int}, B::Matrix{Float64})
    n = size(A, 1); w = size(B, 1)
    @inbounds for c in 1:w
        for i in 1:n; X[i, c] = B[c, perm[i]]; end
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

# C = -(A * B), written into C: the loops of smallmul, then the sign.
function _smallmul_neg!(C::Matrix{Float64}, A::Matrix{Float64}, B::Matrix{Float64})
    m, k = size(A); w = size(B, 2)
    fill!(C, 0.0)
    @inbounds for c in 1:w, kk in 1:k
        b = B[kk, c]
        b == 0 && continue
        for i in 1:m; C[i, c] += A[i, kk] * b; end
    end
    C .*= -1
    return C
end

function _tree_kkt_factor!(slot::FactorSlot, lay::TreeKKTLayout, K::SparseMatrixCSC)
    fac, C = slot.fac, slot.C
    nz = nonzeros(K)
    @inbounds for j in lay.post
        Fj = fac.F[j]; A = Fj.A
        fill!(A, 0.0)
        for (r, c, p) in lay.oo[j]; A[r, c] = nz[p]; end
        for c in lay.children[j]
            pos = lay.iface_pos[c]; Cc = C[c]
            for b in eachindex(pos), a in eachindex(pos)
                A[pos[a], pos[b]] += Cc[a, b]
            end
        end
        j == lay.root && copyto!(fac.Mroot, A)
        _smalllu_inplace!(A, Fj.perm)
        if j != lay.root
            Bio = fac.Aio[j]
            fill!(Bio, 0.0)
            for (c, r, p) in lay.oi[j]; Bio[r, c] = nz[p]; end
            _smallsolve_adjoint!(fac.W[j], A, Fj.perm, Bio)
            _smallmul_neg!(C[j], Bio, fac.W[j])
        end
    end
    return fac
end

# tree_kkt_factor, into the stage's own storage when it has one. Without the
# KKT pattern cache a stage cannot be recognised and the factor is allocated.
function tree_kkt_factor_pooled(lay::TreeKKTLayout, K::SparseMatrixCSC)
    key = _factor_key(K)
    isnothing(key) && return tree_kkt_factor(lay, K)
    slot = lock(() -> get(_FACTOR_SLOTS, key, nothing), _FACTOR_SLOTS_LOCK)
    (!isnothing(slot) && slot.lay === lay) && return _tree_kkt_factor!(slot, lay, K)
    fac = tree_kkt_factor(lay, K)
    N = length(lay.own)
    C = Vector{Matrix{Float64}}(undef, N)
    for j in 1:N
        q = length(lay.iface[j])
        C[j] = zeros(q, q)
    end
    lock(() -> (_FACTOR_SLOTS[key] = FactorSlot(lay, fac, C)), _FACTOR_SLOTS_LOCK)
    return fac
end
