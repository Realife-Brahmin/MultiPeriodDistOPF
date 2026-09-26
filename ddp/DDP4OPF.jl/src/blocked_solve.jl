# Blocked multi-right-hand-side solve over UMFPACK's own factors.
#
# UMFPACK's solve accepts one right-hand side per call, so ldiv!(F, R) on the
# (n_x+1)-column stage RHS reads the whole factor once per column. Here the
# factors are taken from F (L*U = (Rs .* K)[p, q]) and the columns are solved w
# at a time: each factor entry is read once per block and applied to all w
# columns as one contiguous update. Same factors, same arithmetic, fewer passes
# over memory. See ddp/notes/KKT_ORDERING_AND_MA57.md.
#
# Opt in with FILTERDDP_BLOCKED_SOLVE=<w> (16-32 measured best). Unset or 0 keeps
# ldiv!(F, R) exactly.

function _blocked_lower!(Bt::Matrix{Float64}, L::SparseMatrixCSC{Float64, Int64})
    cp, rv, nz = L.colptr, L.rowval, L.nzval
    w = size(Bt, 1)
    @inbounds for j in 1:size(Bt, 2)
        k1, k2 = cp[j], cp[j + 1] - 1
        d = nz[k1]                      # CSC rows are sorted: the diagonal comes first
        if d != 1.0
            @simd for r in 1:w
                Bt[r, j] /= d
            end
        end
        for k in k1 + 1:k2
            i = rv[k]; v = nz[k]
            @simd for r in 1:w
                Bt[r, i] -= v * Bt[r, j]
            end
        end
    end
    return Bt
end

function _blocked_upper!(Bt::Matrix{Float64}, U::SparseMatrixCSC{Float64, Int64})
    cp, rv, nz = U.colptr, U.rowval, U.nzval
    w = size(Bt, 1)
    @inbounds for j in size(Bt, 2):-1:1
        k1, k2 = cp[j], cp[j + 1] - 1
        d = nz[k2]                      # ...and in U it comes last
        @simd for r in 1:w
            Bt[r, j] /= d
        end
        for k in k1:k2 - 1
            i = rv[k]; v = nz[k]
            @simd for r in 1:w
                Bt[r, i] -= v * Bt[r, j]
            end
        end
    end
    return Bt
end

function _solve_column_block!(R, work, cols, L, U, p, q, Rs)
    n = size(R, 1)
    Bb = length(cols) == size(work, 1) ? work : Matrix{Float64}(undef, length(cols), n)
    @inbounds for i in 1:n
        pi_ = p[i]; s = Rs[pi_]
        for (r, c) in enumerate(cols)
            Bb[r, i] = s * R[pi_, c]
        end
    end
    _blocked_lower!(Bb, L)
    _blocked_upper!(Bb, U)
    @inbounds for i in 1:n
        qi = q[i]
        for (r, c) in enumerate(cols)
            R[qi, c] = Bb[r, i]
        end
    end
    return R
end

# Overwrites R with K \ R. Safe in place: each block of columns is fully read
# into its work block before any of those columns is written back. With more
# than one Julia thread (JULIA_NUM_THREADS), the column blocks are split across
# threads: blocks touch disjoint columns of R and only read the factors.
function _blocked_umfpack_solve!(F, R::AbstractMatrix{Float64}, w::Int)
    L = F.L; U = F.U; p = F.p; q = F.q; Rs = F.Rs
    n, m = size(R)
    starts = collect(1:w:m)
    nt = min(Threads.nthreads(), length(starts))
    if nt <= 1
        work = Matrix{Float64}(undef, min(w, m), n)
        for c0 in starts
            _solve_column_block!(R, work, c0:min(c0 + w - 1, m), L, U, p, q, Rs)
        end
    else
        # One buffer per task. Reusing the name `work` from the branch above
        # made it a single captured variable shared by every thread (garbage
        # results, ~1e130); `local` and a distinct name keep it per task.
        Threads.@threads for k in 1:nt
            local thread_work = Matrix{Float64}(undef, min(w, m), n)
            for c0 in starts[k:nt:end]
                _solve_column_block!(R, thread_work, c0:min(c0 + w - 1, m), L, U, p, q, Rs)
            end
        end
    end
    return R
end

function _kkt_wide_solve!(F, R)
    w = parse(Int, get(ENV, "FILTERDDP_BLOCKED_SOLVE", "0"))
    (w > 0 && F isa SparseArrays.UMFPACK.UmfpackLU{Float64} && eltype(R) == Float64) ||
        return ldiv!(F, R)
    return _blocked_umfpack_solve!(F, R, w)
end
