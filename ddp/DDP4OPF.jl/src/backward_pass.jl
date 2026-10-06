_derivative_matrix(A, m, n) = issparse(A) ? sparse(A) : reshape(vec(A), m, n)

# ---------------------------------------------------------------- frozen KKT --
# Fast-decoupled-style reuse of the per-stage KKT factorisation. FDPF freezes the
# Jacobian because susceptance is genuinely constant; here K = [H cu'; cu 0] is
# NOT constant -- it carries the barrier terms Sigma_L + Sigma_U, which scale
# like mu/s^2 and move sharply as mu drops and iterates approach bounds. So this
# is a quasi-Newton approximation whose worth is an empirical question: it saves
# assembly+factorisation per stage per iteration, and pays for it in step
# quality. Opt in with FILTERDDP_FREEZE_KKT=N (refactor every N backward passes);
# N <= 1 is the default and reproduces exact behaviour bit for bit.
const _FROZEN_KKT = Dict{Int,Any}()
const _FROZEN_AT  = Dict{Int,Int}()
const _FROZEN_STATS = Dict{Symbol,Int}(:factorisations => 0, :reuses => 0)
const _KKT_PATTERN_CACHE = Dict{Tuple{UInt,Int},Any}()

_freeze_period() = parse(Int, get(ENV, "FILTERDDP_FREEZE_KKT", "1"))

function _frozen_reset!()
    empty!(_FROZEN_KKT); empty!(_FROZEN_AT)
    empty!(_KKT_PATTERN_CACHE)
    empty!(_LAGGED_CURV)
    _FROZEN_STATS[:factorisations] = 0; _FROZEN_STATS[:reuses] = 0
    return nothing
end

# Opt-in UMFPACK strategy/ordering, for the ordering experiments in
# ddp/notes/KKT_ORDERING_AND_MA57.md. Values are UMFPACK's own codes, e.g.
# FILTERDDP_UMFPACK_STRATEGY=3 (symmetric), FILTERDDP_UMFPACK_ORDERING=3 (METIS).
# Unset (the default) calls lu(K) exactly as before.
function _kkt_lu(K)
    strategy = get(ENV, "FILTERDDP_UMFPACK_STRATEGY", "")
    ordering = get(ENV, "FILTERDDP_UMFPACK_ORDERING", "")
    isempty(strategy) && isempty(ordering) && return lu(K)
    control = SparseArrays.UMFPACK.get_umfpack_control(Float64, Int64)
    isempty(strategy) || (control[SparseArrays.LibSuiteSparse.UMFPACK_STRATEGY + 1] = parse(Float64, strategy))
    isempty(ordering) || (control[SparseArrays.LibSuiteSparse.UMFPACK_ORDERING + 1] = parse(Float64, ordering))
    return lu(K; control=control)
end

function _frozen_lu(t::Int, K, iter::Int)
    period = _freeze_period()
    if period <= 1
        _FROZEN_STATS[:factorisations] += 1
        return _kkt_lu(K)
    end
    if !haskey(_FROZEN_KKT, t) || (iter - _FROZEN_AT[t]) >= period
        _FROZEN_KKT[t] = _kkt_lu(K)
        _FROZEN_AT[t] = iter
        _FROZEN_STATS[:factorisations] += 1
    else
        _FROZEN_STATS[:reuses] += 1
    end
    return _FROZEN_KKT[t]
end


# ----------------------------------------------------- diagonal Hessian --
# Agenda question (R. Gupta, 2026-09-18): "can we not just use a diagonalized
# matrix?". FILTERDDP_DIAG_HESSIAN=1 replaces the stage Hessian block of
# K = [H cu'; cu 0] with its diagonal, leaving the constraint Jacobian intact.
#
# Two things make this a real approximation rather than a cheap trick:
#   * On ieee123 T=3 about 83% of nnz(H) IS the dense nB x nB block
#     fu' * Vxx * fu -- the intertemporal curvature. Diagonalising deletes
#     exactly the information that makes this a second-order method.
#   * diag(H) is singular in general: at the terminal stage nnz(H) = 536
#     against nu = 791, because line flows, voltages and currents pick up a
#     diagonal entry only through cuu, and unbounded controls get no barrier
#     term at all. Hence FILTERDDP_DIAG_HESSIAN_FLOOR (default 1e-8), which
#     floors each diagonal entry from below.
_diag_hessian_enabled() = get(ENV, "FILTERDDP_DIAG_HESSIAN", "0") != "0"
_diag_hessian_floor() = parse(Float64, get(ENV, "FILTERDDP_DIAG_HESSIAN_FLOOR", "1e-8"))
_direct_diag_hessian_enabled() = get(ENV, "FILTERDDP_DIRECT_DIAG_HESSIAN", "0") != "0"

function _diagonalise_hessian(H)
    floor_val = _diag_hessian_floor()
    d = diag(H)
    @inbounds for i in eachindex(d)
        d[i] = max(d[i], floor_val)
    end
    return issparse(H) ? spdiagm(0 => d) : diagm(d)
end

_kkt_pattern_cache_enabled() = get(ENV, "FILTERDDP_CACHE_KKT_PATTERN", "0") != "0"

# Diagnostic (agenda of 2026-10-02: "can the KKT refresh be pinpointed?").
# With FILTERDDP_STALE_JACOBIAN_PERIOD = p > 1, the KKT matrix uses each
# stage's constraint Jacobian from the last iteration that was a multiple of p;
# the Hessian, barrier terms, gradient (cu' * phi) and residuals stay current,
# so only the matrix is approximate (an inexact Newton step). Every linear
# constraint row is constant anyway, so this makes exactly the SOCP rows stale.
# Off (p <= 1) by default: the Jacobian is used as is.
const _STALE_CU = Dict{Int, Any}()
function _stale_jacobian(t::Int, k::Int, cu)
    p = parse(Int, get(ENV, "FILTERDDP_STALE_JACOBIAN_PERIOD", "0"))
    p <= 1 && return cu
    if k == 0 || k % p == 0 || !haskey(_STALE_CU, t)
        _STALE_CU[t] = copy(cu)
    end
    return _STALE_CU[t]
end

# ------------------------------------------------ lagged curvature --
# Parallel-in-time experiment (agenda of 2026-10-07). In the diagonal-Hessian
# stage KKT the ONLY entries that depend on stage t+1 are the battery-power
# diagonal curvature diag(fu' * Vxx_{t+1} * fu) (_diagonal_quadratic_form
# below); every other entry of K_t is a function of the current iterate at
# stage t. With FILTERDDP_LAGGED_CURVATURE=1 those entries are taken from the
# previous backward sweep, so every K_t is known before the sweep starts and
# all T factorizations run in parallel on Julia threads (BLAS pinned to one
# thread meanwhile). Right-hand sides, gains and the value recursion stay
# exact; only the matrix is approximate, an inexact Newton step like the
# diagonal Hessian itself. The first sweep of a solve has nothing to lag and
# runs sequentially. Requires the direct diagonal path and dynamics with
# fuu = 0 (checked).
_lagged_curvature_enabled() = get(ENV, "FILTERDDP_LAGGED_CURVATURE", "0") != "0" &&
    _diag_hessian_enabled() && _direct_diag_hessian_enabled()
const _LAGGED_CURV = Dict{Int, Vector{Float64}}()

function _lagged_prefactor(ocp, traj, reg, nx::Int, nu::Int, nc::Int)
    N = ocp.N
    all(t -> haskey(_LAGGED_CURV, t), 1:N) || return nothing
    cl = ocp.control_limits
    dynamics = ocp.dynamics
    factors = Vector{Any}(nothing, N)
    blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        Threads.@threads :dynamic for t in 1:N
            factors[t] = try
                x, u, ϕ, zl, zu = traj[t].x, traj[t].u, traj[t].ϕ, traj[t].zl, traj[t].zu
                objective_t = t == N ? ocp.term_objective : stage_obj(ocp, t)
                constraints = stage_con(ocp, t)
                luu = _derivative_matrix(objective_t.luu(x, u), nu, nu)
                cu = _derivative_matrix(constraints.cu(x, u), nc, nu)
                cuu = _derivative_matrix(constraints.cuu(x, u, ϕ), nu, nu)
                fuu = _derivative_matrix(dynamics.fuu(x, u, zeros(nx)), nu, nu)
                nnz(sparse(fuu)) == 0 || error("FILTERDDP_LAGGED_CURVATURE needs fuu = 0")
                # Same operations, same order, as the stage loop below.
                ul = u - cl.l
                uu = cl.u - u
                inv_ul = inv.(ul) .* cl.maskl
                inv_uu = inv.(uu) .* cl.masku
                Σ_L = inv_ul .* zl
                Σ_U = inv_uu .* zu
                Ĥ = sparse(luu) + spdiagm(0 => Σ_L + Σ_U + _LAGGED_CURV[t]) + sparse(fuu)
                Ĥ = Ĥ + cuu
                if !iszero(reg)
                    @inbounds for i in axes(Ĥ, 1)
                        Ĥ[i, i] += reg
                    end
                end
                Ĥ = sparse(Symmetric(_diagonalise_hessian(Ĥ)))
                cu_sparse = sparse(cu)
                _kkt_lu([Ĥ sparse(cu_sparse'); cu_sparse spzeros(eltype(Ĥ), nc, nc)])
            catch err
                err
            end
        end
    finally
        BLAS.set_num_threads(blas_threads)
    end
    return factors
end

# ------------------------------------------------ battery-block value update --
# Lead of 2026-09-30 (ddp/notes/PARALLEL_IN_TIME.md, Section 4). With the
# factor-backed policy the n_x feedback columns of K \ rhs are used only in
# the rows E = {battery powers, energy constraints}, which are also the only
# nonzero rows of those right-hand-side columns, and beta' Qu + omega' c equals
# B' alpha + cx' psi by symmetry of K. So the value update needs only
# (K^{-1})_EE = S^{-1}, S the Schur complement of K onto E (2 n_B x 2 n_B),
# and no n_x-column solve. FILTERDDP_BATTERY_SCHUR=1 computes it through
# BATTERY_SCHUR_HOOK[] (K, E) -> S, which the driver sets (MUMPS, so that this
# package takes no new dependency). Exact up to rounding; the stage
# factorization is still used for the feedforward column and the policy.
const BATTERY_SCHUR_HOOK = Ref{Any}(nothing)
_battery_schur_enabled() = get(ENV, "FILTERDDP_BATTERY_SCHUR", "0") != "0" &&
    !isnothing(BATTERY_SCHUR_HOOK[])

# FILTERDDP_TREE_KKT=1: TREE_KKT_HOOK[] (K) -> F replaces the stage's sparse
# LU by a solver that exploits the radial network (tree_kkt.jl in the driver:
# leaves-to-substation block elimination, battery block by feeder plus a
# low-rank substation term). F must support ldiv!(F, b) for the feedforward
# column and the forward-pass policy, and battery_block_rows for the feedback
# rows, so no second factorization is needed. Same battery-block value update
# as above; the terminal stage (l_ux != 0) keeps the sparse LU.
const TREE_KKT_HOOK = Ref{Any}(nothing)
_structured_dynamics() = get(ENV, "FILTERDDP_STRUCTURED_DYNAMICS", "0") != "0"

# Exact stage Hessian with the tree solver. The battery curvature
# fu' * Vxx * fu is dense, but it lies entirely in the P_B x P_B block, which
# the tree solver keeps (it eliminates the network onto the battery rows). So
# with FILTERDDP_TREE_KKT=1 and without FILTERDDP_DIAG_HESSIAN it is not
# assembled into K: K carries the rest of the exact Hessian (cost, barrier and
# constraint curvature, off-diagonal terms included) and the dense block is
# handed to the solver, which adds it to the battery Schur complement. Exact,
# and the network elimination is the same as for the diagonal Hessian.
#
# The Hessian is built with a pattern that does not depend on the values
# (explicit zeros kept), so the KKT pattern cache and the tree layout stay
# valid across iterations.
function _fixed_pattern_hessian(nu::Int, sigma::AbstractVector, blocks...)
    I = collect(1:nu); J = collect(1:nu); V = Vector{Float64}(sigma)
    for M in blocks
        i, j, v = findnz(sparse(M))
        append!(I, i); append!(J, j); append!(V, v)
    end
    return sparse(I, J, V, nu, nu)
end
_tree_kkt_enabled() = get(ENV, "FILTERDDP_TREE_KKT", "0") != "0" && !isnothing(TREE_KKT_HOOK[])

# FILTERDDP_PARSIM=1: accounting for one worker per period (with the tree
# solver). Nothing is reordered and no result changes; each stage's backward
# work is timed in three parts, by what it needs from stage t+1:
#   pre   nothing: derivatives, barrier terms, KKT assembly, network
#         factorization and the Schur preparation, the first network solve;
#   seq   the value function of t+1: its products, the battery block
#         (assembly, factorization, feedback rows), the value update;
#   post  this stage's battery solution only: the second network solve, the
#         multiplier updates and the stored policy.
# One worker per period would take max(pre) + sum(seq) + max(post) for the
# sweep. The value function reaches a stage matrix only through the battery
# block (tree_kkt.jl, _tree_kkt_prepare), which is what makes `pre` free of
# it; ddp/examples/power_system/tree_kkt_phase_check.jl checks that.
_parsim() = get(ENV, "FILTERDDP_PARSIM", "0") != "0"
stage_phase_times(F) = nothing
const _PARSIM_BETA_B = Ref{Any}(nothing)        # battery rows of the feedback, last stage solved
const _PARSIM_SKIP_NS = Ref{UInt64}(0)          # accounting-only work, left out of the times
const _PARSIM_FWD = zeros(7)                    # forward pass: see affine_linesearch.jl
# Seconds of the stage being solved, parts of `seq`: [1] Q_u and C, [2] curvature
# diagonal, [3] B, [4] right-hand side of the battery rows, [5] battery rows,
# [6] value increments, [7] value update; and [8] the feedforward solve.
const _PARSIM_SUB = zeros(8)

# FILTERDDP_LEAN_VALUE=1: the value-function part of a stage, the part that
# cannot be done before the sweep, without its avoidable dense work. With
# the battery dynamics f_x = I and f_u a scaled selector, so
#   C  = l_xx + V_xx + f_xx          (no products with f_x),
#   B  = rows of V_xx scaled        (no sparse-dense products),
# and the feedback rows are solved in the tree solver's grouped order, from
# a right-hand side written there directly (battery_block_rows_grouped):
# four 2 n_B x n_x copies fewer per stage, and the feeder blocks solved on
# the Julia threads. Every entry is computed by the same operations as
# before, so the iterates are unchanged.
_lean_value() = get(ENV, "FILTERDDP_LEAN_VALUE", "0") != "0"
battery_block_rows_grouped(F, E, B, cxE) = nothing
# FILTERDDP_VXX_SEPARABILITY=1 (diagnostic, with the lean path): how much of a
# stage's value-Hessian increment couples batteries on different feeders.
battery_groups(F, nB) = nothing
function _vxx_separability(F, Vinc, nB)
    grp = battery_groups(F, nB)
    (isnothing(grp) || size(Vinc, 1) != nB) && return nothing
    w = 0.0; c = 0.0; mw = 0.0; mc = 0.0
    @inbounds for j in 1:nB, i in 1:nB
        v = abs(Vinc[i, j])
        if grp[i] == grp[j]
            w += v^2; mw = max(mw, v)
        else
            c += v^2; mc = max(mc, v)
        end
    end
    @printf("FILTERDDP_VXX_SEP groups=%d within_fro=%.6e cross_fro=%.6e within_max=%.6e cross_max=%.6e
",
            maximum(grp), sqrt(w), sqrt(c), mw, mc)
    return nothing
end
function _is_identity(A)
    A isa SparseMatrixCSC || return false
    n = size(A, 1)
    (size(A, 2) == n && nnz(A) == n) || return false
    rv = rowvals(A); nzv = nonzeros(A)
    @inbounds for j in 1:n
        q = A.colptr[j]
        (A.colptr[j+1] == q + 1 && rv[q] == j && nzv[q] == 1) || return false
    end
    return true
end
# (lxx + V) + fxx
function _lean_C(lxx, V::Matrix, fxx)
    lxx isa SparseMatrixCSC || return lxx + V + fxx
    C = copy(V)
    rv = rowvals(lxx); nzv = nonzeros(lxx)
    @inbounds for j in axes(lxx, 2), q in nzrange(lxx, j)
        C[rv[q], j] = nzv[q] + C[rv[q], j]
    end
    if fxx isa SparseMatrixCSC
        rv = rowvals(fxx); nzv = nonzeros(fxx)
        @inbounds for j in axes(fxx, 2), q in nzrange(fxx, j)
            C[rv[q], j] += nzv[q]
        end
    else
        C .+= fxx
    end
    return C
end
# A[i, j] += X[row[j], i] * val[j] for the columns j that have a row, in
# cache-sized blocks: X is read down its columns and A written down its own.
function _add_scaled_rows_transposed!(A::Matrix, X::Matrix, row::Vector{Int}, val::Vector)
    n = size(A, 1); m = size(A, 2); bs = 64
    @inbounds for j0 in 1:bs:m, i0 in 1:bs:n
        for j in j0:min(j0 + bs - 1, m)
            r = row[j]; r == 0 && continue
            v = val[j]
            for i in i0:min(i0 + bs - 1, n)
                A[i, j] += X[r, i] * v
            end
        end
    end
    return A
end
# fa' * V for fa with exactly one entry per column
function _lean_B(fa::SparseMatrixCSC, V::Matrix)
    nB = size(fa, 2); nx = size(V, 2)
    all(k -> fa.colptr[k+1] == fa.colptr[k] + 1, 1:nB) || return Matrix(fa' * V)
    rv = rowvals(fa); nzv = nonzeros(fa)
    B = Matrix{eltype(V)}(undef, nB, nx)
    @inbounds for j in 1:nx, k in 1:nB
        B[k, j] = nzv[k] * V[rv[k], j]
    end
    return B
end
# Diagonal of a sparse matrix as a dense vector, zero where no entry is stored
function _dense_diag(A::SparseMatrixCSC)
    n = min(size(A, 1), size(A, 2))
    d = zeros(eltype(A), n)
    rv = rowvals(A); nzv = nonzeros(A)
    @inbounds for j in 1:n, q in nzrange(A, j)
        rv[q] == j && (d[j] = nzv[q])
    end
    return d
end
_dense_diag(A) = Vector(diag(A))
# The sparse diagonal matrix _diagonalise_hessian returns, from the dense diagonal
function _floored_diagonal(d::Vector)
    floor_val = _diag_hessian_floor()
    n = length(d)
    v = Vector{eltype(d)}(undef, n)
    @inbounds for i in 1:n
        v[i] = max(d[i], floor_val)
    end
    return SparseMatrixCSC(n, n, collect(1:n+1), collect(1:n), v)
end

# Rows E of K \ R for right-hand sides R supported on E (given as R[E, :]).
# Default: dense solve with the Schur complement from BATTERY_SCHUR_HOOK.
battery_block_rows(F, K, E, R) = lu!(BATTERY_SCHUR_HOOK[](K, E)) \ R

function _battery_schur_value(K, F, rhs, nu::Int, active_B_rows, B_active, cx)
    t_a = time_ns()
    ldiv!(F, @view(rhs[:, 1:1]))
    α = copy(@view rhs[1:nu, 1])
    ψ = copy(@view rhs[nu+1:end, 1])
    t_b = time_ns()
    cx_s = sparse(cx)
    energy_rows = sort!(unique(cx_s.rowval))
    E = vcat(active_B_rows, nu .+ energy_rows)
    cxE = cx_s[energy_rows, :]                       # sparse: one entry per battery
    nB = length(active_B_rows)
    t_g0 = time_ns()
    grouped = _lean_value() ? battery_block_rows_grouped(F, E, B_active, cxE) : nothing
    t_g1 = time_ns()
    if isnothing(grouped)
        R = vcat(-B_active, -Matrix(cxE))
        t_c = time_ns()
        X = battery_block_rows(F, K, E, R)
        t_d = time_ns()
        if _parsim()
            t0 = time_ns(); _PARSIM_BETA_B[] = X[1:nB, :]; _PARSIM_SKIP_NS[] += time_ns() - t0
        end
        Vxx_inc = (@view X[1:nB, :])' * B_active
        Vxx_inc .+= X[nB+1:end, :]' * cxE
    else
        # Xp: the same rows in the solver's grouped order, row invp[i] holding row i of E
        Xp, invp = grouped
        t_c = t_g0; t_d = t_g1                # the right-hand side is written inside the solve
        nx = size(B_active, 2)
        XB = Matrix{eltype(Xp)}(undef, nB, nx)
        nEn = size(cxE, 1)
        XE = Matrix{eltype(Xp)}(undef, nEn, nx)
        @inbounds for j in 1:nx
            for k in 1:nB; XB[k, j] = Xp[invp[k], j]; end
            for k in 1:nEn; XE[k, j] = Xp[invp[nB + k], j]; end
        end
        _parsim() && (_PARSIM_BETA_B[] = XB)
        Vxx_inc = XB' * B_active
        # + XE' * cxE, one entry of cxE per column
        rv = rowvals(cxE); nzv = nonzeros(cxE)
        if all(j -> cxE.colptr[j+1] - cxE.colptr[j] <= 1, 1:nx)
            row = zeros(Int, nx); val = zeros(eltype(nzv), nx)
            @inbounds for j in 1:nx, q in nzrange(cxE, j)
                row[j] = rv[q]; val[j] = nzv[q]
            end
            _add_scaled_rows_transposed!(Vxx_inc, XE, row, val)
        else
            Vxx_inc .+= XE' * cxE
        end
        get(ENV, "FILTERDDP_VXX_SEPARABILITY", "0") != "0" && _vxx_separability(F, Vxx_inc, nB)
    end
    Vx_inc = B_active' * α[active_B_rows] + cx_s' * ψ
    _PARSIM_SUB[8] = (t_b - t_a) / 1e9; _PARSIM_SUB[4] = (t_c - t_b) / 1e9
    _PARSIM_SUB[5] = (t_d - t_c) / 1e9; _PARSIM_SUB[6] = (time_ns() - t_d) / 1e9
    return α, ψ, Vxx_inc, Vx_inc
end

function _stage_factor(t::Int, K, iter::Int, lagged)
    isnothing(lagged) && return _frozen_lu(t, K, iter)
    F = lagged[t]
    F isa Exception && throw(F)
    return F
end

function _same_sparse_pattern(A, colptr, rowval)
    return A.colptr == colptr && A.rowval == rowval
end

function _cached_kkt!(key::Tuple{UInt,Int}, H::SparseMatrixCSC{T,Int},
                      cu::SparseMatrixCSC{T,Int}, nu::Int, nc::Int) where {T}
    entry = get(_KKT_PATTERN_CACHE, key, nothing)
    if isnothing(entry) || !_same_sparse_pattern(H, entry.H_colptr, entry.H_rowval) ||
            !_same_sparse_pattern(cu, entry.cu_colptr, entry.cu_rowval)
        K = [H sparse(cu'); cu spzeros(T, nc, nc)]
        locations = Dict{Tuple{Int,Int},Int}()
        for col in axes(K, 2), p in nzrange(K, col)
            locations[(K.rowval[p], col)] = p
        end
        hmap = Vector{Int}(undef, nnz(H))
        for col in axes(H, 2), p in nzrange(H, col)
            hmap[p] = locations[(H.rowval[p], col)]
        end
        lower = Vector{Int}(undef, nnz(cu))
        upper = similar(lower)
        for col in axes(cu, 2), p in nzrange(cu, col)
            row = cu.rowval[p]
            lower[p] = locations[(nu + row, col)]
            upper[p] = locations[(col, nu + row)]
        end
        entry = (K=K, H_colptr=copy(H.colptr), H_rowval=copy(H.rowval),
                 cu_colptr=copy(cu.colptr), cu_rowval=copy(cu.rowval),
                 hmap=hmap, lower=lower, upper=upper)
        _KKT_PATTERN_CACHE[key] = entry
    end
    K = entry.K
    fill!(K.nzval, zero(T))
    @views K.nzval[entry.hmap] .= H.nzval
    @views K.nzval[entry.lower] .= cu.nzval
    @views K.nzval[entry.upper] .= cu.nzval
    return K
end

# Compute diag(A' * M * A) without materialising that full product.  Network
# dynamics touch only the battery-power columns of A, so this avoids building
# and immediately discarding the dense battery block in diagonal-Hessian mode.
function _diagonal_quadratic_form(A::SparseMatrixCSC{T}, M, n::Int) where {T}
    d = zeros(promote_type(T, eltype(M)), n)
    @inbounds for j in axes(A, 2)
        lo = A.colptr[j]
        hi = A.colptr[j + 1] - 1
        lo > hi && continue
        rows = @view A.rowval[lo:hi]
        vals = @view A.nzval[lo:hi]
        acc = zero(eltype(d))
        for a in eachindex(rows), b in eachindex(rows)
            acc += vals[a] * M[rows[a], rows[b]] * vals[b]
        end
        d[j] = acc
    end
    return d
end


function backward_pass!(solver::Solver{T, nx, nu, nc, nux, ncx}, ocp::OCP{T, nx, nu, nc},
            traj::Vector{TrajectoryElement{T, nx, nu, nc}}, data::SolverData{T}, options::Options{T}; verbose::Bool=false
            ) where {T, nx, nu, nc, nux, ncx}
    # A stale factorisation must never leak between solves.
    data.k == 0 && _frozen_reset!()
    reg::T = 0.0
    μ = data.μ
    δ_c = 0.
    reg = 0.0

    dynamics = ocp.dynamics
    cl = ocp.control_limits
    ni = (cl.nl + cl.nu) * ocp.N
    rhs_workspace = nothing
    feasibility_diagnostic = get(ENV, "FILTERDDP_FEASIBILITY_DIAGNOSTIC", "0") == "1"
    equality_ssq = T(0); equality_count = 0; equality_max = T(0)
    equality_worst_stage = 0; equality_worst_index = 0
    dynamics_ssq = T(0); dynamics_count = 0; dynamics_max = T(0)
    dynamics_worst_stage = 0; dynamics_worst_index = 0
    bound_ssq = T(0); bound_count = 0; bound_max = T(0)
    bound_worst_stage = 0; bound_worst_index = 0; bound_worst_kind = "none"
    
    parsim = _parsim()
    ps_tot = zeros(5)                    # pre_sum, sum of per-attempt pre_max, seq_sum, post_sum, sum of post_max
    ps_stages = 0
    while reg <= options.reg_max
        equality_ssq = T(0); equality_count = 0; equality_max = T(0)
        equality_worst_stage = 0; equality_worst_index = 0
        dynamics_ssq = T(0); dynamics_count = 0; dynamics_max = T(0)
        dynamics_worst_stage = 0; dynamics_worst_index = 0
        bound_ssq = T(0); bound_count = 0; bound_max = T(0)
        bound_worst_stage = 0; bound_worst_index = 0; bound_worst_kind = "none"
        data.status = 0
        V̂x = zeros(T, nx)
        V̂xx = zeros(T, nx, nx)
        λ = zeros(T, nx)

        data.barrier_lagrangian_curr = T(0.0)
        data.primal_1_curr = T(0.0)
        data.primal_inf = T(0.0)
        data.cs_inf_μ = T(0.0)
        data.cs_inf_0 = T(0.0)
        data.objective = T(0.0)
        data.dual_inf = T(0.0)
        data.expected_change_L = T(0.0)
        ϕ_norm = T(0.0)
        z_norm = T(0.0)
        ps_att = zeros(5)                # pre_sum, pre_max, seq_sum, post_sum, post_max

        # FILTERDDP_LAGGED_CURVATURE: factor every stage up front, in parallel.
        lagged_factors = nothing
        if _lagged_curvature_enabled()
            prefactor_start_ns = time_ns()
            lagged_factors = _lagged_prefactor(ocp, traj, reg, nx, nu, nc)
            get(ENV, "FILTERDDP_TIMING_DIAGNOSTIC", "0") == "1" && @printf(
                "FILTERDDP_PREFACTOR iteration=%d barrier_iteration=%d lagged=%d wall_s=%.9f threads=%d\n",
                data.k, data.j, !isnothing(lagged_factors),
                (time_ns() - prefactor_start_ns) / 1e9, Threads.nthreads())
        end

        for t = ocp.N:-1:1
            timing_diagnostic = get(ENV, "FILTERDDP_TIMING_DIAGNOSTIC", "0") == "1"
            memory_diagnostic = get(ENV, "FILTERDDP_MEMORY_DIAGNOSTIC", "0") == "1"
            stage_alloc_start = memory_diagnostic ? Base.gc_bytes() : 0
            stage_maxrss_start = memory_diagnostic ? Sys.maxrss() : 0
            derivative_start_ns = time_ns()
            ps_seq_ns = UInt64(0); ps_val_ns = UInt64(0); _PARSIM_SKIP_NS[] = 0; fill!(_PARSIM_SUB, 0.0)
            x, u, ϕ, zl, zu = traj[t].x, traj[t].u, traj[t].ϕ, traj[t].zl, traj[t].zu
            # PATCH: per-stage data when supplied, else the shared function
            constraints = stage_con(ocp, t)
            objective_t = stage_obj(ocp, t)

            # evaluate derivatives

            first_order_ns = 0
            second_order_ns = 0
            dynamics_first_order_ns = 0
            if t == ocp.N
                callback_start_ns = time_ns()
                lx = Vector(ocp.term_objective.lx(x, u))
                lu_ = Vector(ocp.term_objective.lu(x, u))
                data.objective += ocp.term_objective.l(x, u)[1]
                first_order_ns += time_ns() - callback_start_ns
                callback_start_ns = time_ns()
                lxx = _derivative_matrix(ocp.term_objective.lxx(x, u), nx, nx)
                lux = _derivative_matrix(ocp.term_objective.lux(x, u), nu, nx)
                luu = _derivative_matrix(ocp.term_objective.luu(x, u), nu, nu)
                second_order_ns += time_ns() - callback_start_ns
            else
                callback_start_ns = time_ns()
                lx = Vector(objective_t.lx(x, u))
                lu_ = Vector(objective_t.lu(x, u))
                data.objective += objective_t.l(x, u)[1]
                first_order_ns += time_ns() - callback_start_ns
                callback_start_ns = time_ns()
                lxx = _derivative_matrix(objective_t.lxx(x, u), nx, nx)
                lux = _derivative_matrix(objective_t.lux(x, u), nu, nx)
                luu = _derivative_matrix(objective_t.luu(x, u), nu, nu)
                second_order_ns += time_ns() - callback_start_ns
            end
            
            if nc > 0
                callback_start_ns = time_ns()
                c = Vector{T}(constraints.c(x, u))
                cx = _derivative_matrix(constraints.cx(x, u), nc, nx)
                cu = _derivative_matrix(constraints.cu(x, u), nc, nu)
                first_order_ns += time_ns() - callback_start_ns
                callback_start_ns = time_ns()
                cxx = _derivative_matrix(constraints.cxx(x, u, ϕ), nx, nx)
                cux = _derivative_matrix(constraints.cux(x, u, ϕ), nu, nx)
                cuu = _derivative_matrix(constraints.cuu(x, u, ϕ), nu, nu)
                second_order_ns += time_ns() - callback_start_ns
                
                # evaluate constraint violation norms
                data.primal_1_curr += norm(c, 1)
                data.primal_inf = max(data.primal_inf, norm(c, Inf))
                if feasibility_diagnostic
                    equality_ssq += sum(abs2, c)
                    equality_count += length(c)
                    if !isempty(c)
                        value, index = findmax(abs, c)
                        if value > equality_max
                            equality_max = value
                            equality_worst_stage = t
                            equality_worst_index = index
                        end
                    end
                end
            end

            if feasibility_diagnostic
                if t < ocp.N
                    dynamics_residual = traj[t+1].x - dynamics.f(x, u)
                    dynamics_ssq += sum(abs2, dynamics_residual)
                    dynamics_count += length(dynamics_residual)
                    if !isempty(dynamics_residual)
                        value, index = findmax(abs, dynamics_residual)
                        if value > dynamics_max
                            dynamics_max = value
                            dynamics_worst_stage = t
                            dynamics_worst_index = index
                        end
                    end
                end
                @inbounds for i in eachindex(u)
                    if cl.maskl[i]
                        violation = max(cl.l[i] - u[i], zero(T))
                        bound_ssq += violation^2; bound_count += 1
                        if violation > bound_max
                            bound_max = violation; bound_worst_stage = t
                            bound_worst_index = i; bound_worst_kind = "lower"
                        end
                    end
                    if cl.masku[i]
                        violation = max(u[i] - cl.u[i], zero(T))
                        bound_ssq += violation^2; bound_count += 1
                        if violation > bound_max
                            bound_max = violation; bound_worst_stage = t
                            bound_worst_index = i; bound_worst_kind = "upper"
                        end
                    end
                end
            end

            callback_start_ns = time_ns()
            fx = _derivative_matrix(dynamics.fx(x, u), nx, nx)
            fu = _derivative_matrix(dynamics.fu(x, u), nx, nu)
            first_order_ns += time_ns() - callback_start_ns
            callback_start_ns = time_ns()
            fxx = _derivative_matrix(nc == 0 ? dynamics.fxx(x, u, V̂x) : dynamics.fxx(x, u, λ), nx, nx)
            fux = _derivative_matrix(nc == 0 ? dynamics.fux(x, u, V̂x) : dynamics.fux(x, u, λ), nu, nx)
            fuu = _derivative_matrix(nc == 0 ? dynamics.fuu(x, u, V̂x) : dynamics.fuu(x, u, λ), nu, nu)
            second_order_ns += time_ns() - callback_start_ns
            derivative_s = (time_ns() - derivative_start_ns) / 1e9
            first_order_s = first_order_ns / 1e9
            second_order_s = second_order_ns / 1e9
            derivative_alloc_bytes = memory_diagnostic ? Base.gc_bytes() - stage_alloc_start : 0
            algebra_alloc_start = memory_diagnostic ? Base.gc_bytes() : 0
            algebra_start_ns = time_ns()

            # evaluate barrier Lagrangian

            ul = u - cl.l
            uu = cl.u - u
            data.barrier_lagrangian_curr -= μ* (dot(log.(ul), cl.maskl) + dot(log.(uu), cl.masku))

            # evaluate complementary slackness errors
            cs_l = ul .* zl
            cs_u = uu .* zu
            data.cs_inf_0 = max(data.cs_inf_0, norm(cs_l .* cl.maskl, Inf))
            data.cs_inf_μ = max(data.cs_inf_μ, norm((cs_l .- μ).* cl.maskl, Inf))
            data.cs_inf_0 = max(data.cs_inf_0, norm(cs_u .* cl.masku, Inf))
            data.cs_inf_μ = max(data.cs_inf_μ, norm((cs_u .- μ) .* cl.masku, Inf))

            inv_ul = inv.(ul) .* cl.maskl
            inv_uu = inv.(uu) .* cl.masku

            # Qû = Lu' -μŪ^{-1}e + fu' * V̂x
            ps_m1 = time_ns()
            _ps0 = time_ns()
            Qû = lu_ + fu' * V̂x + μ .* (inv_uu - inv_ul)
            # C = Lxx + fx' * Vxx * fx + V̄x ⋅ fxx
            lean_fx = _lean_value() && V̂xx isa Matrix && _is_identity(fx)
            C = lean_fx ? _lean_C(lxx, V̂xx, fxx) : lxx + fx' * V̂xx * fx + fxx
            _ps1 = time_ns() - _ps0; ps_seq_ns += _ps1; _PARSIM_SUB[1] = _ps1 / 1e9
            ps_m2 = time_ns()
    
            sparse_stage = issparse(luu) || issparse(fu) || (nc > 0 && (issparse(cu) || issparse(cuu)))
            structured_B = sparse_stage && issparse(fu) && nnz(lux) == 0 &&
                nnz(fux) == 0 && (nc == 0 || nnz(cux) == 0)
            # FILTERDDP_STRUCTURED_DYNAMICS also treats a stage whose l_ux is
            # confined to the rows fu acts on (the soft terminal SOC penalty
            # couples x only to the battery powers) as structured: B is then
            # still nonzero only in those rows.
            lux_in_B = false
            if !structured_B && _structured_dynamics() && sparse_stage && issparse(fu) &&
                    issparse(lux) && nnz(fux) == 0 && (nc == 0 || nnz(cux) == 0)
                lux_in_B = issubset(unique(findnz(lux)[1]), unique(findnz(fu)[2]))
                structured_B = lux_in_B
            end
            active_B_rows = Int[]
            B_active = Matrix{T}(undef, 0, 0)
            ux_tmp = Matrix{T}(undef, 0, 0)
            # Ĥ = Luu + Σ + fu' * Vxx * fu + V̄x ⋅ fuu
            Σ_L = inv_ul .* zl
            Σ_U = inv_uu .* zu
            # Exact Hessian with the tree solver: the battery curvature goes to
            # the solver, not into K (see _fixed_pattern_hessian).
            exact_batt = _tree_kkt_enabled() && !_diag_hessian_enabled() && sparse_stage &&
                structured_B && nc > 0 && get(ENV, "FILTERDDP_FACTOR_BACKED_POLICY", "0") == "1" &&
                !(haskey(ENV, "FILTERDDP_CAPTURE_KKT") &&
                  t == parse(Int, get(ENV, "FILTERDDP_CAPTURE_STAGE", "1")))
            battery_curvature = nothing
            lean_hd = nothing
            ps_m3 = time_ns()
            if sparse_stage
                fu_sparse = sparse(fu)
                if exact_batt
                    Ĥ = _fixed_pattern_hessian(nu, Σ_L + Σ_U, luu, fuu, cuu)
                elseif _diag_hessian_enabled() && _direct_diag_hessian_enabled()
                    _ps0 = time_ns()
                    curvature_diag = _diagonal_quadratic_form(fu_sparse, V̂xx, nu)
                    _ps1 = time_ns() - _ps0; ps_seq_ns += _ps1; _PARSIM_SUB[2] = _ps1 / 1e9
                    _lagged_curvature_enabled() && (_LAGGED_CURV[t] = curvature_diag)
                    if _lean_value() && luu isa SparseMatrixCSC && fuu isa SparseMatrixCSC
                        lean_hd = _dense_diag(luu) .+ (Σ_L + Σ_U + curvature_diag)
                        lean_hd .+= _dense_diag(fuu)
                        Ĥ = luu                      # not used: the diagonal is in lean_hd
                    else
                        Ĥ = sparse(luu) + spdiagm(0 => Σ_L + Σ_U + curvature_diag) + sparse(fuu)
                    end
                else
                    Ĥ = sparse(luu) + spdiagm(0 => Σ_L + Σ_U) +
                         fu_sparse' * sparse(V̂xx) * fu_sparse + sparse(fuu)
                end
            else
                ux_tmp = fu' * V̂xx
                Ĥ = luu + diagm(Σ_L) + diagm(Σ_U) + ux_tmp * fu + fuu
            end
            ps_m4 = time_ns()
            # B = Lux + fu' * Vxx * fx + V̄x ⋅ fux
            if structured_B
                active_B_rows = sort!(unique(findnz(fu)[2]))
                # FILTERDDP_STRUCTURED_DYNAMICS keeps fu's active columns sparse
                # (one entry per battery for MPOPF), so this is not a dense
                # n_B x n_x x n_x product; same values up to the sign of zeros.
                fu_active = _structured_dynamics() ? fu_sparse[:, active_B_rows] :
                    Matrix(@view fu[:, active_B_rows])
                _ps0 = time_ns()
                B_active = (lean_fx && fu_active isa SparseMatrixCSC) ? _lean_B(fu_active, V̂xx) :
                    Matrix(fu_active' * V̂xx * fx)
                lux_in_B && (B_active .+= lux[active_B_rows, :])
                exact_batt && (battery_curvature = Matrix(fu_active' * V̂xx * fu_active))
                _ps1 = time_ns() - _ps0; ps_seq_ns += _ps1; _PARSIM_SUB[3] = _ps1 / 1e9
            else
                isempty(ux_tmp) && (ux_tmp = fu' * V̂xx)
                B = ux_tmp * fx
                B .+= lux
                B .+= fux
            end

            ps_m5 = time_ns()
            if nc > 0
                data.barrier_lagrangian_curr += dot(c, ϕ)
                Qû = Qû + cu' * ϕ
                (lean_fx && cxx isa SparseMatrixCSC && nnz(cxx) == 0) || (C = C + cxx)
                if !isnothing(lean_hd)
                    lean_hd .+= _dense_diag(cuu)
                elseif !exact_batt
                    Ĥ = Ĥ + cuu
                end
                !structured_B && (B .+= cux)
            end
            
            ps_m6 = time_ns()
            # inertia correction / regularisation
            if !iszero(reg)
                if !isnothing(lean_hd)
                    lean_hd .+= reg
                else
                    @inbounds for i in axes(Ĥ, 1)
                        Ĥ[i, i] += reg
                    end
                end
            end

            # Opt-in diagonal curvature model; see _diagonalise_hessian above.
            # Applied after the inertia correction so that reg still reaches the
            # diagonal, and before sparse_kkt is decided so the branch is unchanged.
            _diag_hessian_enabled() && (Ĥ = isnothing(lean_hd) ? _diagonalise_hessian(Ĥ) : _floored_diagonal(lean_hd))

            # Sparse network models use the full saddle-point system directly,
            # avoiding a dense QR basis and the explicit reduced Hessian Z'HZ.
            sparse_kkt = nc > 0 && (issparse(Ĥ) || issparse(cu))
            algebra_s = (time_ns() - algebra_start_ns) / 1e9
            algebra_alloc_bytes = memory_diagnostic ? Base.gc_bytes() - algebra_alloc_start : 0
            kkt_assembly_s = 0.0
            factor_s = 0.0
            solve_s = 0.0
            kkt_alloc_bytes = 0
            factor_alloc_bytes = 0
            solve_alloc_bytes = 0
            blocked_value = false
            blocked_α = Vector{T}()
            blocked_ψ = Vector{T}()
            blocked_Vxx = Matrix{T}(undef, 0, 0)
            blocked_Vx = Vector{T}()
            if sparse_kkt
                kkt_alloc_start = memory_diagnostic ? Base.gc_bytes() : 0
                kkt_start_ns = time_ns()
                (exact_batt || !isnothing(lean_hd)) || (Ĥ = sparse(Symmetric(Ĥ)))
                cu_sparse = _stale_jacobian(t, data.k, sparse(cu))
                K = _kkt_pattern_cache_enabled() ?
                    _cached_kkt!((objectid(solver), t), Ĥ, cu_sparse, nu, nc) :
                    [Ĥ sparse(cu_sparse'); cu_sparse spzeros(T, nc, nc)]
                # Fill-in diagnostic: nnz of the coefficient before factorisation
                # and of the LU factors after, per stage. Opt-in; this is what
                # distinguishes "the matrix got denser" from "pivoting got worse"
                # as the explanation for factorisation cost rising ~6x on an
                # identically sized K when the batteries actually cycle.
                nnz_diagnostic = get(ENV, "FILTERDDP_NNZ_DIAGNOSTIC", "0") == "1"

                capture_this_kkt = haskey(ENV, "FILTERDDP_CAPTURE_KKT") &&
                    t == parse(Int, get(ENV, "FILTERDDP_CAPTURE_STAGE", "1"))
                blocked_value = get(ENV, "FILTERDDP_BLOCKED_VALUE_RHS", "0") == "1" &&
                    get(ENV, "FILTERDDP_FACTOR_BACKED_POLICY", "0") == "1" &&
                    structured_B && !capture_this_kkt
                # FILTERDDP_BATTERY_SCHUR reuses the blocked_value bookkeeping:
                # it also returns alpha, psi and the two value increments.
                battery_schur = (_battery_schur_enabled() || _tree_kkt_enabled()) &&
                    get(ENV, "FILTERDDP_FACTOR_BACKED_POLICY", "0") == "1" &&
                    structured_B && !capture_this_kkt
                tree_kkt = battery_schur && _tree_kkt_enabled()
                battery_schur && (blocked_value = true)
                block_width = min(parse(Int, get(ENV, "FILTERDDP_VALUE_BLOCK_WIDTH", "128")), nx)
                rhs_width = battery_schur ? 1 : blocked_value ? block_width : nx + 1
                if isnothing(rhs_workspace)
                    rhs_workspace = Matrix{T}(undef, nu + nc, rhs_width)
                end
                rhs = rhs_workspace
                @views @. rhs[1:nu, 1] = -Qû
                @views @. rhs[nu+1:end, 1] = -c
                if !blocked_value
                    @views begin
                    if structured_B
                        fill!(rhs[1:nu, 2:nx+1], zero(T))
                        @. rhs[active_B_rows, 2:nx+1] = -B_active
                    else
                        @. rhs[1:nu, 2:nx+1] = -B
                    end
                    @. rhs[nu+1:end, 2:nx+1] = -cx
                    end
                end
                captured_rhs = capture_this_kkt ? copy(rhs) : nothing
                # DIAGNOSTIC: snapshot the multi-RHS block *before* ldiv! overwrites
                # it in place, so the periodic capture below records the true input
                # to the KKT solve rather than its solution.
                periodic_capture_active = haskey(ENV, "FILTERDDP_PERIODIC_CAPTURE_DIR") &&
                    data.k % parse(Int, get(ENV, "FILTERDDP_PERIODIC_CAPTURE_STRIDE", "5")) == 0
                periodic_rhs_snapshot = periodic_capture_active ? copy(rhs) : nothing
                kkt_assembly_s = (time_ns() - kkt_start_ns) / 1e9
                kkt_alloc_bytes = memory_diagnostic ? Base.gc_bytes() - kkt_alloc_start : 0
                kkt_solution = Matrix{T}(undef, 0, 0)
                F = nothing
                try
                    if timing_diagnostic || memory_diagnostic
                        factor_alloc_start = memory_diagnostic ? Base.gc_bytes() : 0
                        factor_start_ns = time_ns()
                        F = tree_kkt ? TREE_KKT_HOOK[](K, battery_curvature) : _stage_factor(t, K, data.k, lagged_factors)
                        if nnz_diagnostic
                            # L and U extraction is costly, hence opt-in only
                            nL = nnz(F.L); nU = nnz(F.U)
                            @printf("FILTERDDP_NNZ iteration=%d stage=%d n=%d nnz_K=%d nnz_LU=%d fill_ratio=%.3f
",
                                    data.k, t, size(K, 1), nnz(K), nL + nU, (nL + nU) / max(nnz(K), 1))
                            flush(stdout)
                        end
                        factor_s = (time_ns() - factor_start_ns) / 1e9
                        factor_alloc_bytes = memory_diagnostic ? Base.gc_bytes() - factor_alloc_start : 0
                        solve_alloc_start = memory_diagnostic ? Base.gc_bytes() : 0
                        solve_start_ns = time_ns()
                        if battery_schur
                            blocked_α, blocked_ψ, blocked_Vxx, blocked_Vx =
                                _battery_schur_value(K, F, rhs, nu, active_B_rows, B_active, cx)
                            all(isfinite, blocked_Vxx) || (data.status = 1)
                        elseif blocked_value
                            ldiv!(F, @view(rhs[:, 1:1]))
                            blocked_α = copy(@view rhs[1:nu, 1])
                            blocked_ψ = copy(@view rhs[nu+1:end, 1])
                            blocked_Vxx = zeros(T, nx, nx)
                            blocked_Vx = zeros(T, nx)
                            for first_col in 1:block_width:nx
                                columns = first_col:min(first_col + block_width - 1, nx)
                                width = length(columns)
                                rhs_block = @view rhs[:, 1:width]
                                fill!(rhs_block, zero(T))
                                @views @. rhs[active_B_rows, 1:width] = -B_active[:, columns]
                                @views @. rhs[nu+1:end, 1:width] = -cx[:, columns]
                                ldiv!(F, rhs_block)
                                β_block = @view rhs[1:nu, 1:width]
                                ω_block = @view rhs[nu+1:end, 1:width]
                                @views blocked_Vxx[columns, :] .=
                                    β_block[active_B_rows, :]' * B_active + ω_block' * cx
                                @views blocked_Vx[columns] .= β_block' * Qû + ω_block' * c
                                all(isfinite, rhs_block) || (data.status = 1; break)
                            end
                        else
                            _kkt_wide_solve!(F, rhs)
                            kkt_solution = rhs
                        end
                        solve_s = (time_ns() - solve_start_ns) / 1e9
                        solve_alloc_bytes = memory_diagnostic ? Base.gc_bytes() - solve_alloc_start : 0
                    else
                        F = tree_kkt ? TREE_KKT_HOOK[](K, battery_curvature) : _stage_factor(t, K, data.k, lagged_factors)
                        if nnz_diagnostic
                            # L and U extraction is costly, hence opt-in only
                            nL = nnz(F.L); nU = nnz(F.U)
                            @printf("FILTERDDP_NNZ iteration=%d stage=%d n=%d nnz_K=%d nnz_LU=%d fill_ratio=%.3f
",
                                    data.k, t, size(K, 1), nnz(K), nL + nU, (nL + nU) / max(nnz(K), 1))
                            flush(stdout)
                        end
                        if battery_schur
                            blocked_α, blocked_ψ, blocked_Vxx, blocked_Vx =
                                _battery_schur_value(K, F, rhs, nu, active_B_rows, B_active, cx)
                            all(isfinite, blocked_Vxx) || (data.status = 1)
                        elseif blocked_value
                            ldiv!(F, @view(rhs[:, 1:1]))
                            blocked_α = copy(@view rhs[1:nu, 1])
                            blocked_ψ = copy(@view rhs[nu+1:end, 1])
                            blocked_Vxx = zeros(T, nx, nx)
                            blocked_Vx = zeros(T, nx)
                            for first_col in 1:block_width:nx
                                columns = first_col:min(first_col + block_width - 1, nx)
                                width = length(columns)
                                rhs_block = @view rhs[:, 1:width]
                                fill!(rhs_block, zero(T))
                                @views @. rhs[active_B_rows, 1:width] = -B_active[:, columns]
                                @views @. rhs[nu+1:end, 1:width] = -cx[:, columns]
                                ldiv!(F, rhs_block)
                                β_block = @view rhs[1:nu, 1:width]
                                ω_block = @view rhs[nu+1:end, 1:width]
                                @views blocked_Vxx[columns, :] .=
                                    β_block[active_B_rows, :]' * B_active + ω_block' * cx
                                @views blocked_Vx[columns] .= β_block' * Qû + ω_block' * c
                                all(isfinite, rhs_block) || (data.status = 1; break)
                            end
                        else
                            _kkt_wide_solve!(F, rhs)
                            kkt_solution = rhs
                        end
                    end
                    !blocked_value && !all(isfinite, kkt_solution) && (data.status = 1)
                catch err
                    verbose && @warn "Sparse KKT factorization failed" exception=(err, catch_backtrace())
                    data.status = 1
                end
            elseif nc > 0
                A = cu'
                qrf = qr(A)
                Q = Matrix(qrf.Q * Matrix{T}(I, nu, nu))
                Y = Q[:, 1:nc]
                Z = Q[:, nc+1:nu]
                
                AY = A' * Y
                fk = lu(AY)
                α_β_y = fk \ [-c -cx]

                Ĥ = Symmetric(Ĥ)
                M = Symmetric(Z' * Ĥ * Z)
                ck = cholesky(M; check=false)
            else
                ck = cholesky(Symmetric(Ĥ); check=false)
            end
            !sparse_kkt && ck.info != 0 && (data.status = 1)
            
            if data.status == 1
                if iszero(reg) # initial setting of regularisation
                    reg = (data.reg_last == 0.0) ? options.reg_1 : max(options.reg_min, options.κ_w_m * data.reg_last)
                elseif get(ENV, "FILTERDDP_PERIODIC_CAPTURE_KKT_ONLY", "0") == "1"
                    capture_file = joinpath(capture_dir,
                        @sprintf("iter%04d_stage%03d.jls", data.k, t))
                    Serialization.serialize(capture_file, (
                        iteration=data.k, stage=t, nx=nx, nu=nu, nc=nc,
                        K=copy(K), rhs=copy(@view periodic_rhs_snapshot[:, 1:1]),
                        barrier_mu=μ, reg=reg))
                else
                    reg = (data.reg_last == 0.0) ? options.κ_̄w_p * reg : options.κ_w_p * reg
                end
                break
            end

            update_start_ns = time_ns()
            update_alloc_start = memory_diagnostic ? Base.gc_bytes() : 0
            if sparse_kkt
                update_rule = solver.update[t]
                copyto!(update_rule.α, blocked_value ? blocked_α : @view(kkt_solution[1:nu, 1]))
                copyto!(update_rule.ψ, blocked_value ? blocked_ψ : @view(kkt_solution[nu+1:end, 1]))
                α = update_rule.α
                ψ = update_rule.ψ
                factor_backed_policy = get(ENV, "FILTERDDP_FACTOR_BACKED_POLICY", "0") == "1" && structured_B
                if factor_backed_policy
                    β = blocked_value ? zeros(T, 0, 0) : @view(kkt_solution[1:nu, 2:nx+1])
                    ω = blocked_value ? zeros(T, 0, 0) : @view(kkt_solution[nu+1:end, 2:nx+1])
                    update_rule.β = zeros(T, 0, 0)
                    update_rule.ω = zeros(T, 0, 0)
                    update_rule.factor_policy = FactorBackedPolicy{T}(
                        F, copy(active_B_rows), copy(B_active), sparse(cx),
                        zeros(T, nu + nc), zeros(T, nx), zeros(T, 0, 0))
                else
                    size(update_rule.β) == (nu, nx) || (update_rule.β = zeros(T, nu, nx))
                    size(update_rule.ω) == (nc, nx) || (update_rule.ω = zeros(T, nc, nx))
                    copyto!(update_rule.β, @view kkt_solution[1:nu, 2:nx+1])
                    copyto!(update_rule.ω, @view kkt_solution[nu+1:end, 2:nx+1])
                    update_rule.factor_policy = nothing
                    β = update_rule.β
                    ω = update_rule.ω
                end
            elseif nc > 0
                α_β_z = ck \ (Z' * ([-Qû -B] - Ĥ * Y * α_β_y))
                α_β = Y * α_β_y + Z * α_β_z
                α = α_β[:, 1]
                β = α_β[:, 2:nx+1]

                ψ_ω = ((Y' * ([-Qû -B] - Ĥ * α_β))' / fk)'
                ψ = ψ_ω[:, 1]
                ω = ψ_ω[:, 2:nx+1]
            else
                α_β = ck \ [-Qû -B]
                α = α_β[:, 1]
                β = α_β[:, 2:nx+1]
                ψ = zeros(T, nc)
                ω = zeros(T, nc, nx)
            end

            if sparse_kkt && capture_this_kkt
                Serialization.serialize(ENV["FILTERDDP_CAPTURE_KKT"], (
                    K=K, rhs=captured_rhs, stage=t, nx=nx, nu=nu, nc=nc,
                    x=copy(x), u=copy(u), phi=copy(ϕ), zl=copy(zl), zu=copy(zu),
                    barrier_mu=μ, future_Vx=copy(V̂x), future_Vxx=copy(V̂xx),
                    next_state=copy(dynamics.f(x, u)), alpha=copy(α),
                    beta=copy(β), psi=copy(ψ), omega=copy(ω)))
            end

            # DIAGNOSTIC (not part of the paper's optimization stack): dump every
            # matrix/vector FilterDDP carries per stage, every FILTERDDP_PERIODIC_
            # CAPTURE_STRIDE-th outer iteration, for offline shape/compressibility
            # inspection. No effect unless FILTERDDP_PERIODIC_CAPTURE_DIR is set.
            #
            # Two modes, since raw dumps scale as O(nu*nx + nc*nx) per snapshot --
            # fine at ieee123/ieee2522 scale (KB-MB), multiple GB per snapshot at
            # large10k scale (nu~54665, nc~42303, nx~1020), which filled a disk
            # during the first large10k attempt. Set FILTERDDP_PERIODIC_CAPTURE_
            # INLINE_STATS=1 to compute the same shape/SVD-rank statistics that
            # analyze_periodic_capture.jl would have computed offline, right here
            # where the arrays already live in memory, and append one small CSV
            # row per snapshot instead of serializing the arrays themselves.
            if sparse_kkt && periodic_capture_active
                capture_dir = ENV["FILTERDDP_PERIODIC_CAPTURE_DIR"]
                mkpath(capture_dir)
                if get(ENV, "FILTERDDP_PERIODIC_CAPTURE_INLINE_STATS", "0") == "1"
                    function _optimal_rank(sv, tol)
                        isempty(sv) && return 0
                        total = sum(abs2, sv)
                        total == 0 && return 0
                        for r in 0:length(sv)
                            tail = r == length(sv) ? 0.0 : sum(abs2, @view sv[r+1:end])
                            sqrt(tail / total) <= tol && return r
                        end
                        return length(sv)
                    end
                    rhs_state = Matrix(@view periodic_rhs_snapshot[:, 2:end])
                    beta_d = Matrix(β)
                    omega_d = Matrix(ω)
                    sv_rhs = isempty(rhs_state) ? Float64[] : svdvals(rhs_state)
                    sv_beta = isempty(beta_d) ? Float64[] : svdvals(beta_d)
                    sv_omega = isempty(omega_d) ? Float64[] : svdvals(omega_d)
                    csv_path = joinpath(capture_dir, "periodic_capture_summary.csv")
                    header_needed = !isfile(csv_path)
                    open(csv_path, "a") do io
                        header_needed && println(io,
                            "iteration,stage,nx,nu,nc,K_shape,K_nnz,K_density,rhs_state_shape,beta_shape,omega_shape,alpha_len,psi_len,SigmaL_len,SigmaU_len,Vxx_shape,Vx_len,rhs_sv_max,rhs_sv_min,rhs_rank_1pct,rhs_rank_5pct,rhs_rank_10pct,beta_rank_1pct,beta_rank_5pct,omega_rank_1pct,omega_rank_5pct,rhs_cond,beta_cond,barrier_mu,reg")
                        @printf(io, "%d,%d,%d,%d,%d,%dx%d,%d,%.6e,%dx%d,%dx%d,%dx%d,%d,%d,%d,%d,%dx%d,%d,%.6e,%.6e,%d,%d,%d,%d,%d,%d,%d,%.3e,%.3e,%.3e,%.3e\n",
                            data.k, t, nx, nu, nc,
                            size(K,1), size(K,2), nnz(K), nnz(K) / (size(K,1)^2),
                            size(rhs_state,1), size(rhs_state,2),
                            size(beta_d,1), size(beta_d,2), size(omega_d,1), size(omega_d,2),
                            length(α), length(ψ), length(Σ_L), length(Σ_U),
                            size(V̂xx,1), size(V̂xx,2), length(V̂x),
                            isempty(sv_rhs) ? NaN : sv_rhs[1], isempty(sv_rhs) ? NaN : sv_rhs[end],
                            _optimal_rank(sv_rhs, 0.01), _optimal_rank(sv_rhs, 0.05), _optimal_rank(sv_rhs, 0.10),
                            _optimal_rank(sv_beta, 0.01), _optimal_rank(sv_beta, 0.05),
                            _optimal_rank(sv_omega, 0.01), _optimal_rank(sv_omega, 0.05),
                            (isempty(sv_rhs) || sv_rhs[end] == 0) ? Inf : sv_rhs[1]/sv_rhs[end],
                            (isempty(sv_beta) || sv_beta[end] == 0) ? Inf : sv_beta[1]/sv_beta[end],
                            μ, reg)
                    end
                else
                    capture_file = joinpath(capture_dir,
                        @sprintf("iter%04d_stage%03d.jls", data.k, t))
                    Serialization.serialize(capture_file, (
                        iteration=data.k, stage=t, nx=nx, nu=nu, nc=nc,
                        x=copy(x), u=copy(u),
                        K=copy(K),
                        rhs_multi=periodic_rhs_snapshot,
                        kkt_solution=copy(kkt_solution),
                        alpha=copy(α), beta=copy(Matrix(β)), psi=copy(ψ),
                        omega=copy(Matrix(ω)),
                        Sigma_L=copy(Σ_L), Sigma_U=copy(Σ_U),
                        Vxx_incoming=copy(V̂xx), Vx_incoming=copy(V̂x),
                        barrier_mu=μ, reg=reg))
                end
            end

            # update parameters of update rule for ineq. dual variables

            if sparse_kkt
                @. update_rule.χl = inv_ul * μ - zl - Σ_L * α
                @. update_rule.χu = inv_uu * μ - zu + Σ_U * α
                copyto!(update_rule.Σ_L, Σ_L)
                copyto!(update_rule.Σ_U, Σ_U)
                χl = update_rule.χl
                χu = update_rule.χu
            else
                χl = inv_ul .* μ - zl - Σ_L .* α
                χu = inv_uu .* μ - zu + Σ_U .* α
                solver.update[t] = UpdateRule{T, nx, nu, nc, nux, ncx}(α, ψ, β, ω, χl, χu, Σ_L, Σ_U)
            end

            # evaluate optimality dual error
            Lu = lu_ - zl + zu + fu' * λ
            if nc > 0
                Lu = Lu + cu' * ϕ
            end

            data.dual_inf = max(data.dual_inf, norm(Lu, Inf))
            z_norm += sum(zl)
            z_norm += sum(zu)
            ϕ_norm += norm(ϕ, 1)

            # Update return V derivatives for next timestep Vxx = C + β' * B + ω' cx
            _ps0 = time_ns()
            if blocked_value
                V̂xx = lean_fx ? (C .+= blocked_Vxx) : C + blocked_Vxx
                V̂x = lx + blocked_Vx + fx' * V̂x
            elseif structured_B
                beta_active = Matrix(@view β[active_B_rows, :])
                beta_B = beta_active' * B_active
                V̂xx = C + beta_B
                V̂x = lx + β' * Qû + fx' * V̂x
            else
                beta_B = β' * B
                V̂xx = C + beta_B
                V̂x = lx + β' * Qû + fx' * V̂x
            end
            λ = lx + fx' * λ
            if nc > 0
                if !blocked_value
                    V̂xx = V̂xx + ω' * cx
                    V̂x = V̂x + ω' * c
                end
                V̂x = V̂x + cx' * ϕ
                λ = λ + cx' * ϕ
            end
            ps_val_ns += time_ns() - _ps0; _PARSIM_SUB[7] = ps_val_ns / 1e9

            # evaluate sufficient decrease condition in forward pass
            data.expected_change_L += dot(Qû, α)
            nc > 0 && (data.expected_change_L += dot(c, ψ))
            update_s = (time_ns() - update_start_ns) / 1e9
            if parsim
                timing_diagnostic || error("FILTERDDP_PARSIM needs FILTERDDP_TIMING_DIAGNOSTIC=1")
                ps_skip = _PARSIM_SKIP_NS[] / 1e9
                ps_stage = (time_ns() - derivative_start_ns) / 1e9 - ps_skip
                ps_pt = sparse_kkt ? stage_phase_times(F) : nothing
                if isnothing(ps_pt)                 # not the tree solver: all of it in sequence
                    ps_seq = ps_stage; ps_post = 0.0
                else
                    ps_seq = (ps_seq_ns + ps_val_ns) / 1e9 + ps_pt[2] +
                        (solve_s - ps_pt[3] - ps_pt[5] - ps_skip)
                    ps_post = ps_pt[5] + (update_s - ps_val_ns / 1e9)
                    fp = solver.update[t].factor_policy
                    isnothing(fp) || isnothing(_PARSIM_BETA_B[]) || (fp.βB = _PARSIM_BETA_B[])
                    _PARSIM_BETA_B[] = nothing
                end
                ps_pre = ps_stage - ps_seq - ps_post
                ps_att[1] += ps_pre; ps_att[2] = max(ps_att[2], ps_pre); ps_att[3] += ps_seq
                ps_att[4] += ps_post; ps_att[5] = max(ps_att[5], ps_post)
                ps_stages += 1
                # Two calls: a format with more than 32 arguments is not specialised and
                # cost about 2 s on the first pass of every timed solve.
                if !isnothing(ps_pt)
                    @printf("FILTERDDP_PARSIM_STAGE iteration=%d stage=%d pre_s=%.6f seq_s=%.6f post_s=%.6f qc_s=%.6f curv_s=%.6f b_s=%.6f finish_s=%.6f block_s=%.6f rhs_s=%.6f rows_s=%.6f vinc_s=%.6f value_s=%.6f",
                        data.k, t, ps_pre, ps_seq, ps_post, _PARSIM_SUB[1], _PARSIM_SUB[2], _PARSIM_SUB[3], ps_pt[2], ps_pt[4],
                        _PARSIM_SUB[4], _PARSIM_SUB[5], _PARSIM_SUB[6] - ps_skip, _PARSIM_SUB[7])
                    @printf(" prepare_s=%.6f net1_s=%.6f net2_s=%.6f deriv_s=%.6f algebra_s=%.6f assembly_s=%.6f update_s=%.6f barrier_s=%.6f struct_s=%.6f hess_s=%.6f bsel_s=%.6f nc_s=%.6f reg_s=%.6f
",
                        ps_pt[1], ps_pt[3], ps_pt[5], derivative_s, algebra_s, kkt_assembly_s, update_s,
                        (ps_m1 - algebra_start_ns) / 1e9, (ps_m3 - ps_m2) / 1e9, (ps_m4 - ps_m3) / 1e9 - _PARSIM_SUB[2],
                        (ps_m5 - ps_m4) / 1e9 - _PARSIM_SUB[3], (ps_m6 - ps_m5) / 1e9,
                        algebra_s - (ps_m6 - algebra_start_ns) / 1e9)
                end
            end
            update_alloc_bytes = memory_diagnostic ? Base.gc_bytes() - update_alloc_start : 0
            timing_diagnostic && @printf(
                "FILTERDDP_TIMING iteration=%d barrier_iteration=%d stage=%d derivative_s=%.9f first_order_s=%.9f second_order_s=%.9f algebra_s=%.9f kkt_assembly_s=%.9f factor_s=%.9f solve_s=%.9f update_s=%.9f K_nnz=%d rhs_cols=%d Vxx_nnz=%d\n",
                data.k, data.j, t, derivative_s, first_order_s, second_order_s, algebra_s, kkt_assembly_s, factor_s, solve_s,
                update_s, sparse_kkt ? nnz(K) : 0, sparse_kkt ? size(rhs, 2) : 0,
                count(!iszero, V̂xx))
            memory_diagnostic && @printf(
                "FILTERDDP_MEMORY stage=%d derivative_alloc_MiB=%.3f algebra_alloc_MiB=%.3f kkt_alloc_MiB=%.3f factor_alloc_MiB=%.3f solve_alloc_MiB=%.3f update_alloc_MiB=%.3f K_MiB=%.3f rhs_MiB=%.3f solution_MiB=%.3f beta_MiB=%.3f omega_MiB=%.3f bound_sens_MiB=%.3f Vxx_MiB=%.3f update_rule_MiB=%.3f maxrss_delta_MiB=%.3f\n",
                t, derivative_alloc_bytes / 2.0^20, algebra_alloc_bytes / 2.0^20,
                kkt_alloc_bytes / 2.0^20, factor_alloc_bytes / 2.0^20,
                solve_alloc_bytes / 2.0^20, update_alloc_bytes / 2.0^20,
                sparse_kkt ? Base.summarysize(K) / 2.0^20 : 0.0,
                sparse_kkt ? Base.summarysize(rhs) / 2.0^20 : 0.0,
                sparse_kkt && kkt_solution !== rhs ? Base.summarysize(kkt_solution) / 2.0^20 : 0.0,
                Base.summarysize(β) / 2.0^20, Base.summarysize(ω) / 2.0^20,
                (Base.summarysize(Σ_L) + Base.summarysize(Σ_U)) / 2.0^20,
                Base.summarysize(V̂xx) / 2.0^20,
                Base.summarysize(solver.update[t]) / 2.0^20,
                max(Sys.maxrss() - stage_maxrss_start, 0) / 2.0^20)
        end
        ps_tot .+= ps_att
        scaling_dual = max(options.s_max, (ϕ_norm + z_norm) / max(ni + nc * ocp.N, 1.0))  / options.s_max
        scaling_cs = max(options.s_max, z_norm / max(ni, 1.0))  / options.s_max
        data.dual_inf /= scaling_dual
        data.cs_inf_0 /= scaling_cs
        data.cs_inf_μ /= scaling_cs
        data.barrier_lagrangian_curr += data.objective
        data.status == 0 && break
    end
    data.reg_last = reg
    parsim && @printf(
        "FILTERDDP_PARSIM_BACKWARD iteration=%d stages=%d pre_sum_s=%.9f pre_max_s=%.9f seq_sum_s=%.9f post_sum_s=%.9f post_max_s=%.9f\n",
        data.k, ps_stages, ps_tot[1], ps_tot[2], ps_tot[3], ps_tot[4], ps_tot[5])
    if feasibility_diagnostic && data.status == 0
        equality_rms = sqrt(equality_ssq / max(equality_count, 1))
        dynamics_rms = sqrt(dynamics_ssq / max(dynamics_count, 1))
        bound_rms = sqrt(bound_ssq / max(bound_count, 1))
        @printf("FILTERDDP_FEASIBILITY iteration=%d barrier_iteration=%d equality_count=%d equality_rms=%.12e equality_max=%.12e equality_worst_stage=%d equality_worst_index=%d dynamics_count=%d dynamics_rms=%.12e dynamics_max=%.12e dynamics_worst_stage=%d dynamics_worst_index=%d bound_count=%d bound_rms=%.12e bound_max=%.12e bound_worst_stage=%d bound_worst_index=%d bound_worst_kind=%s\n",
            data.k, data.j, equality_count, equality_rms, equality_max,
            equality_worst_stage, equality_worst_index,
            dynamics_count, dynamics_rms, dynamics_max,
            dynamics_worst_stage, dynamics_worst_index,
            bound_count, bound_rms, bound_max, bound_worst_stage,
            bound_worst_index, bound_worst_kind)
        flush(stdout)
    end
    data.status != 0 && (verbose && (@warn "Backward pass failure, unable to find an iteration matrix with correct inertia."))
end
