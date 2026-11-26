# The backward pass with the battery sweep done by pairwise combination of the
# periods, inside the solve (FILTERDDP_PAIRED_SWEEP). Exact battery-block
# Hessian only.
#
# pint_scan.jl showed, beside an ordinary solve, that the value recursion of
# the battery energies can be done in about log2(T) rounds of pairwise
# combinations instead of T steps in sequence (Sarkka and Garcia-Fernandez,
# 2023). Here the solve itself uses it. A pass is then
#
#   1. every stage at once: what threaded_pass.jl does in its phase 1, and the
#      stage written as a function G_t(a, b) of the battery energies before
#      and after it (pint_elem). Nothing from stage t+1.
#   2. the rounds: in round r, every G is combined with the first one of the
#      block to its right; the combinations of a round do not depend on each
#      other. After ceil(log2(T+1)) rounds every stage has the value function
#      of the stages after it.
#   3. every stage at once: its battery block with that value function, its
#      battery solution, the second network solve, the multiplier updates and
#      the stored policy.
#
# Nothing is left that has to go through the periods in order except the
# multiplier recursion of the dual residual (vectors of n_B numbers).
#
# FILTERDDP_PAIRED_SWEEP=count runs all of this in sequence on one BLAS thread
# and times every piece, so that the time with one worker per period and one
# per combination can be counted as for FILTERDDP_PARSIM: the slowest stage of
# steps 1 and 3 and, for each round, its slowest combination. It prints a
# FILTERDDP_PARSIM_BACKWARD line with two more fields, pair_sum_s (all
# combinations) and pair_max_s (the slowest of each round, summed), which
# clean_hp_table.py and parsim_from_log.py take off like the other parallel
# parts. Each combination is timed with the collector held off and the
# collector is run between combinations, inside the pass's wall time: its time
# stays in the count as sequential work.
# FILTERDDP_PAIRED_SWEEP=1 runs steps 1 to 3 on Julia threads instead.
# FILTERDDP_PAIRED_CHECK=1 also does the ordinary sweep and prints how far the
# value functions of the two are apart.
#
# The value functions differ from the sweep's in rounding (1e-9 or so), so the
# iterates are not bit for bit those of the sequential run.

isdefined(@__MODULE__, :PintElem) || include(joinpath(@__DIR__, "pint_scan.jl"))
isdefined(@__MODULE__, :threaded_backward_pass!) || include(joinpath(@__DIR__, "threaded_pass.jl"))

_paired_mode() = get(ENV, "FILTERDDP_PAIRED_SWEEP", "0")
# FILTERDDP_PAIRED_ROUNDS=0: steps 1 and 3 as here, and between them the ordinary
# sweep, stage after stage (_tp_seq of threaded_pass.jl), timed like the rest.
# The solve is the ordinary exact one; its count is then on the accounting of
# the pairwise one (stage times without collections), which the solver's own
# FILTERDDP_PARSIM lines are not.
_paired_rounds() = get(ENV, "FILTERDDP_PAIRED_ROUNDS", "1") != "0"

# FILTERDDP_PAIRED_INPLACE=1: the stages' functions live in storage kept from
# pass to pass, and a combination writes its result over its left operand, with
# three n_B x n_B blocks of scratch per worker. Without it every pass allocates
# three such blocks per stage and six per combination (1.1 GB a pass on med2522
# T=96, 18 GB on large10k), which the collector then has to take back.
_paired_inplace() = get(ENV, "FILTERDDP_PAIRED_INPLACE", "0") != "0"

struct PairWork
    K::Matrix{Float64}; Y1::Matrix{Float64}; Y2::Matrix{Float64}; g::Vector{Float64}
end
PairWork(n::Int) = PairWork(Matrix{Float64}(undef, n, n), Matrix{Float64}(undef, n, n),
                            Matrix{Float64}(undef, n, n), Vector{Float64}(undef, n))
const _PP_SLOTS = PintElem[]                 # one per stage
const _PP_WORK = PairWork[]                  # one per worker
const _PP_END = Ref{Any}(nothing)            # the value function after the last stage: zero

function _pp_storage!(N::Int, n::Int, workers::Int)
    if length(_PP_SLOTS) != N || size(_PP_SLOTS[1].Maa, 1) != n
        empty!(_PP_SLOTS); empty!(_PP_WORK)
        for _ in 1:N
            push!(_PP_SLOTS, PintElem(Matrix{Float64}(undef, n, n), Matrix{Float64}(undef, n, n),
                                      Matrix{Float64}(undef, n, n), Vector{Float64}(undef, n), Vector{Float64}(undef, n)))
        end
        _PP_END[] = PintVal(zeros(n, n), zeros(n))
    end
    while length(_PP_WORK) < workers; push!(_PP_WORK, PairWork(n)); end
    return nothing
end

# pint_elem, written into el: the same sums in the same order.
function _pp_elem!(el::PintElem, Sh0::Matrix{Float64}, bat, f::Vector{Float64}, rhat::Vector{Float64},
                   Rhat::SparseMatrixCSC, D::SparseMatrixCSC, dvec::Vector{Float64})
    n = length(f)
    p = 1 ./ f
    Maa, Mab, Mbb = el.Maa, el.Mab, el.Mbb
    @inbounds for l in 1:n
        bl = bat[l]; pl = p[bl]
        for k in 1:n
            bk = bat[k]
            Mbb[bk, bl] = p[bk] * Sh0[k, l] * pl
        end
    end
    rv = rowvals(Rhat); nzv = nonzeros(Rhat)
    copyto!(Mab, Mbb)                                        # Np + Rp, Rp[i, j] = Rhat[j, i] p[j]
    @inbounds for i in 1:size(Rhat, 2), q in nzrange(Rhat, i)
        j = rv[q]
        Mab[i, j] += nzv[q] * p[j]
    end
    copyto!(Maa, Mab)                                        # ... + Rp' + D
    @inbounds for i in 1:size(Rhat, 2), q in nzrange(Rhat, i)
        j = rv[q]
        Maa[j, i] += nzv[q] * p[j]
    end
    rd = rowvals(D); nzd = nonzeros(D)
    @inbounds for j in 1:size(D, 2), q in nzrange(D, j); Maa[rd[q], j] += nzd[q]; end
    Mab .*= -1
    @inbounds for i in 1:n
        pr = p[i] * rhat[i]
        el.ma[i] = pr + dvec[i]
        el.mb[i] = -pr
    end
    return el
end

# pint_combine with the result in x's storage. y is only read: in a round the
# left operands are distinct and none of them is a right operand.
function _pp_combine!(x::PintElem, y::PintElem, w::PairWork)
    K = w.K; g = w.g
    K .= x.Mbb .+ y.Maa
    g .= x.mb .+ y.ma
    C = cholesky!(Symmetric(K); check=false)
    if !issuccess(C)
        r = pint_combine(x, y)
        copyto!(x.Maa, r.Maa); copyto!(x.Mab, r.Mab); copyto!(x.Mbb, r.Mbb); copyto!(x.ma, r.ma); copyto!(x.mb, r.mb)
        return x
    end
    U = C.U
    transpose!(w.Y1, x.Mab); ldiv!(U', w.Y1)
    copyto!(w.Y2, y.Mab); ldiv!(U', w.Y2)
    ldiv!(U', g)
    BLAS.syrk!('U', 'T', -1.0, w.Y1, 1.0, x.Maa); LinearAlgebra.copytri!(x.Maa, 'U')
    copyto!(x.Mbb, y.Mbb); BLAS.syrk!('U', 'T', -1.0, w.Y2, 1.0, x.Mbb); LinearAlgebra.copytri!(x.Mbb, 'U')
    mul!(x.Mab, w.Y1', w.Y2, -1.0, 0.0)
    mul!(x.ma, w.Y1', g, -1.0, 1.0)
    copyto!(x.mb, y.mb); mul!(x.mb, w.Y2', g, -1.0, 1.0)
    return x
end
function _pp_combine!(x::PintElem, y::PintVal, w::PairWork)
    K = w.K; g = w.g
    K .= x.Mbb .+ y.V
    g .= x.mb .+ y.v
    C = cholesky!(Symmetric(K); check=false)
    if !issuccess(C)
        r = pint_combine(x, y)
        copyto!(x.Maa, r.V); copyto!(x.ma, r.v)
        return PintVal(x.Maa, x.ma)
    end
    U = C.U
    transpose!(w.Y1, x.Mab); ldiv!(U', w.Y1)
    ldiv!(U', g)
    BLAS.syrk!('U', 'T', -1.0, w.Y1, 1.0, x.Maa); LinearAlgebra.copytri!(x.Maa, 'U')
    mul!(x.ma, w.Y1', g, -1.0, 1.0)
    return PintVal(x.Maa, x.ma)
end

# Collector time (s), bytes allocated (GiB) and collections so far. Its time
# follows the volume allocated, not how often it runs: held off for a whole
# pass and run once (med2522 T=24, 0.8 GiB a pass), it took 0.15 s against
# 0.2 s in its own 5 to 8 runs, and the stages got slower by as much
# (logs *_chk_gc_budget*, 2026-10-11). So there is no switch for that here.
_pp_gcstat() = (g = Base.gc_num(); (g.total_time / 1e9, Base.gc_bytes() / 2.0^30, Int(g.pause), Int(g.full_sweep)))

# The stage as a function of its two end states, from its phase-1 record.
function _pp_elem(S, t::Int, nx::Int)
    lay, stat, b0 = S.lay, S.stat, S.b0
    nB = length(S.active); nE = length(lay.E); nN = nE - nB
    fa = S.fu_active
    (nx == nB && all(k -> fa.colptr[k+1] == fa.colptr[k] + 1 && rowvals(fa)[fa.colptr[k]] == k, 1:nB)) || return nothing
    f = [nonzeros(fa)[fa.colptr[k]] for k in 1:nB]
    # reduced right-hand side of the battery rows, without the value gradient
    b = S.b; nz = S.nzc; y = S.y1
    rE = Vector{Float64}(undef, nE)
    @inbounds for a in 1:nE; rE[a] = b[lay.E[a]]; end
    @inbounds for (a, tt, p) in stat.kef; rE[a] -= nz[p] * y[lay.Fc[tt]]; end
    @inbounds for (si, i) in enumerate(lay.iso)
        w = b[i] / nz[stat.iso_d[si]]
        for (a, p) in stat.iso_e[si]; rE[a] -= nz[p] * w; end
    end
    # the halved block without curvature, in battery order, and the energy rows' part
    bat = b0.bat
    inplace = _paired_inplace()
    S0 = inplace ? b0.Sh0 : Matrix{Float64}(undef, nB, nB)
    inplace || @inbounds for l in 1:nB, k in 1:nB; S0[bat[k], bat[l]] = b0.Sh0[k, l]; end
    ninv = zeros(nN)
    I_ = Int[]; J_ = Int[]; V_ = Float64[]
    for (g, rg) in enumerate(stat.granges)
        npg = b0.nP[g]; PNg = b0.PN[g]; rv = rowvals(PNg); nzv = nonzeros(PNg)
        for l in 1:length(rg) - npg
            e = stat.perm[rg[npg + l]] - nB
            ninv[e] = b0.Ninv[g][l]
            for q in nzrange(PNg, l)
                push!(I_, stat.perm[rg[rv[q]]]); push!(J_, e); push!(V_, nzv[q])
            end
        end
    end
    SPN = sparse(I_, J_, V_, nB, nN)
    r0P = rE[1:nB]
    r0N = rE[nB+1:end]
    Nd = Diagonal(ninv)
    rhat = r0P - SPN * (ninv .* r0N)
    L = sparse(S.lux)[S.active, :]
    Rhat = sparse(-L + SPN * (Nd * S.cxE))
    D = sparse(sparse(S.lxx) - S.cxE' * (Nd * S.cxE))
    dvec = S.lx + S.cxϕ + S.cxE' * (ninv .* r0N)
    inplace && return _pp_elem!(_PP_SLOTS[t], S0, bat, f, rhat, Rhat, D, dvec)
    return pint_elem(PintLeaf(t, S0, rhat, Rhat, D, dvec, f, Float64[], zeros(0, 0), Float64[], 0.0))
end

# Step 3 of a stage up to its battery solution, given the value function after it.
function _pp_solve(S, V̂x, V̂xx)
    lay, stat, fac = S.lay, S.stat, S.fac
    Qû = S.lu_ + S.fu' * V̂x + S.mubar
    B_active = SPARDDP._lean_B(S.fu_active, V̂xx)
    S.lux_in_B && (B_active .+= S.lux[S.active, :])
    Qû = Qû + S.cuϕ
    Vc = SPARDDP._lean_curvature(S.fu_active, V̂xx)
    st = _tree_kkt_finishV(stat, fac, S.b0, Vc)
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
    return (; s, xE, Qû, B_active, XB=zeros(0, 0))
end

# The rounds. xs[1:T] are the stages' functions, xs[T+1] the value function
# after the last stage; on return xs[t] is the value function before stage t.
function _pp_rounds!(xs::Vector{Any}, timed::Bool, inplace::Bool=false)
    n = length(xs)
    rounds = 0; merges = 0; pair_sum = 0.0; pair_max = 0.0
    every = _pint_gc_every(size(xs[end].V, 1)); done = 0
    h = 1
    while h < n
        pairs = Tuple{Int, Int}[]
        for i in 1:2h:n
            mid = i + h
            mid > n && continue
            for t in i:mid-1; push!(pairs, (t, mid)); end
        end
        rounds += 1; merges += length(pairs)
        if timed && inplace
            lvl = 0.0; work1 = _PP_WORK[1]
            for (t, mid) in pairs
                dt, xs[t] = _pint_timed(() -> _pp_combine!(xs[t], xs[mid], work1))
                lvl = max(lvl, dt); pair_sum += dt
            end
            pair_max += lvl
        elseif timed
            lvl = 0.0
            for (t, mid) in pairs
                dt, xs[t] = _pint_timed(() -> pint_combine(xs[t], xs[mid]))
                lvl = max(lvl, dt); pair_sum += dt
                done += 1; done % every == 0 && GC.gc(false)
            end
            pair_max += lvl
        elseif inplace
            nw = min(length(_PP_WORK), length(pairs))
            @sync for c in 1:nw
                Threads.@spawn begin
                    local wk = _PP_WORK[$c]
                    for q in $c:nw:length(pairs)
                        local tq, mq = pairs[q]
                        xs[tq] = _pp_combine!(xs[tq], xs[mq], wk)
                    end
                end
            end
        else
            out = Vector{Any}(undef, length(pairs))
            Threads.@threads for q in eachindex(pairs)
                t, mid = pairs[q]
                out[q] = pint_combine(xs[t], xs[mid])
            end
            for (q, (t, _)) in enumerate(pairs); xs[t] = out[q]; end
        end
        h *= 2
    end
    return (; rounds, merges, pair_sum, pair_max)
end

_pp_rel(A, B) = norm(A - B) / max(norm(B), floatmin())

function paired_backward_pass!(solver::SPARDDP.Solver{T, nx, nu, nc}, ocp, traj, d, options, data, idx) where {T, nx, nu, nc}
    T === Float64 || return false
    isnothing(SPARDDP._reuse_candidate(solver, d)) || return false      # the solver's reuse pass
    mode = _paired_mode()
    (mode != "0" && !SPARDDP._diag_hessian_enabled() && SPARDDP._lean_exact_enabled()) || return false
    (SPARDDP._lean_value() && _half_block() && SPARDDP._kkt_pattern_cache_enabled()) || return false
    timed = mode == "count"
    N = ocp.N
    μ = d.μ
    cl = ocp.control_limits
    t_start = time_ns()
    serial = timed || _TP_LAST_DIAG[] !== false || isnothing(_TP_LAYOUT[]) ||
        any(t -> !haskey(SPARDDP._KKT_PATTERN_CACHE, (objectid(solver), t)), 1:N)
    _TP_LAST_DIAG[] = false
    blas_threads = BLAS.get_num_threads()
    blas_default = _BLAS_DEFAULT[]
    # One BLAS thread per worker throughout, and no switching of it from inside tasks.
    BLAS.set_num_threads(1); _BLAS_DEFAULT[] = 1
    pre_t = zeros(N); post_t = zeros(N)
    S = Vector{Any}(undef, N); E = Vector{Any}(undef, N); R = Vector{Any}(undef, N); P = Vector{Any}(undef, N)
    info = nothing; t_1 = t_2 = t_3 = t_start; err_V = NaN; err_v = NaN
    keep_reuse = SPARDDP._reuse_barrier()
    keep_rows = SPARDDP._parsim() || !isnothing(SPARDDP.AFFINE_DIRECTIONS_HOOK[])
    ordinary = !_paired_rounds()
    inplace = _paired_inplace() && !ordinary
    inplace && _pp_storage!(N, nx, timed ? 1 : Threads.nthreads())
    seq_sum = 0.0
    g_0 = g_1 = g_2 = g_3 = _pp_gcstat()
    try
        # 1. stage data and the stage as a function of its end states
        step1 = t -> begin
            t0 = time_ns(); c0 = Base.gc_num().total_time
            S[t] = _tp_pre(solver, ocp, traj, t, d.k, μ, false, nx, nu, nc, data, idx, serial)
            E[t] = isnothing(S[t]) ? nothing : (ordinary ? true : _pp_elem(S[t], t, nx))
            pre_t[t] = (time_ns() - t0) / 1e9
            # counting: the stage's time without the collections that fell in it
            timed && (pre_t[t] -= (Base.gc_num().total_time - c0) / 1e9)
        end
        if serial
            for t in N:-1:1; step1(t); end
        else
            Threads.@threads for t in 1:N; step1(t); end
        end
        (any(isnothing, S) || any(isnothing, E)) && return false
        t_1 = time_ns(); g_1 = _pp_gcstat()

        # 2. the rounds, or the ordinary sweep
        xs = Vector{Any}(undef, N + 1)
        if ordinary
            V̂x = zeros(nx); V̂xx = zeros(nx, nx)
            for t in N:-1:1
                t0 = time_ns(); c0 = Base.gc_num().total_time
                out = _tp_seq(S[t], V̂x, V̂xx, false, nu, nc)
                isnothing(out) && return false
                R[t], V̂x, V̂xx = out
                seq_sum += (time_ns() - t0) / 1e9 - (Base.gc_num().total_time - c0) / 1e9
            end
            info = (; rounds=0, merges=0, pair_sum=0.0, pair_max=0.0)
        else
            for t in 1:N; xs[t] = E[t]; end
            xs[N+1] = inplace ? _PP_END[] : PintVal(zeros(nx, nx), zeros(nx))
            info = _pp_rounds!(xs, timed, inplace)
        end
        t_2 = time_ns(); g_2 = _pp_gcstat()

        if !ordinary && get(ENV, "FILTERDDP_PAIRED_CHECK", "0") != "0"   # the ordinary sweep, to compare
            V̂x = zeros(nx); V̂xx = zeros(nx, nx); err_V = 0.0; err_v = 0.0
            for t in N:-1:1
                out = _tp_seq(S[t], V̂x, V̂xx, false, nu, nc)
                isnothing(out) && return false
                _, V̂x, V̂xx = out
                err_V = max(err_V, _pp_rel(xs[t].V, V̂xx)); err_v = max(err_v, _pp_rel(xs[t].v, V̂x))
            end
            t_2 = time_ns(); g_2 = _pp_gcstat()                   # not part of the pass
        end

        # the multiplier recursion of the dual residual (vectors only)
        λs = Vector{Vector{Float64}}(undef, N)
        λ = zeros(nx)
        for t in N:-1:1
            λs[t] = λ
            λ = S[t].lx + S[t].fx' * λ
            λ = λ + S[t].cxϕ
        end

        # 3. every stage with the value function after it
        step3 = t -> begin
            t0 = time_ns(); c0 = Base.gc_num().total_time
            if !ordinary
                Vn = xs[t+1]
                R[t] = _pp_solve(S[t], Vn.v, Vn.V)
            end
            P[t] = isnothing(R[t]) ? nothing : _tp_post(solver, S[t], R[t], t, μ, λs[t], nx, nu, nc, keep_reuse, keep_rows, cl)
            post_t[t] = (time_ns() - t0) / 1e9
            timed && (post_t[t] -= (Base.gc_num().total_time - c0) / 1e9)
        end
        if timed
            for t in N:-1:1; step3(t); end
        else
            Threads.@threads for t in 1:N; step3(t); end
        end
        (any(isnothing, R) || any(isnothing, P)) && return false
        t_3 = time_ns(); g_3 = _pp_gcstat()
    finally
        _BLAS_DEFAULT[] = blas_default
        BLAS.set_num_threads(blas_threads)
    end

    # the pass's scalars, in the order of the sequential pass
    d.status = 0
    d.barrier_lagrangian_curr = 0.0
    d.primal_1_curr = 0.0
    d.primal_inf = 0.0
    d.cs_inf_μ = 0.0
    d.cs_inf_0 = 0.0
    d.objective = 0.0
    d.dual_inf = 0.0
    d.expected_change_L = 0.0
    ϕ_norm = 0.0
    z_norm = 0.0
    for t in N:-1:1
        s, p = S[t], P[t]
        d.objective += s.lval
        d.primal_1_curr += s.c1
        d.primal_inf = max(d.primal_inf, s.cinf)
        d.barrier_lagrangian_curr -= μ * s.logsum
        d.cs_inf_0 = max(d.cs_inf_0, s.cs0)
        d.cs_inf_μ = max(d.cs_inf_μ, s.csμ)
        d.barrier_lagrangian_curr += s.cϕ
        d.dual_inf = max(d.dual_inf, p.dual)
        z_norm += p.zl_sum
        z_norm += p.zu_sum
        ϕ_norm += p.ϕ1
        d.expected_change_L += p.e1
        d.expected_change_L += p.e2
    end
    ni = (cl.nl + cl.nu) * N
    scaling_dual = max(options.s_max, (ϕ_norm + z_norm) / max(ni + nc * N, 1.0)) / options.s_max
    scaling_cs = max(options.s_max, z_norm / max(ni, 1.0)) / options.s_max
    d.dual_inf /= scaling_dual
    d.cs_inf_0 /= scaling_cs
    d.cs_inf_μ /= scaling_cs
    d.barrier_lagrangian_curr += d.objective
    d.reg_last = 0.0
    if keep_reuse
        stages = Vector{SPARDDP.ReuseStage{Float64}}(undef, N)
        for t in 1:N; stages[t] = P[t].reuse; end
        SPARDDP._REUSE_CACHE[] = SPARDDP.ReusePass{Float64}(objectid(solver), d.k, false, stages, d.primal_1_curr,
            d.primal_inf, d.objective, d.dual_inf, d.cs_inf_0, scaling_cs)
    end
    if get(ENV, "FILTERDDP_TIMING_DIAGNOSTIC", "0") == "1"
        total = (time_ns() - t_start) / 1e9
        @printf("FILTERDDP_PAIRED_PASS iteration=%d barrier_iteration=%d stages=%d mode=%s threads=%d rounds=%d combinations=%d step1_s=%.6f rounds_s=%.6f step3_s=%.6f total_s=%.6f scan_err_Vxx=%.3e scan_err_Vx=%.3e chol_fail=%d\n",
                d.k, d.j, N, mode, Threads.nthreads(), info.rounds, info.merges, (t_1 - t_start) / 1e9,
                (t_2 - t_1) / 1e9, (t_3 - t_2) / 1e9, total, err_V, err_v, _PINT_CHOL_FAIL[])
        # One worker per period and per combination: the same line as FILTERDDP_PARSIM's,
        # with the combinations as a second group of parallel work.
        # Stage times are without the collections that fell in them; those are
        # pre_gc_s and post_gc_s (pair_gc_s between the combinations). A worker of
        # its own collects its own garbage, 1/stages of it; parsim_from_log.py
        # and clean_hp_table.py count it that way, or all of it in sequence.
        timed && @printf(
            "FILTERDDP_PARSIM_BACKWARD iteration=%d stages=%d pre_sum_s=%.9f pre_max_s=%.9f seq_sum_s=%.9f post_sum_s=%.9f post_max_s=%.9f pair_sum_s=%.9f pair_max_s=%.9f pre_gc_s=%.9f post_gc_s=%.9f pair_gc_s=%.9f\n",
            d.k, N, sum(pre_t), maximum(pre_t), seq_sum, sum(post_t), maximum(post_t), info.pair_sum, info.pair_max,
            g_1[1] - g_0[1], g_3[1] - g_2[1], g_2[1] - g_1[1])
        @printf("FILTERDDP_PAIRED_GC iteration=%d inplace=%d step1_gc_s=%.6f rounds_gc_s=%.6f step3_gc_s=%.6f step1_gib=%.3f rounds_gib=%.3f step3_gib=%.3f collections=%d full=%d\n",
                d.k, inplace, g_1[1] - g_0[1], g_2[1] - g_1[1], g_3[1] - g_2[1], g_1[2] - g_0[2], g_2[2] - g_1[2],
                g_3[2] - g_2[2], g_3[3] - g_0[3], g_3[4] - g_0[4])
    end
    return true
end

function install_paired_pass(data, idx)
    _TP_LAYOUT[] = nothing
    _TP_LAST_DIAG[] = nothing
    SPARDDP.BACKWARD_PASS_HOOK[] = (solver, ocp, traj, d, options) ->
        paired_backward_pass!(solver, ocp, traj, d, options, data, idx)
    if _paired_mode() != "count" && get(ENV, "FILTERDDP_THREADED_FORWARD", "1") != "0"
        SPARDDP.AFFINE_DIRECTIONS_HOOK[] = threaded_affine_directions
        SPARDDP.AFFINE_ROLLOUT_HOOK[] = threaded_affine_rollout!
    end
    println("PAIRED_SWEEP installed: mode=", _paired_mode(), ", ", Threads.nthreads(), " Julia threads")
    return nothing
end
