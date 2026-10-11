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
    S0 = Matrix{Float64}(undef, nB, nB)
    @inbounds for l in 1:nB, k in 1:nB; S0[bat[k], bat[l]] = b0.Sh0[k, l]; end
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
    return pint_elem(PintLeaf(t, S0, rhat, Rhat, D, dvec, f, Float64[], zeros(0, 0), Float64[], 0.0))
end

# Step 3 of a stage up to its battery solution, given the value function after it.
function _pp_solve(S, V̂x, V̂xx)
    lay, stat, fac = S.lay, S.stat, S.fac
    Qû = S.lu_ + S.fu' * V̂x + S.mubar
    B_active = DDP4OPF._lean_B(S.fu_active, V̂xx)
    S.lux_in_B && (B_active .+= S.lux[S.active, :])
    Qû = Qû + S.cuϕ
    Vc = DDP4OPF._lean_curvature(S.fu_active, V̂xx)
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
function _pp_rounds!(xs::Vector{Any}, timed::Bool)
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
        if timed
            lvl = 0.0
            for (t, mid) in pairs
                dt, xs[t] = _pint_timed(() -> pint_combine(xs[t], xs[mid]))
                lvl = max(lvl, dt); pair_sum += dt
                done += 1; done % every == 0 && GC.gc(false)
            end
            pair_max += lvl
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

function paired_backward_pass!(solver::DDP4OPF.Solver{T, nx, nu, nc}, ocp, traj, d, options, data, idx) where {T, nx, nu, nc}
    T === Float64 || return false
    isnothing(DDP4OPF._reuse_candidate(solver, d)) || return false      # the solver's reuse pass
    mode = _paired_mode()
    (mode != "0" && !DDP4OPF._diag_hessian_enabled() && DDP4OPF._lean_exact_enabled()) || return false
    (DDP4OPF._lean_value() && _half_block() && DDP4OPF._kkt_pattern_cache_enabled()) || return false
    timed = mode == "count"
    N = ocp.N
    μ = d.μ
    cl = ocp.control_limits
    t_start = time_ns()
    serial = timed || _TP_LAST_DIAG[] !== false || isnothing(_TP_LAYOUT[]) ||
        any(t -> !haskey(DDP4OPF._KKT_PATTERN_CACHE, (objectid(solver), t)), 1:N)
    _TP_LAST_DIAG[] = false
    blas_threads = BLAS.get_num_threads()
    blas_default = _BLAS_DEFAULT[]
    # One BLAS thread per worker throughout, and no switching of it from inside tasks.
    BLAS.set_num_threads(1); _BLAS_DEFAULT[] = 1
    pre_t = zeros(N); post_t = zeros(N)
    S = Vector{Any}(undef, N); E = Vector{Any}(undef, N); R = Vector{Any}(undef, N); P = Vector{Any}(undef, N)
    info = nothing; t_1 = t_2 = t_3 = t_start; err_V = NaN; err_v = NaN
    keep_reuse = DDP4OPF._reuse_barrier()
    keep_rows = DDP4OPF._parsim() || !isnothing(DDP4OPF.AFFINE_DIRECTIONS_HOOK[])
    try
        # 1. stage data and the stage as a function of its end states
        step1 = t -> begin
            t0 = time_ns()
            S[t] = _tp_pre(solver, ocp, traj, t, d.k, μ, false, nx, nu, nc, data, idx, serial)
            E[t] = isnothing(S[t]) ? nothing : _pp_elem(S[t], t, nx)
            pre_t[t] = (time_ns() - t0) / 1e9
        end
        if serial
            for t in N:-1:1; step1(t); end
        else
            Threads.@threads for t in 1:N; step1(t); end
        end
        (any(isnothing, S) || any(isnothing, E)) && return false
        t_1 = time_ns()

        # 2. the rounds
        xs = Vector{Any}(undef, N + 1)
        for t in 1:N; xs[t] = E[t]; end
        xs[N+1] = PintVal(zeros(nx, nx), zeros(nx))
        info = _pp_rounds!(xs, timed)
        t_2 = time_ns()

        if get(ENV, "FILTERDDP_PAIRED_CHECK", "0") != "0"        # the ordinary sweep, to compare
            V̂x = zeros(nx); V̂xx = zeros(nx, nx); err_V = 0.0; err_v = 0.0
            for t in N:-1:1
                out = _tp_seq(S[t], V̂x, V̂xx, false, nu, nc)
                isnothing(out) && return false
                _, V̂x, V̂xx = out
                err_V = max(err_V, _pp_rel(xs[t].V, V̂xx)); err_v = max(err_v, _pp_rel(xs[t].v, V̂x))
            end
            t_2 = time_ns()                                       # not part of the pass
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
            t0 = time_ns()
            Vn = xs[t+1]
            R[t] = _pp_solve(S[t], Vn.v, Vn.V)
            P[t] = isnothing(R[t]) ? nothing : _tp_post(solver, S[t], R[t], t, μ, λs[t], nx, nu, nc, keep_reuse, keep_rows, cl)
            post_t[t] = (time_ns() - t0) / 1e9
        end
        if timed
            for t in N:-1:1; step3(t); end
        else
            Threads.@threads for t in 1:N; step3(t); end
        end
        (any(isnothing, R) || any(isnothing, P)) && return false
        t_3 = time_ns()
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
        stages = Vector{DDP4OPF.ReuseStage{Float64}}(undef, N)
        for t in 1:N; stages[t] = P[t].reuse; end
        DDP4OPF._REUSE_CACHE[] = DDP4OPF.ReusePass{Float64}(objectid(solver), d.k, false, stages, d.primal_1_curr,
            d.primal_inf, d.objective, d.dual_inf, d.cs_inf_0, scaling_cs)
    end
    if get(ENV, "FILTERDDP_TIMING_DIAGNOSTIC", "0") == "1"
        total = (time_ns() - t_start) / 1e9
        @printf("FILTERDDP_PAIRED_PASS iteration=%d barrier_iteration=%d stages=%d mode=%s threads=%d rounds=%d combinations=%d step1_s=%.6f rounds_s=%.6f step3_s=%.6f total_s=%.6f scan_err_Vxx=%.3e scan_err_Vx=%.3e chol_fail=%d\n",
                d.k, d.j, N, mode, Threads.nthreads(), info.rounds, info.merges, (t_1 - t_start) / 1e9,
                (t_2 - t_1) / 1e9, (t_3 - t_2) / 1e9, total, err_V, err_v, _PINT_CHOL_FAIL[])
        # One worker per period and per combination: the same line as FILTERDDP_PARSIM's,
        # with the combinations as a second group of parallel work.
        timed && @printf(
            "FILTERDDP_PARSIM_BACKWARD iteration=%d stages=%d pre_sum_s=%.9f pre_max_s=%.9f seq_sum_s=%.9f post_sum_s=%.9f post_max_s=%.9f pair_sum_s=%.9f pair_max_s=%.9f\n",
            d.k, N, sum(pre_t), maximum(pre_t), 0.0, sum(post_t), maximum(post_t), info.pair_sum, info.pair_max)
    end
    return true
end

function install_paired_pass(data, idx)
    _TP_LAYOUT[] = nothing
    _TP_LAST_DIAG[] = nothing
    DDP4OPF.BACKWARD_PASS_HOOK[] = (solver, ocp, traj, d, options) ->
        paired_backward_pass!(solver, ocp, traj, d, options, data, idx)
    if _paired_mode() != "count" && get(ENV, "FILTERDDP_THREADED_FORWARD", "1") != "0"
        DDP4OPF.AFFINE_DIRECTIONS_HOOK[] = threaded_affine_directions
        DDP4OPF.AFFINE_ROLLOUT_HOOK[] = threaded_affine_rollout!
    end
    println("PAIRED_SWEEP installed: mode=", _paired_mode(), ", ", Threads.nthreads(), " Julia threads")
    return nothing
end
