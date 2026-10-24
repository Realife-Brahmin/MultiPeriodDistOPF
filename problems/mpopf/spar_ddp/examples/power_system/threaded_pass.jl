# One worker per period, for real: FilterDDP's backward pass in three phases.
#
# The one-worker-per-period figures (FILTERDDP_PARSIM, parsim_from_log.py) are
# a count taken from a sequential run. This file does the pass that way, on
# Julia threads, so that the count can be checked against a clock on a machine
# with enough cores. By what a stage needs from the next one:
#
#   phase 1, all stages at once: derivatives, barrier terms, the stage matrix,
#            the network factorization and its preparation for the battery
#            block, the first network solve. Nothing from stage t+1.
#   phase 2, in sequence from the last stage: the terms that carry the value
#            function of t+1 (battery right-hand side, battery block and its
#            factorization, battery solution, value update).
#   phase 3, all stages at once: the second network solve, the multiplier
#            updates and the stored policy. This stage's battery solution only.
#
# The operations and their order within a stage are those of
# SPARDDP.backward_pass! for the configuration of run_parsim.sh (tree solver,
# lean value path, halved block, factor-backed policy, cached KKT pattern;
# diagonal Hessian, or the exact one with FILTERDDP_LEAN_EXACT=1), so the
# iterates are the same (check with compare_runs.py). Anything outside that
# configuration, a failed factorization or a non-finite number: the function
# returns false and the solver's own pass runs instead.
#
# install_threaded_pass(data, idx, nu) sets SPARDDP.BACKWARD_PASS_HOOK. One
# FILTERDDP_THREADED_PASS line per pass gives the wall time of each phase.
# No FILTERDDP_FEASIBILITY / FILTERDDP_PARSIM lines are printed by these passes.

using SPARDDP
using LinearAlgebra
using SparseArrays
using Printf

const _TP_LAYOUT = Ref{Any}(nothing)        # tree layout of the current K pattern
const _TP_LAST_DIAG = Ref{Any}(nothing)     # Hessian model of the previous pass

function _tp_layout!(data, idx, nu, K)
    st = _TP_LAYOUT[]
    if isnothing(st) || st.colptr != K.colptr || st.rowval != K.rowval
        lay = tree_kkt_layout(data, idx, nu, K)
        st = (lay=lay, stat=tree_kkt_static(lay, K), colptr=copy(K.colptr), rowval=copy(K.rowval))
        _TP_LAYOUT[] = st
    end
    return st
end

# Phase 1 of stage t. Returns a NamedTuple, or nothing if the stage is outside
# the supported configuration.
function _tp_pre(solver, ocp, traj, t::Int, k::Int, μ, diag::Bool, nx::Int, nu::Int, nc::Int, data, idx, build_layout::Bool)
    DM = SPARDDP._derivative_matrix
    x, u, ϕ, zl, zu = traj[t].x, traj[t].u, traj[t].ϕ, traj[t].zl, traj[t].zu
    cl = ocp.control_limits
    con = SPARDDP.stage_con(ocp, t)
    obj = t == ocp.N ? ocp.term_objective : SPARDDP.stage_obj(ocp, t)
    dyn = ocp.dynamics

    lx = Vector(obj.lx(x, u)); lu_ = Vector(obj.lu(x, u)); lval = obj.l(x, u)[1]
    lxx = DM(obj.lxx(x, u), nx, nx); lux = DM(obj.lux(x, u), nu, nx); luu = DM(obj.luu(x, u), nu, nu)
    c = Vector{Float64}(con.c(x, u))
    cx = DM(con.cx(x, u), nc, nx); cu = DM(con.cu(x, u), nc, nu)
    cxx = DM(con.cxx(x, u, ϕ), nx, nx); cux = DM(con.cux(x, u, ϕ), nu, nx); cuu = DM(con.cuu(x, u, ϕ), nu, nu)
    fx = DM(dyn.fx(x, u), nx, nx); fu = DM(dyn.fu(x, u), nx, nu)
    # Second derivatives of the dynamics are linear in their multiplier; a
    # nonzero one would make them depend on the sweep. Not supported.
    one_x = ones(nx)
    fxx = DM(dyn.fxx(x, u, one_x), nx, nx); fux = DM(dyn.fux(x, u, one_x), nu, nx); fuu = DM(dyn.fuu(x, u, one_x), nu, nu)
    (iszero(fxx) && iszero(fux) && iszero(fuu)) || return nothing
    (SPARDDP._is_identity(fx) && issparse(fu) && issparse(luu) && issparse(fuu) && iszero(cux)) || return nothing
    fu_sparse = sparse(fu)
    active = sort!(unique(findnz(fu)[2]))
    lux_in_B = false
    if !iszero(lux)
        (issparse(lux) && issubset(unique(findnz(lux)[1]), active)) || return nothing
        lux_in_B = true
    end
    fu_active = fu_sparse[:, active]

    ul = u - cl.l
    uu = cl.u - u
    logsum = dot(log.(ul), cl.maskl) + dot(log.(uu), cl.masku)
    cs_l = ul .* zl
    cs_u = uu .* zu
    cs0 = max(norm(cs_l .* cl.maskl, Inf), norm(cs_u .* cl.masku, Inf))
    csμ = max(norm((cs_l .- μ) .* cl.maskl, Inf), norm((cs_u .- μ) .* cl.masku, Inf))
    inv_ul = inv.(ul) .* cl.maskl
    inv_uu = inv.(uu) .* cl.masku
    mubar = μ .* (inv_uu - inv_ul)
    Σ_L = inv_ul .* zl
    Σ_U = inv_uu .* zu
    cuϕ = cu' * ϕ
    cxϕ = cx' * ϕ
    cϕ = dot(c, ϕ)

    luu_d = fuu_d = cuu_d = nothing
    if diag
        # the battery-power entries get their curvature in phase 2
        luu_d = SPARDDP._dense_diag(luu); fuu_d = SPARDDP._dense_diag(fuu); cuu_d = SPARDDP._dense_diag(cuu)
        hd = luu_d .+ (Σ_L + Σ_U + zeros(nu))
        hd .+= fuu_d
        hd .+= cuu_d
        Ĥ = SPARDDP._floored_diagonal(hd)
    else
        Ĥ = SPARDDP._fixed_pattern_hessian(nu, Σ_L + Σ_U, luu, fuu, cuu)
    end
    cu_sparse = SPARDDP._stale_jacobian(t, k, sparse(cu))
    key = (objectid(solver), t)
    K = SPARDDP._cached_kkt!(key, Ĥ, cu_sparse, nu, nc)
    L = build_layout ? _tp_layout!(data, idx, nu, K) : _TP_LAYOUT[]
    (isnothing(L) || L.colptr != K.colptr || L.rowval != K.rowval) && return nothing
    lay, stat = L.lay, L.stat

    fac = network_kkt_factor(lay, K)
    Bd, Up = _blas1(() -> _tree_kkt_prepare(lay, stat, fac, K))
    b0 = nothing
    if !diag
        b0 = _blas1(() -> _tree_kkt_finish0(lay, stat, fac, K, Bd, Up, length(active)))
        isnothing(b0) && return nothing
    end
    nzc = copy(nonzeros(K))
    pre = reverse(lay.post)
    z = zeros(stat.nN); y1 = zeros(lay.n); y2 = zeros(lay.n); t2 = zeros(lay.n)

    # Right-hand side without the value gradient (it only enters the battery
    # powers, which the network solve does not read), and the first network solve.
    Qû0 = lu_ + fu' * zeros(nx) + mubar
    Qû0 = Qû0 + cuϕ
    b = Vector{Float64}(undef, nu + nc)
    @views @. b[1:nu] = -Qû0
    @views @. b[nu+1:end] = -c
    _tp_network_solve!(y1, lay, fac, stat.ooff, pre, z, b)

    cx_s = sparse(cx)
    energy_rows = sort!(unique(cx_s.rowval))
    cxE = cx_s[energy_rows, :]
    vcat(active, nu .+ energy_rows) == lay.E || return nothing
    hmap = diag ? SPARDDP._KKT_PATTERN_CACHE[key].hmap : nothing

    return (; x, u, ϕ, zl, zu, lx, lu_, lval, lxx, lux, c, cx, cu, cxx, fx, fu, fxx, fu_sparse, active, fu_active,
            lux_in_B, logsum, cs0, csμ, inv_ul, inv_uu, mubar, Σ_L, Σ_U, cuϕ, cxϕ, cϕ, luu_d, fuu_d, cuu_d,
            K, lay, stat, fac, Bd, Up, b0, nzc, pre, z, y1, y2, t2, b, cx_s, energy_rows, cxE, hmap,
            c1=norm(c, 1), cinf=norm(c, Inf))
end

# x[network rows] = N \ b[network rows]: tree_kkt_network_solve! of tree_kkt.jl,
# callable before the stage's solver object exists.
function _tp_network_solve!(x, lay, fac, ooff, pre, z, b)
    @inbounds for j in lay.post
        own = lay.own[j]; m = length(own); o = ooff[j]
        F = fac.F[j]; A = F.A; perm = F.perm
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
                z[po+pos[q]] -= v
            end
        end
    end
    @inbounds for j in pre
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

# Phase 2 of stage t: everything that reads the value function of t+1.
# Returns (stage results, V̂x, V̂xx) or nothing on failure.
function _tp_seq(S, V̂x, V̂xx, diag::Bool, nu::Int, nc::Int)
    lay, stat, fac = S.lay, S.stat, S.fac
    Qû = S.lu_ + S.fu' * V̂x + S.mubar
    C = SPARDDP._lean_C(S.lxx, V̂xx, S.fxx)
    B_active = SPARDDP._lean_B(S.fu_active, V̂xx)
    S.lux_in_B && (B_active .+= S.lux[S.active, :])
    Qû = Qû + S.cuϕ
    C = C + S.cxx
    nB = length(S.active)
    if diag
        curv = SPARDDP._diagonal_quadratic_form(S.fu_sparse, V̂xx, nu)
        floor_val = SPARDDP._diag_hessian_floor()
        nzK = nonzeros(S.K)
        @inbounds for j in S.active
            h = S.luu_d[j] + (S.Σ_L[j] + S.Σ_U[j] + curv[j])
            h += S.fuu_d[j]
            h += S.cuu_d[j]
            v = max(h, floor_val)
            nzK[S.hmap[j]] = v
            S.nzc[S.hmap[j]] = v
        end
        st = _blas1(() -> _tree_kkt_finish(lay, stat, fac, S.K, S.Bd, S.Up, false, nothing))
    else
        Vc = SPARDDP._lean_curvature(S.fu_active, V̂xx)
        st = _blas1(() -> _tree_kkt_finishV(stat, fac, S.b0, Vc))
    end
    s = TreeKKTSolver(lay, stat, fac, st, S.nzc, S.pre, S.z, S.y1, S.y2, S.t2, zeros(5))

    # the middle of tree_kkt_solve!: battery right-hand side and battery solution
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

    hv = tree_kkt_value_halved(s, B_active, S.cxE)
    isnothing(hv) && return nothing
    XB, Vinc = hv
    all(isfinite, Vinc) || return nothing
    ψE = zeros(nc)
    @inbounds for (q, r) in enumerate(S.energy_rows); ψE[r] = xE[nB+q, 1]; end
    Vx_inc = B_active' * xE[1:nB, 1] + S.cx_s' * ψE
    V̂xx_new = (C .+= Vinc)
    V̂x_new = S.lx + Vx_inc + S.fx' * V̂x
    V̂x_new = V̂x_new + S.cxϕ
    return (; s, xE, Qû, B_active, XB), V̂x_new, V̂xx_new
end

# Phase 3 of stage t: the rest of tree_kkt_solve!, the multiplier updates and
# the stored policy. λnext is the multiplier recursion's value from stage t+1.
function _tp_post(solver, S, R, t::Int, μ, λnext, nx::Int, nu::Int, nc::Int, keep_reuse::Bool, keep_rows::Bool, cl)
    lay, stat = S.lay, S.stat
    b = S.b; nz = S.nzc; y = S.y1; y2 = S.y2; t2 = S.t2; xE = R.xE
    nE = length(lay.E)
    fill!(t2, 0.0)
    @inbounds for (a, f, p) in stat.kef; t2[lay.Fc[f]] += nz[p] * xE[a, 1]; end
    _tp_network_solve!(y2, lay, S.fac, stat.ooff, S.pre, S.z, t2)
    @inbounds for (si, i) in enumerate(lay.iso)
        v = b[i]
        for (a, p) in stat.iso_e[si]; v -= nz[p] * xE[a, 1]; end
        b[i] = v / nz[stat.iso_d[si]]
    end
    @inbounds for j in eachindex(lay.own), i in lay.own[j]; b[i] = y[i] - y2[i]; end
    @inbounds for a in 1:nE; b[lay.E[a]] = xE[a, 1]; end
    all(isfinite, b) || return nothing

    rule = solver.update[t]
    copyto!(rule.α, @view b[1:nu])
    copyto!(rule.ψ, @view b[nu+1:end])
    α = rule.α
    ψ = rule.ψ
    rule.β = zeros(Float64, 0, 0)
    rule.ω = zeros(Float64, 0, 0)
    # Battery rows of the feedback, for the forward pass's state recursion. The
    # exact Hessian's sweep does not form them (it needs only the value update).
    βB = (keep_rows && isempty(R.XB)) ? tree_kkt_feedback_rows(R.s, R.B_active, S.cxE) : R.XB
    rule.factor_policy = SPARDDP.FactorBackedPolicy{Float64}(
        R.s, copy(S.active), copy(R.B_active), sparse(S.cx),
        zeros(Float64, nu + nc), zeros(Float64, nx), keep_rows ? βB : zeros(Float64, 0, 0))
    inv_ul, inv_uu, Σ_L, Σ_U, zl, zu = S.inv_ul, S.inv_uu, S.Σ_L, S.Σ_U, S.zl, S.zu
    @. rule.χl = inv_ul * μ - zl - Σ_L * α
    @. rule.χu = inv_uu * μ - zu + Σ_U * α
    copyto!(rule.Σ_L, Σ_L)
    copyto!(rule.Σ_U, Σ_U)

    Lu = S.lu_ - zl + zu + S.fu' * λnext
    Lu = Lu + S.cuϕ
    reuse = keep_reuse ? SPARDDP._reuse_stage(S.c, S.lu_, S.cu, S.ϕ, S.lx, S.cx, S.fu, S.u, cl) : nothing
    return (; dual=norm(Lu, Inf), zl_sum=sum(zl), zu_sum=sum(zu), ϕ1=norm(S.ϕ, 1),
            e1=dot(R.Qû, α), e2=dot(S.c, ψ), reuse)
end

function threaded_backward_pass!(solver::SPARDDP.Solver{T, nx, nu, nc}, ocp, traj, d, options, data, idx) where {T, nx, nu, nc}
    T === Float64 || return false
    # the pass after a barrier update keeps its factors: the solver's own code does that
    isnothing(SPARDDP._reuse_candidate(solver, d)) || return false
    diag = SPARDDP._diag_hessian_enabled()
    (diag ? SPARDDP._direct_diag_hessian_enabled() : SPARDDP._lean_exact_enabled()) || return false
    (SPARDDP._lean_value() && _half_block() && SPARDDP._kkt_pattern_cache_enabled()) || return false
    N = ocp.N
    μ = d.μ
    cl = ocp.control_limits
    t_start = time_ns()

    # Phase 1. The KKT pattern cache and the tree layout are shared tables:
    # the pass that fills or rebuilds them runs its stages one after another.
    serial = _TP_LAST_DIAG[] !== diag || isnothing(_TP_LAYOUT[]) ||
        any(t -> !haskey(SPARDDP._KKT_PATTERN_CACHE, (objectid(solver), t)), 1:N)
    _TP_LAST_DIAG[] = diag
    S = Vector{Any}(undef, N)
    blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        if serial
            for t in N:-1:1
                S[t] = _tp_pre(solver, ocp, traj, t, d.k, μ, diag, nx, nu, nc, data, idx, true)
            end
        else
            Threads.@threads for t in 1:N
                S[t] = _tp_pre(solver, ocp, traj, t, d.k, μ, diag, nx, nu, nc, data, idx, false)
            end
        end
    finally
        BLAS.set_num_threads(blas_threads)
    end
    any(isnothing, S) && return false
    t_1 = time_ns()

    # Phase 2, with the multiplier recursion of the dual residual alongside.
    R = Vector{Any}(undef, N)
    λs = Vector{Vector{Float64}}(undef, N)          # λ as stage t receives it
    V̂x = zeros(nx); V̂xx = zeros(nx, nx); λ = zeros(nx)
    for t in N:-1:1
        out = _tp_seq(S[t], V̂x, V̂xx, diag, nu, nc)
        isnothing(out) && return false
        R[t], V̂x, V̂xx = out
        λs[t] = λ
        λ = S[t].lx + S[t].fx' * λ
        λ = λ + S[t].cxϕ
    end
    t_2 = time_ns()

    # Phase 3.
    keep_reuse = SPARDDP._reuse_barrier()
    keep_rows = isdefined(SPARDDP, :AFFINE_DIRECTIONS_HOOK) && !isnothing(SPARDDP.AFFINE_DIRECTIONS_HOOK[])
    P = Vector{Any}(undef, N)
    BLAS.set_num_threads(1)
    try
        Threads.@threads for t in 1:N
            P[t] = _tp_post(solver, S[t], R[t], t, μ, λs[t], nx, nu, nc, keep_reuse, keep_rows, cl)
        end
    finally
        BLAS.set_num_threads(blas_threads)
    end
    any(isnothing, P) && return false
    t_3 = time_ns()

    # The pass's scalars, accumulated in the order of the sequential pass.
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
        SPARDDP._REUSE_CACHE[] = SPARDDP.ReusePass{Float64}(objectid(solver), d.k, diag, stages, d.primal_1_curr,
            d.primal_inf, d.objective, d.dual_inf, d.cs_inf_0, scaling_cs)
    end
    get(ENV, "FILTERDDP_TIMING_DIAGNOSTIC", "0") == "1" && @printf(
        "FILTERDDP_THREADED_PASS iteration=%d barrier_iteration=%d stages=%d threads=%d serial_phase1=%d phase1_s=%.6f phase2_s=%.6f phase3_s=%.6f total_s=%.6f\n",
        d.k, d.j, N, Threads.nthreads(), serial, (t_1 - t_start) / 1e9, (t_2 - t_1) / 1e9, (t_3 - t_2) / 1e9,
        (time_ns() - t_start) / 1e9)
    return true
end

# ------------------------------------------------------------ forward pass --
# The line search for affine dynamics (affine_linesearch.jl), with the periods
# at once. The state depends on the controls only through the battery powers,
# whose feedback rows the backward pass stored (policy.βB), so the states of
# all periods follow from a recursion of n_B x n_x products. Given its state,
# a period's policy solve, directions and ratio test need no other period,
# and neither does its part of a trial point. The states come from the
# recursion, not from the full policy solves, so the directions agree with the
# sequential ones to rounding, not bit for bit.

function threaded_affine_directions(solver, ocp, d, τ::T) where T
    N = ocp.N
    cl = ocp.control_limits
    xs = Vector{Vector{T}}(undef, N)
    xs[1] = solver.nominal[1].x
    for t in 1:N-1
        nom = solver.nominal[t]; rule = solver.update[t]; pol = rule.factor_policy
        (isnothing(pol) || isempty(pol.βB)) && return nothing      # no stored rows: the sequential rollout
        rows = pol.active_rows
        ur = copy(nom.u)
        @views ur[rows] .+= rule.α[rows] .+ pol.βB * (xs[t] .- nom.x)
        xs[t+1] = ocp.dynamics.f(xs[t], ur)
    end
    dx = Vector{Vector{T}}(undef, N); du = Vector{Vector{T}}(undef, N)
    dϕ = Vector{Vector{T}}(undef, N); dzl = Vector{Vector{T}}(undef, N); dzu = Vector{Vector{T}}(undef, N)
    γs = ones(T, N); γds = ones(T, N)
    split = SPARDDP._split_step()
    blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        Threads.@threads for t in 1:N
            nom = solver.nominal[t]; rule = solver.update[t]
            δx = xs[t] - nom.x
            βδx, ωδx = SPARDDP.policy_actions!(rule, δx)
            dx[t] = δx
            du[t] = rule.α .+ βδx
            dϕ[t] = rule.ψ .+ ωδx
            dzl[t] = rule.χl .- rule.Σ_L .* βδx
            dzu[t] = rule.χu .+ rule.Σ_U .* βδx
            dd, el, eu = du[t], dzl[t], dzu[t]
            γ = one(T); γd = one(T)
            @inbounds for i in eachindex(dd)
                if cl.maskl[i]
                    dd[i] < 0 && (γ = min(γ, τ * (nom.u[i] - cl.l[i]) / -dd[i]))
                    if el[i] < 0
                        g = τ * nom.zl[i] / -el[i]
                        split ? (γd = min(γd, g)) : (γ = min(γ, g))
                    end
                end
                if cl.masku[i]
                    dd[i] > 0 && (γ = min(γ, τ * (cl.u[i] - nom.u[i]) / dd[i]))
                    if eu[i] < 0
                        g = τ * nom.zu[i] / -eu[i]
                        split ? (γd = min(γd, g)) : (γ = min(γ, g))
                    end
                end
            end
            γs[t] = γ; γds[t] = γd
        end
    finally
        BLAS.set_num_threads(blas_threads)
    end
    γmax = minimum(γs)
    return SPARDDP.AffineDirections{T}(dx, du, dϕ, dzl, dzu, γmax, split ? minimum(γds) : T(NaN))
end

function threaded_affine_rollout!(solver::SPARDDP.Solver{T, nx, nu, nc}, ocp, d, dir, γ::T) where {T, nx, nu, nc}
    μ = d.μ
    cl = ocp.control_limits
    N = ocp.N
    γz = isnan(dir.γdual) ? γ : dir.γdual
    c1 = zeros(T, N); cϕ = zeros(T, N); logl = zeros(T, N); logu = zeros(T, N); lv = zeros(T, N)
    bad = falses(N)
    Threads.@threads for t in 1:N
        nom = solver.nominal[t]
        x = nom.x .+ γ .* dir.dx[t]
        u = nom.u .+ γ .* dir.du[t]
        ϕ = nom.ϕ .+ γ .* dir.dϕ[t]
        zl = nom.zl .+ γz .* dir.dzl[t]
        zu = nom.zu .+ γz .* dir.dzu[t]
        solver.current[t] = SPARDDP.TrajectoryElement{T, nx, nu, nc}(x, u, ϕ, zl, zu)
        if nc > 0
            c = SPARDDP.stage_con(ocp, t).c(x, u)
            c1[t] = norm(c, 1)
            cϕ[t] = dot(c, ϕ)
        end
        ul = u - cl.l
        uu = cl.u - u
        if any((ul .<= 0) .* cl.maskl) || any((uu .<= 0) .* cl.masku) ||
           any((zl .<= 0) .* cl.maskl) || any((zu .<= 0) .* cl.masku)
            bad[t] = true
        else
            logl[t] = dot(log.(ul), cl.maskl)
            logu[t] = dot(log.(uu), cl.masku)
            lv[t] = t == N ? ocp.term_objective.l(x, u)[1] : SPARDDP.stage_obj(ocp, t).l(x, u)[1]
        end
    end
    d.status = 0
    d.primal_1_next = T(0.0)
    d.barrier_lagrangian_next = T(0.0)
    for t in 1:N                                   # sums in the order of the sequential rollout
        if nc > 0
            d.primal_1_next += c1[t]
            d.barrier_lagrangian_next += cϕ[t]
        end
        if bad[t]
            d.status = 2
            return true
        end
        d.barrier_lagrangian_next -= μ * logl[t]
        d.barrier_lagrangian_next -= μ * logu[t]
        d.barrier_lagrangian_next += lv[t]
    end
    return true
end

function install_threaded_pass(data, idx)
    _TP_LAYOUT[] = nothing
    _TP_LAST_DIAG[] = nothing
    SPARDDP.BACKWARD_PASS_HOOK[] = (solver, ocp, traj, d, options) ->
        threaded_backward_pass!(solver, ocp, traj, d, options, data, idx)
    if isdefined(SPARDDP, :AFFINE_DIRECTIONS_HOOK) && get(ENV, "FILTERDDP_THREADED_FORWARD", "1") != "0"
        SPARDDP.AFFINE_DIRECTIONS_HOOK[] = threaded_affine_directions
        SPARDDP.AFFINE_ROLLOUT_HOOK[] = threaded_affine_rollout!
    end
    println("THREADED_PASS installed: ", Threads.nthreads(), " Julia threads")
    return nothing
end
