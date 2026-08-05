# Stage 6: add inequality constraints to the Stage 5 copper-plate battery problem,
# one layer at a time, using FilterDDP's intended mechanisms.
#
#   6a  grid-power bounds      -> native ControlLimits on Psub  (incl. Psub >= 0)
#   6b  battery-power bounds   -> native ControlLimits on P_B
#   6c  energy bounds          -> STATE bound, needs the slack reformulation
#   6d  all three together
#   6e  the same, at T = 6
#   6f  tolerance sweep on whichever case fails to certify
#   6g  behaviour on a provably infeasible bound set
#
# Objective is the tADMM paper's: linear energy cost at price c^t plus a quadratic
# charge on battery throughput. See copper_plate_battery.jl for the full model.
#
# Independent reference: the OCP reduces exactly to a strictly convex QP in the T
# grid-import powers, solved here by ACTIVE-SET ENUMERATION — every subset of
# inequality rows of size <= n is tried, the equality-constrained KKT system is
# solved, and the unique subset that is primal feasible with nonnegative
# multipliers is the global optimum. Shares no code with FilterDDP.
#
# Reduction (z = Psub, P_B^t = p_L^t - z_t, s = Σ P_B^t = D - Σ z):
#   J(z) = Σ c^t z_t Δt + C_B Δt Σ (p_L^t - z_t)² + w Δt² (D - Σ z)²
#        = ½ zᵀQz + qᵀz + const,
#   Q = 2 C_B Δt I + 2 w Δt² 11ᵀ,  q_t = c^t Δt - 2 C_B Δt p_L^t - 2 w Δt² D
#
# The battery is lossless (no round-trip efficiency term); see the header of
# copper_plate_battery.jl for why a single bidirectional η is not used.

using FilterDDP
using LinearAlgebra
using Printf
using StaticArrays

const T_ = Float64

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

const PL6 = T_[1.0, 1.6, 1.2, 1.8, 0.9, 1.4]      # demand
const C6  = T_[0.30, 0.55, 0.40, 0.62, 0.25, 0.48] # energy price
const C_B = T_(0.5)
const Δt  = T_(1.0)
const B_0 = T_(2.0)
const W   = T_(5.0)

# Lagrange interpolation through nodes τ = 1..n, exact at the nodes and smooth.
# NOTE: degree n-1. Fine for n = 3 and n = 6; see notes/MPOPF_APPLICABILITY.md
# for why this does NOT scale to the horizons an MPOPF actually needs.
function lagrange(v, τ)
    n = length(v)
    acc = zero(τ) * v[1]
    for i = 1:n
        term = v[i]
        for j = 1:n
            j == i && continue
            term = term * (τ - j) / (i - j)
        end
        acc = acc + term
    end
    return acc
end

# ---------------------------------------------------------------------------
# Independent reference: strictly convex QP with linear inequalities G z <= h
# ---------------------------------------------------------------------------

"""
Global optimum of min ½zᵀQz + qᵀz s.t. Gz <= h, for Q ≻ 0, by active-set enumeration.

Depth-first over subsets of active rows of cardinality <= n, evaluating each in
place and returning the first KKT point that is both primal and dual feasible.
Q ≻ 0 makes the minimiser unique, so the first such point is the global optimum.
Returns `nothing` when no KKT point exists, i.e. the feasible set is empty.
"""
function solve_qp(Q, q, G, h; tol = 1e-9)
    n = length(q)
    m = size(G, 1)

    # An inconsistent active set (e.g. both z_i <= hi and z_i >= lo) gives a
    # SINGULAR KKT matrix. Julia's `\` does not always throw on those; it can
    # return a least-squares vector that is coincidentally primal/dual feasible
    # and would then be accepted as "the optimum". So the linear system must be
    # checked as actually solved, not merely as producing a feasible-looking z.
    function try_set(S)
        local z, λS
        try
            if isempty(S)
                z = -(Q \ q); λS = T_[]
                norm(Q * z + q) > 1e-8 * max(1.0, norm(q)) && return nothing
            else
                GS = G[S, :]
                K = [Q GS'; GS zeros(T_, length(S), length(S))]
                rhs = [-q; h[S]]
                sol = K \ rhs
                norm(K * sol - rhs) > 1e-8 * max(1.0, norm(rhs)) && return nothing
                z = sol[1:n]; λS = sol[n+1:end]
            end
        catch
            return nothing
        end
        any(!isfinite, z) && return nothing
        (!isempty(λS) && any(!isfinite, λS)) && return nothing
        all(G * z .<= h .+ tol) || return nothing          # primal feasible
        all(λS .>= -tol) || return nothing                 # dual feasible
        return (z = z, obj = 0.5 * z' * Q * z + q' * z, active = copy(S), λ = λS)
    end

    r = try_set(Int[])
    r === nothing || return r
    kmax = min(n, m)
    S = Int[]
    function dfs(start::Int)
        length(S) >= kmax && return nothing
        for i = start:m
            push!(S, i)
            r = try_set(S)
            r === nothing || (pop!(S); return r)
            r = dfs(i + 1)
            r === nothing || (pop!(S); return r)
            pop!(S)
        end
        return nothing
    end
    return dfs(1)
end

function reference(T::Int; ps_lo = -Inf, ps_hi = Inf, pb_lo = -Inf, pb_hi = Inf,
                   B_lo = -Inf, B_hi = Inf)
    pL = PL6[1:T]; c = C6[1:T]
    D = sum(pL)
    Q = 2 * C_B * Δt * Matrix{T_}(I, T, T) + 2 * W * Δt^2 * ones(T_, T, T)
    q = c .* Δt .- 2 * C_B * Δt .* pL .- 2 * W * Δt^2 * D

    rows = Vector{Vector{T_}}(); rhs = T_[]; tag = String[]
    addrow!(r, s, name) = (push!(rows, r); push!(rhs, s); push!(tag, name))
    for t = 1:T
        et = zeros(T_, T); et[t] = 1.0
        isfinite(ps_hi) && addrow!(copy(et), ps_hi, "Psub[$t] <= $ps_hi")
        isfinite(ps_lo) && addrow!(-copy(et), -ps_lo, "Psub[$t] >= $ps_lo")
        isfinite(pb_hi) && addrow!(-copy(et), pb_hi - pL[t], "P_B[$t] <= $pb_hi")
        isfinite(pb_lo) && addrow!(copy(et), pL[t] - pb_lo, "P_B[$t] >= $pb_lo")
        # B_{t+1} = B_0 - Δt (D_t - S_t),  S_t = Σ_{k<=t} z_k
        st = zeros(T_, T); st[1:t] .= 1.0
        Dt = sum(pL[1:t])
        isfinite(B_hi) && addrow!(copy(st), (B_hi - B_0) / (Δt) + Dt, "B[$(t+1)] <= $B_hi")
        isfinite(B_lo) && addrow!(-copy(st), -((B_lo - B_0) / (Δt) + Dt), "B[$(t+1)] >= $B_lo")
    end
    G = isempty(rows) ? zeros(T_, 0, T) : reduce(vcat, [r' for r in rows])
    sol = solve_qp(Q, q, G, rhs)
    sol === nothing && return nothing
    psub = sol.z
    pb = pL .- psub
    B = zeros(T_, T + 1); B[1] = B_0
    for t = 1:T
        B[t+1] = B[t] - pb[t] * Δt
    end
    J = sum(c .* psub .* Δt) + C_B * sum(pb .^ 2) * Δt + W * (B[T+1] - B_0)^2
    return (psub = psub, pb = pb, B = B, J = J, active = tag[sol.active])
end

# ---------------------------------------------------------------------------
# FilterDDP models
# ---------------------------------------------------------------------------

const STATUS = Dict(0 => "Optimal", 1 => "Backward pass failed", 7 => "Line search failed",
                    8 => "Max iterations")
status_str(s) = get(STATUS, s, "status $s")

"""Cases without energy bounds: nu = 2, nc = 1."""
function solve_ddp_bounds(T::Int; ps_lo = -Inf, ps_hi = Inf, pb_lo = -Inf, pb_hi = Inf,
                          tol = 1e-10)
    nx, nu = 2, 2
    pL = PL6[1:T]; c = C6[1:T]

    dyn = Dynamics((x, u) -> [x[1] - u[2] * Δt, x[2] + 1.0], nx, nu)
    l = (x, u) -> lagrange(c, x[2]) * u[1] * Δt + C_B * u[2]^2 * Δt
    lN = (x, u) -> c[T] * u[1] * Δt + C_B * u[2]^2 * Δt +
                   W * (x[1] - u[2] * Δt - B_0)^2
    cons = EqualityConstraints((x, u) -> [u[1] + u[2] - lagrange(pL, x[2])], nx, nu)
    cl = ControlLimits(SVector{nu,T_}([ps_lo, pb_lo]), SVector{nu,T_}([ps_hi, pb_hi]))

    ocp = build_ocp(T, Objective(l, nx, nu), Objective(lN, nx, nu), dyn, cons, cl)
    solver = Solver(ocp; options = Options{T_}(verbose = false, optimality_tolerance = tol))
    x1 = SVector{nx,T_}([B_0, 1.0])
    ps0 = isfinite(ps_lo) && isfinite(ps_hi) ? (ps_lo + ps_hi) / 2 : 1.0
    pb0 = isfinite(pb_lo) && isfinite(pb_hi) ? (pb_lo + pb_hi) / 2 : 0.0
    ū = [SVector{nu,T_}([ps0, pb0]) for _ = 1:T]
    solve!(solver, x1, ū)
    x, u = get_trajectory(solver)
    psub = [u[t][1] for t = 1:T]; pb = [u[t][2] for t = 1:T]
    B = [x[t][1] for t = 1:T]; push!(B, B[T] - pb[T] * Δt)
    return (solver = solver, psub = psub, pb = pb, B = B,
            zl = [solver.nominal[t].zl for t = 1:T],
            zu = [solver.nominal[t].zu for t = 1:T])
end

"""Cases with energy bounds: slack control, nu = 3, nc = 2."""
function solve_ddp_energy(T::Int; ps_lo = -Inf, ps_hi = Inf, pb_lo = -Inf, pb_hi = Inf,
                          B_lo, B_hi, tol = 1e-10)
    nx, nu = 2, 3          # u = (Psub, P_B, s)
    pL = PL6[1:T]; c = C6[1:T]

    dyn = Dynamics((x, u) -> [x[1] - u[2] * Δt, x[2] + 1.0], nx, nu)
    l = (x, u) -> lagrange(c, x[2]) * u[1] * Δt + C_B * u[2]^2 * Δt
    lN = (x, u) -> c[T] * u[1] * Δt + C_B * u[2]^2 * Δt +
                   W * (x[1] - u[2] * Δt - B_0)^2
    # c1: power balance.  c2: slack defn for the NEXT energy, B_{t+1} - B_lo - s = 0
    #     with s ∈ [0, B_hi - B_lo], which is exactly B_lo <= B_{t+1} <= B_hi.
    cons = EqualityConstraints(
        (x, u) -> [u[1] + u[2] - lagrange(pL, x[2]),
                   (x[1] - u[2] * Δt) - B_lo - u[3]], nx, nu)
    cl = ControlLimits(SVector{nu,T_}([ps_lo, pb_lo, 0.0]),
                       SVector{nu,T_}([ps_hi, pb_hi, B_hi - B_lo]))

    ocp = build_ocp(T, Objective(l, nx, nu), Objective(lN, nx, nu), dyn, cons, cl)
    solver = Solver(ocp; options = Options{T_}(verbose = false, optimality_tolerance = tol))
    x1 = SVector{nx,T_}([B_0, 1.0])
    ū = [SVector{nu,T_}([1.0, 0.0, (B_hi - B_lo) / 2]) for _ = 1:T]
    solve!(solver, x1, ū)
    x, u = get_trajectory(solver)
    psub = [u[t][1] for t = 1:T]; pb = [u[t][2] for t = 1:T]
    B = [x[t][1] for t = 1:T]; push!(B, B[T] - pb[T] * Δt)
    return (solver = solver, psub = psub, pb = pb, B = B, s = [u[t][3] for t = 1:T],
            zl = [solver.nominal[t].zl for t = 1:T],
            zu = [solver.nominal[t].zu for t = 1:T])
end

# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

function report(name, T, got, ref; energy = false)
    pL = PL6[1:T]; c = C6[1:T]
    J = sum(c .* got.psub .* Δt) + C_B * sum(got.pb .^ 2) * Δt + W * (got.B[T+1] - B_0)^2
    bal = maximum(abs(got.psub[t] + got.pb[t] - pL[t]) for t = 1:T)
    dyn = maximum(abs(got.B[t+1] - (got.B[t] - got.pb[t] * Δt)) for t = 1:T)
    sd = got.solver.data

    println("\n" * "="^78)
    println(name)
    println("="^78)
    @printf("iterations %d, termination %s (status %d), reg_last %.2e\n",
            sd.k, status_str(sd.status), sd.status, sd.reg_last)

    if ref === nothing
        println("reference: NO KKT POINT EXISTS — this bound set is infeasible.")
        @printf("FilterDDP objective %.12f, max power-balance residual %.3e\n", J, bal)
        return (J = J, gap = NaN, bal = bal, dyn = dyn, perr = NaN,
                status = sd.status, iters = sd.k)
    end
    @printf("objective  FilterDDP %.12f   reference %.12f   gap %.3e\n", J, ref.J, abs(J - ref.J))
    @printf("max |ΔPsub| %.3e   max |ΔP_B| %.3e   max |ΔB| %.3e\n",
            maximum(abs.(got.psub .- ref.psub)), maximum(abs.(got.pb .- ref.pb)),
            maximum(abs.(got.B .- ref.B)))
    @printf("max power-balance residual %.3e   max dynamics residual %.3e\n", bal, dyn)
    @printf("solver primal_inf %.3e  dual_inf %.3e  cs_inf %.3e  μ %.2e\n",
            sd.primal_inf, sd.dual_inf, sd.cs_inf_0, sd.μ)

    println("\n  t      Psub        Psub_ref        P_B         P_B_ref        B          B_ref")
    for t = 1:T
        @printf("%3d  % 10.6f  % 10.6f  % 10.6f  % 10.6f  % 10.6f  % 10.6f\n",
                t, got.psub[t], ref.psub[t], got.pb[t], ref.pb[t], got.B[t], ref.B[t])
    end
    @printf("%3d  %10s  %10s  %10s  %10s  % 10.6f  % 10.6f\n",
            T + 1, "-", "-", "-", "-", got.B[T+1], ref.B[T+1])

    println("\nactive set (reference, from KKT multipliers):")
    isempty(ref.active) ? println("  (none)") : foreach(s -> println("  $s"), ref.active)
    println("bound multipliers from FilterDDP (zl / zu per stage; large => active):")
    for t = 1:T
        @printf("  t=%d  zl = %-26s zu = %s\n", t,
                string(round.(Array(got.zl[t]), digits = 4)),
                string(round.(Array(got.zu[t]), digits = 4)))
    end
    energy && @printf("energy slack s = %s\n", string(round.(got.s, digits = 6)))
    return (J = J, gap = abs(J - ref.J), bal = bal, dyn = dyn,
            perr = maximum(abs.(got.psub .- ref.psub)), status = sd.status, iters = sd.k)
end

function main()
    summary = Vector{Any}()

    # ---- 6a: grid-power lower bound only ---------------------------------
    # There is no upper limit on substation import in this model; the only
    # grid-side restriction is Psub >= 0, i.e. no export upstream, as in the
    # tADMM paper. It never binds here (every Psub exceeds 1), which makes this
    # a check that an INACTIVE bound is identified as such and leaves the
    # equality-only solution untouched.
    cfg = (ps_lo = 0.0,)
    push!(summary, ("6a grid lower bound, T=3",
        report("6a — grid-power lower bound only, T = 3  (Psub >= 0, no export)", 3,
               solve_ddp_bounds(3; cfg...), reference(3; cfg...))))

    # ---- 6b: battery-power bounds only ----------------------------------
    cfg = (pb_lo = -0.5, pb_hi = 0.10)
    push!(summary, ("6b battery-power bounds, T=3",
        report("6b — battery-power bounds only, T = 3  (-0.5 <= P_B <= 0.10)", 3,
               solve_ddp_bounds(3; cfg...), reference(3; cfg...))))

    # ---- 6c: energy bounds only (slack reformulation) -------------------
    cfg = (B_lo = 1.98, B_hi = 2.05)
    push!(summary, ("6c energy bounds, T=3",
        report("6c — energy bounds only, T = 3  (1.98 <= B <= 2.05), slack reformulation", 3,
               solve_ddp_energy(3; cfg...), reference(3; cfg...); energy = true)))

    # ---- 6d: all three ---------------------------------------------------
    # With the substation import unlimited above, the balance is always
    # satisfiable, so the binding constraints come from the battery and energy
    # limits alone.
    cfg = (ps_lo = 0.0, pb_lo = -0.5, pb_hi = 0.10, B_lo = 1.98, B_hi = 2.05)
    push!(summary, ("6d all bounds, T=3",
        report("6d — all bounds together, T = 3", 3,
               solve_ddp_energy(3; cfg...), reference(3; cfg...); energy = true)))

    # ---- 6e: T = 6 -------------------------------------------------------
    cfg = (ps_lo = 0.0, pb_lo = -0.5, pb_hi = 0.30, B_lo = 1.85, B_hi = 2.10)
    push!(summary, ("6e all bounds, T=6",
        report("6e — all bounds together, T = 6", 6,
               solve_ddp_energy(6; cfg...), reference(6; cfg...); energy = true)))

    # ---- 6f: tolerance sweep on case 6d ----------------------------------
    println("\n" * "="^78)
    println("6f — tolerance sweep on case 6d")
    println("="^78)
    cfg = (ps_lo = 0.0, pb_lo = -0.5, pb_hi = 0.10, B_lo = 1.98, B_hi = 2.05)
    ref3 = reference(3; cfg...)
    @printf("%10s %6s %6s %-22s %11s %11s\n", "tol", "iters", "status", "termination", "|ΔJ|", "dual_inf")
    for tol in [1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-12]
        g = solve_ddp_energy(3; cfg..., tol = tol)
        c3 = C6[1:3]
        J = sum(c3 .* g.psub .* Δt) + C_B * sum(g.pb .^ 2) * Δt + W * (g.B[4] - B_0)^2
        @printf("%10.0e %6d %6d %-22s %11.3e %11.3e\n", tol, g.solver.data.k,
                g.solver.data.status, status_str(g.solver.data.status),
                abs(J - ref3.J), g.solver.data.dual_inf)
    end

    # ---- 6g: behaviour on a genuinely INFEASIBLE bound set ----------------
    println("\n" * "="^78)
    println("6g — genuinely infeasible bounds, T = 6")
    println("="^78)
    # With substation import unlimited above, the balance can always be met, so
    # infeasibility has to be built from the storage side: forbid charging and
    # then demand more stored energy than the device starts with.
    icfg = (ps_lo = 0.0, pb_lo = 0.0, pb_hi = 0.30, B_lo = 2.05, B_hi = 2.10)
    @printf("Elementary proof: P_B >= %.2f forbids charging, so B can only fall:\n", icfg.pb_lo)
    @printf("B^2 = B_0 - P_B^1 Δt <= B_0 = %.2f.  But B^{t+1} >= %.2f is required,\n",
            B_0, icfg.B_lo)
    @printf("and %.2f < %.2f, so no feasible trajectory exists.\n", B_0, icfg.B_lo)
    r = reference(6; icfg...)
    println("independent active-set enumeration: ",
            r === nothing ? "no KKT point exists (confirms infeasible)" : "FOUND a point (unexpected)")
    g = solve_ddp_energy(6; icfg...)
    sd = g.solver.data
    @printf("FilterDDP: iters %d, status %d (%s), primal_inf %.3e, dual_inf %.3e\n",
            sd.k, sd.status, status_str(sd.status), sd.primal_inf, sd.dual_inf)
    bal = maximum(abs(g.psub[t] + g.pb[t] - PL6[t]) for t = 1:6)
    @printf("max power-balance residual of the returned point: %.3e\n", bal)
    println("NOTE: FilterDDP has no feasibility restoration phase, so an infeasible")
    println("      problem surfaces only as a non-zero status plus a residual the")
    println("      caller must check. There is no 'infeasible' status code.")

    println("\n" * "="^78)
    println("STAGE 6 SUMMARY")
    println("="^78)
    @printf("%-32s %6s %6s %11s %11s %11s\n", "case", "iters", "status", "|ΔJ|", "max|ΔPsub|", "max bal")
    for (name, r) in summary
        @printf("%-32s %6d %6d %11.3e %11.3e %11.3e\n", name, r.iters, r.status, r.gap, r.perr, r.bal)
    end
end

main()
