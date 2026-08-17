# Stage 5: minimal copper-plate battery scheduling problem for FilterDDP.jl
#
# Equality-constrained only. No AC power flow, no LinDistFlow, no SOCP, one bus,
# one signed battery power variable, no bounds (bounds are Stage 6).
#
# ---------------------------------------------------------------------------
# Formulation  (matches ddp/paper/sections/copper_plate_model.tex eq. (1))
# ---------------------------------------------------------------------------
#
#   min   Σ_t c^t·Psub^t·Δt + C_B·(P_B^t)²·Δt + w·(B^{T+1} - B^1)²
#   s.t.  Psub^t + P_B^t - p_L^t = 0            t = 1..T     (power balance)
#         B^1 = B_0
#         B^{t+1} = B^t - P_B^t·Δt             t = 1..T     (storage, lossless)
#
#   Psub = grid import, P_B = signed battery power, B = stored energy.
#   Sign convention follows the tADMM formulation: P_B > 0 DISCHARGES.
#
# There is deliberately NO round-trip efficiency term. A single η applied to a
# signed P_B cannot be physically right in both directions: charging genuinely
# stores η_c·|P_B|, but discharging P_B to the bus must drain P_B/η_d from the
# device, and η and 1/η coincide only at η = 1. Representing losses honestly
# needs split P_c, P_d ≥ 0 with B^{t+1} = B^t + η_c·P_c·Δt − P_d·Δt/η_d. The
# tADMM paper models a lossless battery, so this benchmark does too, and the
# symbol is dropped rather than carried at 1.0 where it would invite someone to
# set it to 0.95 later and get a quietly wrong model.
#
# The objective is the one used in the IAS-Trans tADMM paper: a linear energy
# cost at the time-varying price c^t, plus a small quadratic charge on battery
# throughput. There is deliberately NO quadratic term on Psub — C_B alone makes
# the problem strictly convex, which is all the closed-form reference needs.
#
# Horizon T = 6, mapped onto FilterDDP stages t = 1..N with N = T.
#   state    x_t = (B_t)           the stored energy, and nothing else
#   control  u_t = (Psub_t, P_B_t)
#
# Note 1 — there is no time index in the state.  Earlier revisions carried one.
# FilterDDP shipped with ONE `Objective` and ONE `EqualityConstraints` for the
# whole horizon, so per-period data (p_L^t, c^t) had nowhere to live but the
# state, and was recovered from a Lagrange interpolant of a counter τ. That was
# never a limitation of the algorithm — Remark 1 of the global-convergence paper
# says the objective, dynamics and constraint functions "can in general be
# time-varying" — only of the implementation, so it is lifted by
# ddp/patches/per_stage_data.patch. Each stage now holds its own c^t and p_L^t,
# the interpolant is gone, and nx drops from 2 to 1. The interpolant is parked in
# tadmm_profiles.jl for the case where upstream cannot be patched.
#
# Note 2 — why the terminal energy target is a COST and not a CONSTRAINT.
# We would like B^{T+1} = B_0 as a hard terminal equality. FilterDDP cannot express
# it: `c` is the same function at every stage, so a terminal-only equality is not
# representable, and a constraint written to vanish at non-terminal stages would
# have ∇_u c = 0 there, making `lu(AY)` singular in the backward pass
# (src/backward_pass.jl:124-127). See notes/FILTERDDP_API.md §7. We therefore
# penalise (B^{T+1} - B_0)² with weight w in the TERMINAL cost. Because
# B^{T+1} = B^N - P_B^N·Δt and both B^N and P_B^N are available at stage N, this
# is an exact algebraic restatement, not an approximation of the dynamics.
#
# Note 3 — what w > 0 is for.  The tADMM paper has no terminal energy condition:
# its periods are coupled by the SOC *bounds* plus the dynamics. Bounds are
# inequalities, and this base instance is deliberately equality-only, so without w
# the problem would separate into T independent single-period problems with
# P_B^t = c^t/(2·C_B). That separation is demonstrated, not assumed: see the w = 0
# probe at the bottom of this file. Once Stage 6 adds SOC bounds, w can be dropped
# and the coupling comes from the bounds exactly as in the paper.

using FilterDDP
using LinearAlgebra
using Printf
using StaticArrays

const T_ = Float64

# ---------------------------------------------------------------------------
# Problem data (tiny, deterministic)
# ---------------------------------------------------------------------------

const N      = 6                       # = T, number of periods
const nx     = 1                       # (B)  -- no time index since the per-stage patch
const nu     = 2                       # (Psub, P_B)

# Demand and price are the tADMM profiles at T = 6; see tadmm_profiles.jl.
include("tadmm_profiles.jl")

const pL   = T_.(tadmm_pL(N))          # demand   [p.u.]
const cvec = T_.(tadmm_cost(N))        # energy price per period (> 0)
const C_B  = T_(0.05)                   # battery throughput coefficient (> 0)
const Δt   = T_(1.0)                   # period length
const B_0  = T_(2.0)                   # initial stored energy
const w    = T_(5.0)                   # terminal-energy weight (see Note 2/3)


# ---------------------------------------------------------------------------
# Reference solution 1: closed form
# ---------------------------------------------------------------------------
# Eliminating Psub^t = p_L^t - P_B^t and B^{T+1} - B^1 = -Δt·s (s = Σ P_B^t):
#   J(P_B) = const - Δt Σ c^t P_B^t + C_B Δt Σ (P_B^t)² + w Δt² s²
# Stationarity, divided through by Δt:  -c^k + 2 C_B P_B^k + 2 w Δt s = 0, so
#   P_B^k = (c^k - 2 w Δt s) / (2 C_B),   s = Σ c^t / (2 C_B + 2 T w Δt).

function reference_closed_form()
    s = sum(cvec) / (2 * C_B + 2 * N * w * Δt)
    pb = (cvec .- 2 * w * Δt * s) ./ (2 * C_B)
    psub = pL .- pb
    B = zeros(T_, N + 1)
    B[1] = B_0
    for t = 1:N
        B[t+1] = B[t] - pb[t] * Δt
    end
    J = sum(cvec .* psub .* Δt) + C_B * sum(pb .^ 2) * Δt + w * (B[N+1] - B_0)^2
    return (psub = psub, pb = pb, B = B, J = J, s = s)
end

# ---------------------------------------------------------------------------
# Reference solution 2: independent dense KKT solve of the same QP
# ---------------------------------------------------------------------------
# Built from scratch over z = [psub1,pb1,...,psubN,pbN,B1,...,B_{N+1}] with all
# constraints written explicitly, so it shares no algebra with reference 1.

function reference_kkt()
    n = 3N + 1                          # 2N controls + (N+1) stored energies
    ips(t) = 2t - 1
    ipb(t) = 2t
    iB(t) = 2N + t                      # B1..B_{N+1} -> 2N+1 .. 3N+1

    Q = zeros(T_, n, n)                 # J = ½ zᵀQz + qᵀz + const
    q = zeros(T_, n)
    for t = 1:N
        q[ips(t)] = cvec[t] * Δt        # linear energy cost
        Q[ipb(t), ipb(t)] = 2 * C_B * Δt  # battery throughput
    end
    Q[iB(N + 1), iB(N + 1)] = 2 * w
    q[iB(N + 1)] = -2 * w * B_0

    rows = Vector{Vector{T_}}()
    rhs = Vector{T_}()
    for t = 1:N                          # balance: psub_t + pb_t = p_L_t
        r = zeros(T_, n); r[ips(t)] = 1.0; r[ipb(t)] = 1.0
        push!(rows, r); push!(rhs, pL[t])
    end
    for t = 1:N                          # dynamics: B_{t+1} - B_t + pb_t Δt = 0
        r = zeros(T_, n); r[iB(t+1)] = 1.0; r[iB(t)] = -1.0; r[ipb(t)] = Δt
        push!(rows, r); push!(rhs, 0.0)
    end
    r = zeros(T_, n); r[iB(1)] = 1.0     # initial condition
    push!(rows, r); push!(rhs, B_0)

    A = reduce(vcat, [r' for r in rows])
    m = size(A, 1)
    K = [Q A'; A zeros(T_, m, m)]
    sol = K \ [-q; rhs]
    z = sol[1:n]
    J = 0.5 * z' * Q * z + q' * z + w * B_0^2   # + dropped constant
    return (psub = [z[ips(t)] for t = 1:N], pb = [z[ipb(t)] for t = 1:N],
            B = [z[iB(t)] for t = 1:N+1], J = J, λ = sol[n+1:end])
end

# ---------------------------------------------------------------------------
# FilterDDP model
# ---------------------------------------------------------------------------

function build_solver(; w_term::T_ = w, tol::T_ = 1e-10, verbose::Bool = false,
                        maxit::Int = 1000)
    # State is the stored energy alone. Each stage carries its own price and
    # demand (per_stage_data.patch), so nothing is interpolated.
    dyn = Dynamics((x, u) -> [x[1] - u[2] * Δt], nx, nu)

    objs = Any[Objective((x, u) -> cvec[t] * u[1] * Δt + C_B * u[2]^2 * Δt, nx, nu)
               for t = 1:N]
    cons = Any[EqualityConstraints((x, u) -> [u[1] + u[2] - pL[t]], nx, nu) for t = 1:N]
    lN = (x, u) -> cvec[N] * u[1] * Δt + C_B * u[2]^2 * Δt +
                   w_term * (x[1] - u[2] * Δt - B_0)^2
    term_obj = Objective(lN, nx, nu)

    # no bounds at Stage 5: ±Inf disables the barrier entirely (nl = nu = 0)
    cl = ControlLimits(SVector{nu,T_}(fill(-Inf, nu)), SVector{nu,T_}(fill(Inf, nu)))

    ocp = build_ocp(N, objs[1], term_obj, dyn, cons[1], cl;
                    stage_objectives = objs, stage_constraints = cons)
    options = Options{T_}(verbose = verbose, optimality_tolerance = tol,
                          max_iterations = maxit)
    return Solver(ocp; options = options)
end

const STATUS = Dict(0 => "Optimal solution found",
                    1 => "Backward pass failed (inertia)",
                    7 => "Line search failed in forward pass",
                    8 => "Maximum iterations reached")
status_str(s) = get(STATUS, s, "internal status $s")

"""Recompute every residual from the returned trajectory, independently of solver internals."""
function audit(solver, ref)
    x, u = get_trajectory(solver)
    psub = [u[t][1] for t = 1:N]
    pb = [u[t][2] for t = 1:N]
    B = [x[t][1] for t = 1:N]
    B_end = B[N] - pb[N] * Δt

    bal = maximum(abs(psub[t] + pb[t] - pL[t]) for t = 1:N)
    dyn_res = N == 1 ? 0.0 :
        maximum(abs(B[t+1] - (B[t] - pb[t] * Δt)) for t = 1:N-1)
    J = sum(cvec .* psub .* Δt) + C_B * sum(pb .^ 2) * Δt + w * (B_end - B_0)^2

    return (psub = psub, pb = pb, B = B, B_end = B_end, bal = bal, dyn = dyn_res,
            J = J, term = abs(B_end - B_0),
            gap = abs(J - ref.J), psub_err = maximum(abs.(psub .- ref.psub)),
            pb_err = maximum(abs.(pb .- ref.pb)))
end

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

function main()
    ref1 = reference_closed_form()
    ref2 = reference_kkt()

    println("="^78)
    println("Stage 5 — copper-plate battery, T = 6, equality constraints only")
    println("="^78)
    @printf("data: p_L = %s\n      c   = %s\n", pL, cvec)
    @printf("      C_B = %.3f, B_0 = %.3f, Δt = %.3f, w = %.3f  (lossless battery)\n\n",
            C_B, B_0, Δt, w)

    println("--- independent reference solutions ---")
    fmt(v) = "[" * join((@sprintf("% .8f", x) for x in v), " ") * "]"
    @printf("closed form : Psub = %s  J = %.12f\n", fmt(ref1.psub), ref1.J)
    @printf("dense KKT   : Psub = %s  J = %.12f\n", fmt(ref2.psub), ref2.J)
    @printf("agreement between the two references: max|ΔPsub| = %.3e, |ΔJ| = %.3e\n\n",
            maximum(abs.(ref1.psub .- ref2.psub)), abs(ref1.J - ref2.J))

    # --- three different initial control trajectories -----------------------
    inits = Dict(
        "zeros"       => [SVector{nu,T_}([0.0, 0.0]) for _ = 1:N],
        "flat-demand" => [SVector{nu,T_}([1.0, 0.0]) for _ = 1:N],
        "bad/far"     => [SVector{nu,T_}([-5.0, 8.0]) for _ = 1:N],
    )

    println("--- FilterDDP solves from three initial control trajectories ---")
    @printf("%-14s %6s %6s %-28s %12s %11s %11s %11s\n",
            "init ū", "iters", "status", "termination", "objective",
            "|Δ J|", "max bal", "max dyn")
    for name in ["zeros", "flat-demand", "bad/far"]
        solver = build_solver()
        x1 = SVector{nx,T_}([B_0])
        solve!(solver, x1, inits[name])
        A = audit(solver, ref1)
        d = solver.data
        @printf("%-14s %6d %6d %-28s %12.6f %11.3e %11.3e %11.3e\n",
                name, d.k, d.status, status_str(d.status), A.J, A.gap, A.bal, A.dyn)
    end

    # --- detailed report for the nominal run --------------------------------
    solver = build_solver()
    x1 = SVector{nx,T_}([B_0])
    solve!(solver, x1, inits["zeros"])
    A = audit(solver, ref1)
    d = solver.data

    println("\n--- detailed result (init = zeros) ---")
    @printf("iterations            : %d\n", d.k)
    @printf("termination           : %s (status %d)\n", status_str(d.status), d.status)
    @printf("objective (solver)    : %.12f\n", d.objective)
    @printf("objective (recomputed): %.12f\n", A.J)
    @printf("objective (reference) : %.12f\n", ref1.J)
    @printf("objective gap         : %.3e\n", A.gap)
    @printf("primal inf (solver)   : %.3e\n", d.primal_inf)
    @printf("dual inf   (solver)   : %.3e\n", d.dual_inf)
    @printf("max power-balance res : %.3e\n", A.bal)
    @printf("max dynamics res      : %.3e\n", A.dyn)
    @printf("terminal energy res   : %.3e   (B_end = %.10f, target %.10f)\n",
            A.term, A.B_end, B_0)
    @printf("max |Psub - Psub_ref| : %.3e\n", A.psub_err)
    @printf("max |P_B  - P_B_ref|  : %.3e\n", A.pb_err)
    @printf("final regularisation  : %.3e\n", d.reg_last)

    println("\n  t       Psub          P_B            B        Psub_ref      P_B_ref")
    for t = 1:N
        @printf("%3d  % 11.7f  % 11.7f  % 11.7f  % 11.7f  % 11.7f\n",
                t, A.psub[t], A.pb[t], A.B[t], ref1.psub[t], ref1.pb[t])
    end
    @printf("%3d  %11s  %11s  % 11.7f  %11s  % 11.7f\n",
            N + 1, "-", "-", A.B_end, "-", ref1.B[N+1])

    # --- probe: w = 0 decouples the periods ---------------------------------
    println("\n--- probe: w = 0 (terminal energy term removed) ---")
    println("With w = 0 the horizon separates into T independent single-period problems,")
    println("each with the unique minimiser P_B^t = c^t / (2·C_B).")
    pb_sep = cvec ./ (2 * C_B)
    @printf("predicted P_B = %s\n", string(round.(pb_sep, digits = 8)))
    s0 = build_solver(w_term = 0.0)
    solve!(s0, SVector{nx,T_}([B_0]), inits["zeros"])
    _, u0 = get_trajectory(s0)
    pb0 = [u0[t][2] for t = 1:N]
    @printf("FilterDDP P_B = %s   (iters %d, status %d)\n",
            string(round.(pb0, digits = 8)), s0.data.k, s0.data.status)
    @printf("max deviation from the separable prediction: %.3e\n",
            maximum(abs.(pb0 .- pb_sep)))
    println("=> confirms w > 0 is what couples the periods in the equality-only instance.")
end

main()
