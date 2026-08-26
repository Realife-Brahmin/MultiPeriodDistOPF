# The copper-plate model built with PER-STAGE data.
#
# This is the model as it should have been written all along. Every earlier
# version had to carry a time index tau in the state and recover c^t and p_L^t
# from a Lagrange interpolant of it, because the package exposed one stage cost
# and one constraint function for the whole horizon. That restriction was never
# in the ALGORITHM -- Remark 1 of the global-convergence paper says the
# objective, dynamics and constraint functions "can in general be time-varying"
# -- only in the implementation. ddp/patches/per_stage_data.patch lifts it.
#
# What that removes:
#   * the interpolant (eq:ocp_interp), and with it the degree-(T-1) expression
#     that stack-overflowed Symbolics at T = 48
#   * the time index tau, so the state is just the stored energy: nx = 2 -> 1
#   * every derivative with respect to tau, which was never read anyway
#
# The parked interpolant lives in tadmm_profiles.jl for the case where the
# upstream package cannot be patched.
#
# Run:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/per_stage_model.jl

using FilterDDP
using LinearAlgebra
using Printf
using StaticArrays

const T_ = Float64
include("tadmm_profiles.jl")

const C_B = T_(0.05)
const Δt  = T_(1.0)
const B_0 = T_(2.0)
const W   = T_(5.0)

# ---------------------------------------------------------------------------
# The model, in per-stage form
# ---------------------------------------------------------------------------
# state    x = (B,)                       nx = 1
# control  u = (P_Subs, P_B)              nu = 2      (+ slack s when B is bounded)
# stage t carries its OWN c^t and p_L^t; nothing is interpolated.

"""
Build the copper-plate OCP at horizon `T` with genuinely per-stage cost and
constraint functions. `B_lo`/`B_hi` finite adds the energy slack control.
"""
function build_per_stage(T::Int; ps_lo = -Inf, ps_hi = Inf,
                         pb_lo = -Inf, pb_hi = Inf, B_lo = -Inf, B_hi = Inf)
    pL = tadmm_pL(T); c = tadmm_cost(T)
    energy = isfinite(B_lo) && isfinite(B_hi)
    nx = 1
    nu = energy ? 3 : 2

    dyn = Dynamics((x, u) -> [x[1] - u[2] * Δt], nx, nu)

    # One Objective and one EqualityConstraints PER STAGE, each closing over
    # that stage's own price and demand. No tau, no interpolation.
    stage_objs = Any[]
    stage_cons = Any[]
    for t = 1:T
        ct, pLt = c[t], pL[t]
        push!(stage_objs, Objective((x, u) -> ct * u[1] * Δt + C_B * u[2]^2 * Δt, nx, nu))
        push!(stage_cons, energy ?
            EqualityConstraints((x, u) -> [u[1] + u[2] - pLt,
                                           (x[1] - u[2] * Δt) - B_lo - u[3]], nx, nu) :
            EqualityConstraints((x, u) -> [u[1] + u[2] - pLt], nx, nu))
    end

    # Terminal stage: its own price plus the terminal-energy penalty.
    cT = c[T]
    lN = (x, u) -> cT * u[1] * Δt + C_B * u[2]^2 * Δt +
                   W * (x[1] - u[2] * Δt - B_0)^2

    cl = energy ?
        ControlLimits(SVector{nu,T_}([ps_lo, pb_lo, 0.0]),
                      SVector{nu,T_}([ps_hi, pb_hi, B_hi - B_lo])) :
        ControlLimits(SVector{nu,T_}([ps_lo, pb_lo]), SVector{nu,T_}([ps_hi, pb_hi]))

    return build_ocp(T, stage_objs[1], Objective(lN, nx, nu), dyn, stage_cons[1], cl;
                     stage_objectives = stage_objs, stage_constraints = stage_cons), nx, nu, energy
end

function solve_per_stage(T::Int; tol = 1e-10, kwargs...)
    ocp, nx, nu, energy = build_per_stage(T; kwargs...)
    solver = Solver(ocp; options = Options{T_}(verbose = false, optimality_tolerance = tol))
    u0 = energy ? [0.0, 0.0, 0.0] : [1.0, 0.0]
    ū = [SVector{nu,T_}(T_.(u0)) for _ = 1:T]
    solve!(solver, SVector{nx,T_}([B_0]), ū)
    x, u = get_trajectory(solver)
    psub = [u[t][1] for t = 1:T]; pb = [u[t][2] for t = 1:T]
    B = [x[t][1] for t = 1:T]; push!(B, B[T] - pb[T] * Δt)
    return (solver = solver, psub = psub, pb = pb, B = B)
end

"""Closed form of the base (equality-only) instance, eq:cp_closed_form."""
function closed_form(T::Int)
    pL = tadmm_pL(T); c = tadmm_cost(T)
    s  = sum(c) / (2 * C_B + 2 * T * W * Δt)
    pb = (c .- 2 * W * Δt * s) ./ (2 * C_B)
    B = zeros(T_, T + 1); B[1] = B_0
    for t = 1:T
        B[t+1] = B[t] - pb[t] * Δt
    end
    J = sum(c .* (pL .- pb) .* Δt) + C_B * sum(pb .^ 2) * Δt + W * (B[T+1] - B_0)^2
    return (pb = pb, J = J)
end

function main()
    println("="^84)
    println("Per-stage data: no interpolant, no time index, nx = 1")
    println("="^84)
    @printf("%5s %8s %7s %7s %14s %14s\n",
            "T", "build", "iters", "status", "|dJ| vs CF", "max|dP_B| vs CF")
    println("-"^84)
    for T in [6, 24, 48, 96, 192]
        ref = closed_form(T)
        t0 = time()
        local g, J
        try
            g = solve_per_stage(T)
            c = tadmm_cost(T)
            J = sum(c .* g.psub .* Δt) + C_B * sum(g.pb .^ 2) * Δt +
                W * (g.B[T+1] - B_0)^2
            @printf("%5d %7.0fs %7d %7d %14.3e %14.3e\n",
                    T, time() - t0, g.solver.data.k, g.solver.data.status,
                    abs(J - ref.J), maximum(abs.(g.pb .- ref.pb)))
        catch err
            @printf("%5d %7.0fs %7s %7s %14s %14s   <-- %s\n",
                    T, time() - t0, "ERR", "-", "-", "-",
                    err isa StackOverflowError ? "StackOverflowError" : string(typeof(err)))
        end
        flush(stdout)
    end
    println("\nCompare against the interpolant form, which stack-overflowed at T = 48.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
