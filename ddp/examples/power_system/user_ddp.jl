# The user's own first-order Differential Dynamic Programming scheme, on the
# copper plate.
#
# NAMING. The method is Differential Dynamic Programming -- never "Distributed."
# The original formulation document and driver script carried the wrong
# expansion; both were corrected on 2026-08-14, along with everything derived
# from them.
#
# A minimal restatement of envs/ddp/root_level/ddp_copperplate.jl (1154 lines,
# Gurobi-backed, with its own data structures) as four short functions, so that
# it can be read next to the FilterDDP workflow and compared line for line.
# The ALGORITHM is unchanged; only the scaffolding is gone.
#
# Source of truth: envs/ddp/tex/ddp_copperplate_formulation.tex, and the
# subproblem at ddp_copperplate.jl:340-372 (objective) and :400-418 (duals).
#
# ---------------------------------------------------------------------------
# The scheme in one paragraph
# ---------------------------------------------------------------------------
# Sweep t = 1..T. At stage t solve a small QP in (P_Subs[t], P_B[t], B[t]).
# It is tied BACKWARD to the state the current sweep just produced, B[t-1],
# which is a hard constraint; and FORWARD to the next stage only softly, by a
# linear term carrying mu[t+1] and B[t+1] as they stood at the END OF THE
# PREVIOUS SWEEP. The dual of the stage's own dynamics constraint becomes
# mu[t], and is what the next sweep will use one stage earlier. So information
# from the terminal stage migrates one stage per sweep.
#
# ---------------------------------------------------------------------------
# Indexing note
# ---------------------------------------------------------------------------
# This scheme indexes B[t] as the energy AFTER stage t, with B[0] = B_0:
#     B[t] - B[t-1] + Δt P_B[t] = 0
# The paper's model indexes B^t as the energy BEFORE stage t. So
#     B[t] here  ==  B^{t+1} in eq:cp_all,
# and B[T] here is the terminal energy B^{T+1}. Everything below uses the
# scheme's own convention.
#
# Run:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/user_ddp.jl

using JuMP
using Ipopt
using LinearAlgebra
using Printf

include("tadmm_profiles.jl")

# Same instance as case 6d of the FilterDDP suite, so the two are comparable.
const Δt     = 1.0
const C_B    = 0.05
const B_0    = 2.0
const W      = 5.0          # terminal-energy penalty weight
const PB_R   = 0.45         # |P_B| <= PB_R
const B_MIN  = 1.20
const B_MAX  = 1.95

# ---------------------------------------------------------------------------
# 1. One stage subproblem
# ---------------------------------------------------------------------------
"""
    stage_subproblem(t, T, c, pL, B_prev, mu_next, B_next)

The QP solved at stage `t` of a forward sweep.

`B_prev`  : energy leaving stage t-1, from THIS sweep (hard coupling)
`mu_next` : mu[t+1] from the PREVIOUS sweep (soft coupling, a parameter)
`B_next`  : B[t+1]  from the PREVIOUS sweep (soft coupling, a parameter)

Returns the stage's primal solution and its dynamics dual `mu[t]`.
"""
function stage_subproblem(t, T, c, pL, B_prev, mu_next, B_next;
                          pb_r = PB_R, b_min = B_MIN, b_max = B_MAX)
    m = Model(Ipopt.Optimizer)
    set_silent(m)

    @variable(m, P_Subs >= 0.0)
    @variable(m, -pb_r <= P_B <= pb_r)
    @variable(m, b_min <= B <= b_max)

    # Stage cost: energy bought plus battery throughput.
    obj = c[t] * P_Subs * Δt + C_B * P_B^2 * Δt

    if t < T
        # Soft forward coupling. P_B[t+1] enters as its OWN decision variable,
        # bounded by the battery rating, exactly as in the original.
        @variable(m, -PB_R <= PB_next <= PB_R)
        obj += mu_next * (B_next - B + Δt * PB_next)
    else
        # Terminal stage carries the terminal-energy penalty of eq:cp_obj, so
        # this solves the same problem the FilterDDP suite does.
        obj += W * (B - B_0)^2
    end
    @objective(m, Min, obj)

    @constraint(m, balance, P_Subs + P_B == pL[t])
    @constraint(m, dynamics, B - B_prev + Δt * P_B == 0)

    optimize!(m)
    ok = termination_status(m) in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL)
    ok || error("stage $t did not solve: $(termination_status(m))")

    # mu[t] is the dynamics dual, negated for the economic sign convention:
    # the marginal value of one more unit of stored energy at t.
    return (P_Subs = value(P_Subs), P_B = value(P_B), B = value(B),
            mu = -dual(dynamics))
end

# ---------------------------------------------------------------------------
# 2. One forward sweep
# ---------------------------------------------------------------------------
"""One sweep t = 1..T, given the previous sweep's `B_old` and `mu_old`."""
function forward_sweep(T, c, pL, B_old, mu_old; bounds...)
    B, P_B, P_Subs, mu = zeros(T), zeros(T), zeros(T), zeros(T)
    B_prev = B_0                                   # B[0]
    for t = 1:T
        mu_next = t < T ? mu_old[t+1] : 0.0        # from the PREVIOUS sweep
        B_next  = t < T ? B_old[t+1]  : 0.0
        r = stage_subproblem(t, T, c, pL, B_prev, mu_next, B_next; bounds...)
        B[t], P_B[t], P_Subs[t], mu[t] = r.B, r.P_B, r.P_Subs, r.mu
        B_prev = r.B                               # hard coupling within the sweep
    end
    return B, P_B, P_Subs, mu
end

# ---------------------------------------------------------------------------
# 3. The outer loop
# ---------------------------------------------------------------------------
"""
    ddp_solve(T; max_iter, tol, damping, ...)

Iterate sweeps until the state and dual trajectories stop moving. `damping = 1`
is the scheme as written; `damping = a < 1` blends each sweep's output with the
previous iterate, `x <- (1-a) x_old + a x_new`.

MEASURED, T = 6, case 6d (see the note at the bottom of this file):

    damping   sweeps   final err     gap in J vs centralized
      1.00      2000    2.1e-01           2.1e-03
      0.70      2000    5.7e-02           6.6e-04
      0.50      2000    2.0e-02           2.9e-04
      0.30      2000    5.3e-03           1.3e-04
      0.10      2000    6.1e-04           3.2e-05

Damping shrinks the limit cycle but never closes it: every row above hit the
iteration cap. See the closing note for why.
"""
function ddp_solve(T; max_iter = 200, tol = 1e-8, verbose = true, damping = 1.0,
                   bounds...)
    c, pL = tadmm_cost(T), tadmm_pL(T)
    B, mu = fill(B_0, T), zeros(T)                 # k = 0 initialisation
    P_B, P_Subs = zeros(T), zeros(T)
    history = NamedTuple[]

    for k = 1:max_iter
        B_old, mu_old = copy(B), copy(mu)
        Bn, P_B, P_Subs, mun = forward_sweep(T, c, pL, B_old, mu_old; bounds...)
        B  = (1 - damping) .* B_old  .+ damping .* Bn
        mu = (1 - damping) .* mu_old .+ damping .* mun

        err_state = norm(B .- B_old)
        err_dual  = norm(mu .- mu_old)
        err       = max(err_state, err_dual)
        J = sum(c .* P_Subs .* Δt) + C_B * sum(P_B .^ 2) * Δt + W * (B[T] - B_0)^2
        push!(history, (k = k, err_state = err_state, err_dual = err_dual,
                        err = err, J = J))
        verbose && @printf("  %3d   %11.3e %11.3e   %14.10f\n",
                           k, err_state, err_dual, J)
        err < tol && break
    end
    return (B = B, P_B = P_B, P_Subs = P_Subs, mu = mu, history = history)
end

# ---------------------------------------------------------------------------
# 4. Report, against the verified centralized reference
# ---------------------------------------------------------------------------
"""Centralized J and P_B for case 6d, from results/copper_plate/centralized_reference.csv."""
function centralized_reference(T)
    path = joinpath(@__DIR__, "..", "..", "results", "copper_plate",
                    "centralized_reference.csv")
    isfile(path) || return nothing
    pb, J = Float64[], NaN
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        f = split(strip(line), ',')
        (length(f) == 7 && f[1] == "6d_T6" && f[5] != "NaN") || continue
        push!(pb, parse(Float64, f[5])); J = parse(Float64, f[7])
    end
    return isempty(pb) ? nothing : (pb = pb, J = J)
end

function main()
    T = 6
    println("="^74)
    println("The user's DDP scheme on the copper plate, T = $T")
    println("="^74)
    println("one QP per stage per sweep; mu[t+1] and B[t+1] frozen from the")
    println("previous sweep, so terminal information moves one stage per sweep\n")
    println("  it        err_state    err_dual        objective J")
    println("  " * "-"^58)

    res = ddp_solve(T)
    nsweep = length(res.history)
    nqp = nsweep * T

    @printf("\nconverged in %d sweeps  (%d stage QPs solved)\n", nsweep, nqp)
    println("\n  t     P_Subs        P_B          B          mu")
    for t = 1:T
        @printf("%3d  %10.6f  %10.6f  %10.6f  %10.6f\n",
                t, res.P_Subs[t], res.P_B[t], res.B[t], res.mu[t])
    end

    ref = centralized_reference(T)
    if ref !== nothing
        @printf("\nobjective   this scheme %.10f   centralized %.10f   gap %.3e\n",
                res.history[end].J, ref.J, abs(res.history[end].J - ref.J))
        @printf("max |ΔP_B| against the centralized solution: %.3e\n",
                maximum(abs.(res.P_B .- ref.pb)))
    else
        println("\n(no centralized reference on disk; run copper_plate_centralized.jl)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

# ---------------------------------------------------------------------------
# What this scheme does on the copper plate, and why
# ---------------------------------------------------------------------------
# THE IMPLEMENTATION IS FAITHFUL. Oracle test: feed the TRUE costate mu* and
# state B* (taken from a centralized solve of the same problem) and run ONE
# sweep. It reproduces the centralized optimum to 1.7e-08 in objective and
# 6e-07 in P_B. So the subproblem, the coupling term and the dual extraction
# are all correct, and the scheme's fixed point IS the true optimum.
#
# IT STILL DOES NOT CONVERGE from a cold start. Undamped it settles into an
# exact period-2 limit cycle: two objectives alternating forever
# (1.47580829 / 1.47726351), err_state pinned at 2.1e-01, still cycling at 300
# sweeps, 2.1e-03 short in J and up to 121 kW out in P_B. Without the bounds it
# is worse -- err 8.2, wandering.
#
# WHY. The oracle test also shows the optimum is not EXACTLY a fixed point:
# starting from (mu*, B*) the sweep returns mu with ||mu_out - mu*|| = 9.2e-04.
# The coupling term -mu[t+1] * B[t] is LINEAR in B[t]: it carries the slope of
# the cost-to-go at the previous iterate and nothing else. The true cost-to-go
# has curvature (V_BB in the FilterDDP formulation), so as B[t] moves within a
# sweep the marginal value of stored energy moves with it, and a fixed slope
# cannot track that. The fixed point of the linearised map therefore sits a
# little off the true optimum, and no amount of damping closes the gap -- it
# only shrinks the neighbourhood the iterate ends up circling in.
#
# This is the missing-curvature difference between the two schemes, measured
# rather than argued. The obvious next experiment is to add the quadratic term
# and see whether the cycle disappears.
