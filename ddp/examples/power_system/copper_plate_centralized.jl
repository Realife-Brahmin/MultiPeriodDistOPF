# Centralized reference for the copper-plate storage-scheduling model.
#
# Solves the ORIGINAL problem eq:cp_all of paper/sections/copper_plate_model.tex
# directly in JuMP + Ipopt, in the paper's own variables (P_Subs^t, P_B^t, B^t),
# with no reduction and no reformulation. This is a genuinely independent check on
# the FilterDDP results, because it differs from the DDP model in three ways that
# each could hide a bug:
#
#   1. TRUE per-period data. FilterDDP has one stage cost and one constraint
#      function for the whole horizon, so it must recover c^t and p_L^t from a
#      Lagrange interpolant in the state (eq:ocp_interp). This model indexes the
#      data arrays directly. Agreement therefore validates the interpolation
#      trick itself, not just the solver.
#   2. NO slack reformulation. The energy bound B^min <= B^{t+1} <= B^max is a
#      state bound, which the DDP form cannot express and which the FilterDDP
#      model routes through an extra control s^t (eq:ocp_slack). Here it is
#      imposed directly on B.
#   3. B^{T+1} is an explicit variable, not a quantity reconstructed from the
#      stage-T state and control.
#
# The terminal condition is the same quadratic PENALTY used everywhere else, not
# a hard equality -- the point is to reproduce the problem FilterDDP is actually
# solving, so this deliberately does NOT "fix" that modeling choice.
#
# Run from the unified environment (see envs/ddp2026/README.md):
#     julia --startup-file=no --project=envs/ddp2026 \
#           ddp/examples/power_system/copper_plate_centralized.jl
#
# Writes ddp/results/copper_plate/centralized_reference.csv, which
# verify_against_centralized.jl then diffs against a live FilterDDP solve.

using JuMP
using Ipopt
using LinearAlgebra
using Printf

const T_ = Float64

# ---------------------------------------------------------------------------
# Data -- identical to copper_plate_battery_bounds.jl
# ---------------------------------------------------------------------------

const PL6 = T_[1.0, 1.6, 1.2, 1.8, 0.9, 1.4]        # demand p_L^t
const C6  = T_[0.30, 0.55, 0.40, 0.62, 0.25, 0.48]  # energy price c^t
const C_B = T_(0.5)                                 # throughput coefficient
const Δt  = T_(1.0)                                 # period length
const B_0 = T_(2.0)                                 # initial stored energy
const W   = T_(5.0)                                 # terminal-energy weight

# ---------------------------------------------------------------------------
# The centralized model, eq:cp_all verbatim
# ---------------------------------------------------------------------------

"""
    centralized(T; ps_lo, ps_hi, pb_lo, pb_hi, B_lo, B_hi, tol)

Solve eq:cp_all for a horizon of `T` periods. Omitted bounds default to the
base instance (no inequalities at all). Returns `nothing` when Ipopt proves the
bound set infeasible.
"""
function centralized(T::Int; ps_lo = -Inf, ps_hi = Inf, pb_lo = -Inf, pb_hi = Inf,
                     B_lo = -Inf, B_hi = Inf, tol = 1e-12, verbose = false)
    pL = PL6[1:T]
    c  = C6[1:T]

    m = Model(Ipopt.Optimizer)
    verbose || set_silent(m)
    set_optimizer_attribute(m, "tol", tol)
    set_optimizer_attribute(m, "constr_viol_tol", tol)
    set_optimizer_attribute(m, "acceptable_tol", tol)

    @variable(m, Psub[1:T])
    @variable(m, P_B[1:T])
    @variable(m, B[1:T+1])

    # eq:cp_obj -- linear energy cost + quadratic throughput charge + terminal penalty
    @objective(m, Min,
        sum((c[t] * Psub[t] + C_B * P_B[t]^2) * Δt for t = 1:T)
        + W * (B[T+1] - B[1])^2)

    # eq:cp_balance -- one nodal balance per period
    @constraint(m, balance[t = 1:T], Psub[t] + P_B[t] - pL[t] == 0)

    # eq:cp_soc_init and eq:cp_soc_dyn -- lossless reservoir, the only time coupling
    @constraint(m, soc_init, B[1] == B_0)
    @constraint(m, soc_dyn[t = 1:T], B[t+1] == B[t] - P_B[t] * Δt)

    # eq:cp_pg_lim / eq:cp_pb_lim / eq:cp_soc_lim -- applied only where finite, so
    # that the base instance really is the equality-constrained problem.
    for t = 1:T
        isfinite(ps_lo) && set_lower_bound(Psub[t], ps_lo)
        isfinite(ps_hi) && set_upper_bound(Psub[t], ps_hi)
        isfinite(pb_lo) && set_lower_bound(P_B[t], pb_lo)
        isfinite(pb_hi) && set_upper_bound(P_B[t], pb_hi)
        isfinite(B_lo)  && set_lower_bound(B[t+1], B_lo)
        isfinite(B_hi)  && set_upper_bound(B[t+1], B_hi)
    end

    optimize!(m)
    st = termination_status(m)
    if st != MOI.LOCALLY_SOLVED && st != MOI.OPTIMAL
        return (status = st, psub = nothing, pb = nothing, B = nothing, J = NaN)
    end

    psub = value.(Psub); pb = value.(P_B); Bv = value.(B)
    J = sum(c .* psub .* Δt) + C_B * sum(pb .^ 2) * Δt + W * (Bv[T+1] - Bv[1])^2
    return (status = st, psub = psub, pb = pb, B = Bv, J = J,
            active = active_bounds(psub, pb, Bv, T, ps_lo, ps_hi, pb_lo, pb_hi, B_lo, B_hi))
end

"""Bounds that are tight at the solution, by primal value -- reported so the
active set can be compared against FilterDDP's bound multipliers."""
function active_bounds(psub, pb, B, T, ps_lo, ps_hi, pb_lo, pb_hi, B_lo, B_hi; tol = 1e-7)
    tags = String[]
    for t = 1:T
        isfinite(ps_lo) && abs(psub[t] - ps_lo) < tol && push!(tags, "Psub[$t] >= $ps_lo")
        isfinite(ps_hi) && abs(psub[t] - ps_hi) < tol && push!(tags, "Psub[$t] <= $ps_hi")
        isfinite(pb_lo) && abs(pb[t] - pb_lo)   < tol && push!(tags, "P_B[$t] >= $pb_lo")
        isfinite(pb_hi) && abs(pb[t] - pb_hi)   < tol && push!(tags, "P_B[$t] <= $pb_hi")
        isfinite(B_lo)  && abs(B[t+1] - B_lo)   < tol && push!(tags, "B[$(t+1)] >= $B_lo")
        isfinite(B_hi)  && abs(B[t+1] - B_hi)   < tol && push!(tags, "B[$(t+1)] <= $B_hi")
    end
    return tags
end

# ---------------------------------------------------------------------------
# Closed form, eq:cp_closed_form -- valid for the BASE instance only
# ---------------------------------------------------------------------------

function closed_form(T::Int)
    pL = PL6[1:T]; c = C6[1:T]
    s  = sum(c) / (2 * C_B + 2 * T * W * Δt)
    pb = (c .- 2 * W * Δt * s) ./ (2 * C_B)
    psub = pL .- pb
    B = zeros(T_, T + 1); B[1] = B_0
    for t = 1:T
        B[t+1] = B[t] - pb[t] * Δt
    end
    J = sum(c .* psub .* Δt) + C_B * sum(pb .^ 2) * Δt + W * (B[T+1] - B_0)^2
    return (psub = psub, pb = pb, B = B, J = J)
end

# ---------------------------------------------------------------------------
# Cases -- the same set the FilterDDP scripts run
# ---------------------------------------------------------------------------

const CASES = [
    ("base_T3", 3, NamedTuple()),
    ("base_T6", 6, NamedTuple()),
    ("6a_T3",   3, (ps_lo = 0.0,)),
    ("6b_T3",   3, (pb_lo = -0.5, pb_hi = 0.10)),
    ("6c_T3",   3, (B_lo = 1.98, B_hi = 2.05)),
    ("6d_T3",   3, (ps_lo = 0.0, pb_lo = -0.5, pb_hi = 0.10, B_lo = 1.98, B_hi = 2.05)),
    ("6e_T6",   6, (ps_lo = 0.0, pb_lo = -0.5, pb_hi = 0.30, B_lo = 1.85, B_hi = 2.10)),
    ("6g_T6",   6, (ps_lo = 0.0, pb_lo = 0.0, pb_hi = 0.30, B_lo = 2.05, B_hi = 2.10)),
]

function main()
    outdir = joinpath(@__DIR__, "..", "..", "results", "copper_plate")
    mkpath(outdir)
    csv = joinpath(outdir, "centralized_reference.csv")
    rows = ["case,T,t,psub,pb,B,J"]

    println("="^78)
    println("CENTRALIZED REFERENCE — eq:cp_all solved directly in JuMP + Ipopt")
    println("="^78)

    for (name, T, cfg) in CASES
        r = centralized(T; cfg...)
        println("\n" * "-"^78)
        @printf("%-10s T = %d   Ipopt status: %s\n", name, T, string(r.status))
        println("-"^78)

        if r.psub === nothing
            println("  INFEASIBLE — Ipopt certifies no feasible point exists.")
            println("  (FilterDDP returns status 7 here, which does not distinguish")
            println("   infeasibility from a hard problem — see the Stage 6g note.)")
            push!(rows, "$name,$T,-1,NaN,NaN,NaN,NaN")
            continue
        end

        @printf("  objective J = %.14f\n", r.J)

        # Base instance only: check against the closed form of eq:cp_closed_form.
        if isempty(cfg)
            cf = closed_form(T)
            @printf("  closed form  = %.14f   |ΔJ| = %.3e\n", cf.J, abs(r.J - cf.J))
            @printf("  max |ΔPsub| = %.3e   max |ΔP_B| = %.3e   max |ΔB| = %.3e\n",
                    maximum(abs.(r.psub .- cf.psub)), maximum(abs.(r.pb .- cf.pb)),
                    maximum(abs.(r.B .- cf.B)))
        end

        bal = maximum(abs(r.psub[t] + r.pb[t] - PL6[t]) for t = 1:T)
        dyn = maximum(abs(r.B[t+1] - (r.B[t] - r.pb[t] * Δt)) for t = 1:T)
        @printf("  max balance residual %.3e   max dynamics residual %.3e\n", bal, dyn)

        println("\n    t       Psub          P_B            B")
        for t = 1:T
            @printf("  %3d  % 12.8f  % 12.8f  % 12.8f\n", t, r.psub[t], r.pb[t], r.B[t])
            push!(rows, "$name,$T,$t,$(r.psub[t]),$(r.pb[t]),$(r.B[t]),$(r.J)")
        end
        @printf("  %3d  %12s  %12s  % 12.8f\n", T + 1, "-", "-", r.B[T+1])
        push!(rows, "$name,$T,$(T+1),NaN,NaN,$(r.B[T+1]),$(r.J)")

        println("\n  active bounds at the solution:")
        isempty(r.active) ? println("    (none)") : foreach(s -> println("    $s"), r.active)
    end

    open(csv, "w") do io
        for row in rows
            println(io, row)
        end
    end
    println("\n" * "="^78)
    println("wrote $csv")
    println("="^78)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
