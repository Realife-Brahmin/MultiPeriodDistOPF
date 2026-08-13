# Stage 10: what FilterDDP actually forms, in this problem's own variables.
#
# The paper's optimal-control form is written in the solver's generic notation
# (x, u, l, f, c). This script pins every one of those objects down to the
# copper-plate quantities P_Subs^t, P_B^t, B^t and the KKT multipliers, and
# CHECKS the mapping numerically against the derivative functions FilterDDP
# actually builds -- so the equations written in the paper are verified, not
# asserted.
#
# Every claim below is evaluated at a random interior point of case 6d
# (all bounds, T = 6), against ocp.stage_objective.*, ocp.constraints.* and
# ocp.dynamics.* as constructed by build_ocp.
#
# Run:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/kkt_structure.jl

include("copper_plate_battery_bounds.jl")
using LinearAlgebra

const B_LO = 1.20
const B_HI = 1.95
const TT   = 6

# Build exactly the OCP that case 6d solves (nx = 2, nu = 3, nc = 2).
function build_case_ocp(T::Int; B_lo = B_LO, B_hi = B_HI,
                        ps_lo = 0.0, pb_lo = -0.45, pb_hi = 0.45)
    nx, nu = 2, 3
    pL = tadmm_pL(T); c = tadmm_cost(T)
    dyn  = Dynamics((x, u) -> [x[1] - u[2] * Δt, x[2] + 1.0], nx, nu)
    l    = (x, u) -> lagrange(c, x[2]) * u[1] * Δt + C_B * u[2]^2 * Δt
    lN   = (x, u) -> c[T] * u[1] * Δt + C_B * u[2]^2 * Δt +
                     W * (x[1] - u[2] * Δt - B_0)^2
    cons = EqualityConstraints(
        (x, u) -> [u[1] + u[2] - lagrange(pL, x[2]),
                   (x[1] - u[2] * Δt) - B_lo - u[3]], nx, nu)
    cl = ControlLimits(SVector{nu,T_}([ps_lo, pb_lo, 0.0]),
                       SVector{nu,T_}([Inf, pb_hi, B_hi - B_lo]))
    return build_ocp(T, Objective(l, nx, nu), Objective(lN, nx, nu), dyn, cons, cl)
end

ocp = build_case_ocp(TT)
pL  = tadmm_pL(TT); cvec = tadmm_cost(TT)

# An arbitrary interior evaluation point: B, τ, P_Subs, P_B, s.
x = SVector{2,T_}([1.63, 3.0])
u = SVector{3,T_}([1.71, 0.22, 0.19])
Bv, τ, Psub, PB, sv = x[1], x[2], u[1], u[2], u[3]

# Derivatives of the Lagrange interpolants at τ, by central differences --
# these are the only non-elementary pieces, and they enter through the τ column.
h  = 1e-5
dc  = (lagrange(cvec, τ + h) - lagrange(cvec, τ - h)) / (2h)
dpL = (lagrange(pL,   τ + h) - lagrange(pL,   τ - h)) / (2h)

checks = Tuple{String,Any,Any}[]
add(name, got, want) = push!(checks, (name, got, want))

# ---------------------------------------------------------------------------
# Dynamics:  B^{t+1} = B^t - P_B^t Δt,   τ^{t+1} = τ^t + 1
# ---------------------------------------------------------------------------
add("f_x  = I",              Array(ocp.dynamics.fx(x, u)), [1.0 0.0; 0.0 1.0])
add("f_u  = [0 -Δt 0; 0 0 0]", Array(ocp.dynamics.fu(x, u)), [0.0 -Δt 0.0; 0.0 0.0 0.0])
# Dynamics are linear, so every second derivative vanishes: the V_x-weighted
# f_xx / f_ux / f_uu contractions drop out of the backward pass entirely.
add("f_xx = 0", Array(ocp.dynamics.fxx(x, u, SVector{2,T_}([1.0, 1.0]))), zeros(2, 2))
add("f_uu = 0", Array(ocp.dynamics.fuu(x, u, SVector{2,T_}([1.0, 1.0]))), zeros(3, 3))
add("f_ux = 0", Array(ocp.dynamics.fux(x, u, SVector{2,T_}([1.0, 1.0]))), zeros(3, 2))

# ---------------------------------------------------------------------------
# Stage cost:  l = c(τ) P_Subs Δt + C_B P_B^2 Δt
# ---------------------------------------------------------------------------
add("l_x  = [0, c'(τ) P_Subs Δt]",
    Array(ocp.stage_objective.lx(x, u)), [0.0, dc * Psub * Δt])
add("l_u  = [c(τ)Δt, 2 C_B P_B Δt, 0]",
    Array(ocp.stage_objective.lu(x, u)),
    [lagrange(cvec, τ) * Δt, 2 * C_B * PB * Δt, 0.0])
add("l_uu = diag(0, 2 C_B Δt, 0)",
    Array(ocp.stage_objective.luu(x, u)), diagm([0.0, 2 * C_B * Δt, 0.0]))

# ---------------------------------------------------------------------------
# Constraints:  c1 = P_Subs + P_B - p_L(τ)
#               c2 = (B - P_B Δt) - B^min - s
# ---------------------------------------------------------------------------
add("c    = [balance, slack]", Array(ocp.constraints.c(x, u)),
    [Psub + PB - lagrange(pL, τ), (Bv - PB * Δt) - B_LO - sv])
add("c_u  = [1 1 0; 0 -Δt -1]", Array(ocp.constraints.cu(x, u)),
    [1.0 1.0 0.0; 0.0 -Δt -1.0])
add("c_x  = [0 -p_L'(τ); 1 0]", Array(ocp.constraints.cx(x, u)),
    [0.0 -dpL; 1.0 0.0])
add("c_uu = 0 (linear in u)",
    Array(ocp.constraints.cuu(x, u, SVector{2,T_}([0.7, -0.3]))), zeros(3, 3))
add("c_ux = 0",
    Array(ocp.constraints.cux(x, u, SVector{2,T_}([0.7, -0.3]))), zeros(3, 2))

# ---------------------------------------------------------------------------
# Terminal cost adds  w (B^{T+1} - B_0)^2  with  B^{T+1} = B^N - P_B^N Δt
# ---------------------------------------------------------------------------
Bend = Bv - PB * Δt
add("lN_u = [c^T Δt, 2C_B P_B Δt - 2wΔt(B^{T+1}-B_0), 0]",
    Array(ocp.term_objective.lu(x, u)),
    [cvec[TT] * Δt, 2 * C_B * PB * Δt - 2 * W * Δt * (Bend - B_0), 0.0])
add("lN_xx[1,1] = 2w", Array(ocp.term_objective.lxx(x, u))[1, 1], 2 * W)
add("lN_uu[2,2] = 2C_BΔt + 2wΔt²",
    Array(ocp.term_objective.luu(x, u))[2, 2], 2 * C_B * Δt + 2 * W * Δt^2)

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
println("="^78)
println("Stage 10 — FilterDDP's objects, pinned to the copper-plate variables")
println("="^78)
println("evaluation point:  B = $Bv, τ = $τ, P_Subs = $Psub, P_B = $PB, s = $sv\n")
worst = Ref(0.0)
for (name, got, want) in checks
    err = maximum(abs.(vec([got...]) .- vec([want...])))
    worst[] = max(worst[], err)
    @printf("  %-46s max err %.2e  %s\n", name, err, err < 1e-6 ? "ok" : "MISMATCH")
end
@printf("\nworst disagreement over all %d checks: %.2e\n", length(checks), worst[])
worst[] < 1e-6 || error("a hand-derived expression does not match the built OCP")

# ---------------------------------------------------------------------------
# The reduced system actually solved per stage
# ---------------------------------------------------------------------------
println("\n" * "="^78)
println("The per-stage reduced system")
println("="^78)
cu = Array(ocp.constraints.cu(x, u))
nullbasis = nullspace(cu)
@printf("n_u = 3, n_c = 2  ->  null space of c_u has dimension %d\n", size(nullbasis, 2))
d = vec(nullbasis) / vec(nullbasis)[2]        # normalise so the P_B entry is 1
@printf("free direction (normalised on P_B):  ΔP_Subs = %+.3f, ΔP_B = %+.3f, Δs = %+.3f\n",
        d[1], d[2], d[3])
println("  i.e. the ONLY freedom left at a stage is to move power between the")
println("  substation and the battery, with the energy slack following.")
@printf("residual ‖c_u d‖ = %.2e\n", norm(cu * d))

# The scalar that FilterDDP Choleskys, written out.
println("\nWith f_uu = c_uu = 0, the stage Hessian is")
println("  Ĥ = diag(0, 2 C_B Δt, 0) + diag(Σ_L + Σ_U) + f_u' V_BB f_u + reg·I")
println("so the reduced Hessian along the single free direction is the SCALAR")
println("  dᵀĤd = (Σ_L+Σ_U)_Psub + 2 C_B Δt + Δt² V_BB + (Σ_L+Σ_U)_P_B")
println("         + Δt² (Σ_L+Σ_U)_s ,        Σ_L = z_l/(u-l),  Σ_U = z_u/(u_hi-u)")
println("\nThe 'dense n_u x n_u KKT solve' is therefore a 1x1 division here, and")
println("the Cholesky that can fail is a positivity test on that one number.")
