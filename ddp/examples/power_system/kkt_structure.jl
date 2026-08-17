# Stage 10: what FilterDDP actually forms, in this problem's own variables.
#
# The paper's optimal-control form is written in the solver's generic notation
# (x, u, l, f, c). This script pins every one of those objects down to the
# copper-plate quantities P_Subs^t, P_B^t, B^t and the KKT multipliers, and
# CHECKS the mapping numerically against the derivative functions FilterDDP
# actually builds -- so the equations in the paper are verified, not asserted.
#
# Since ddp/patches/per_stage_data.patch the state is the stored energy alone
# (nx = 1): each stage carries its own c^t and p_L^t, so there is no time index
# and no interpolant, and every derivative with respect to a time index is gone.
#
# Everything below is evaluated at an interior point of case 6d (all bounds,
# T = 6), against the Objective / EqualityConstraints / Dynamics that build_ocp
# actually constructed.
#
# Run:
#   julia --startup-file=no --project=envs/ddp2026 \
#         ddp/examples/power_system/kkt_structure.jl

include("copper_plate_battery_bounds.jl")
using LinearAlgebra

const B_LO   = 1.20
const B_HI   = 1.95
const TT     = 6
const TSTAGE = 3            # the stage whose functions we inspect

ocp, _, _ = build_per_stage_ocp(TT; ps_lo = 0.0, pb_lo = -0.45, pb_hi = 0.45,
                                B_lo = B_LO, B_hi = B_HI)
pL   = tadmm_pL(TT)
cvec = tadmm_cost(TT)

# Interior evaluation point: B, then P_Subs, P_B, s.
x = SVector{1,T_}([1.63])
u = SVector{3,T_}([1.71, 0.22, 0.19])
Bv, Psub, PB, sv = x[1], u[1], u[2], u[3]

obj_t = FilterDDP.stage_obj(ocp, TSTAGE)
con_t = FilterDDP.stage_con(ocp, TSTAGE)
ct, pLt = cvec[TSTAGE], pL[TSTAGE]

checks = Tuple{String,Any,Any}[]
add(name, got, want) = push!(checks, (name, got, want))

# ---------------------------------------------------------------------------
# Dynamics:  B^{t+1} = B^t - P_B^t Δt
# ---------------------------------------------------------------------------
add("f_x  = [1]",         Array(ocp.dynamics.fx(x, u)), [1.0])
add("f_u  = [0  -Δt  0]", Array(ocp.dynamics.fu(x, u)), [0.0 -Δt 0.0])
# Linear dynamics: every second derivative vanishes, so the value-weighted
# f_xx / f_ux / f_uu contractions drop out of the backward pass entirely.
add("f_xx = 0", Array(ocp.dynamics.fxx(x, u, SVector{1,T_}([1.0]))), zeros(1, 1))
add("f_uu = 0", Array(ocp.dynamics.fuu(x, u, SVector{1,T_}([1.0]))), zeros(3, 3))
add("f_ux = 0", Array(ocp.dynamics.fux(x, u, SVector{1,T_}([1.0]))), zeros(3, 1))

# ---------------------------------------------------------------------------
# Stage cost:  l^t = c^t P_Subs Δt + C_B P_B^2 Δt      (c^t is now a NUMBER)
# ---------------------------------------------------------------------------
add("l_x  = [0] (cost does not depend on B)", Array(obj_t.lx(x, u)), [0.0])
add("l_u  = [c^t Δt, 2 C_B P_B Δt, 0]", Array(obj_t.lu(x, u)),
    [ct * Δt, 2 * C_B * PB * Δt, 0.0])
add("l_uu = diag(0, 2 C_B Δt, 0)", Array(obj_t.luu(x, u)),
    diagm([0.0, 2 * C_B * Δt, 0.0]))
add("l_ux = 0", Array(obj_t.lux(x, u)), zeros(3, 1))

# ---------------------------------------------------------------------------
# Constraints:  c1 = P_Subs + P_B - p_L^t
#               c2 = (B - P_B Δt) - B^min - s
# ---------------------------------------------------------------------------
add("c    = [balance, slack]", Array(con_t.c(x, u)),
    [Psub + PB - pLt, (Bv - PB * Δt) - B_LO - sv])
add("c_u  = [1 1 0; 0 -Δt -1]", Array(con_t.cu(x, u)), [1.0 1.0 0.0; 0.0 -Δt -1.0])
add("c_x  = [0; 1]", Array(con_t.cx(x, u)), reshape([0.0, 1.0], 2, 1))
add("c_uu = 0 (linear in u)",
    Array(con_t.cuu(x, u, SVector{2,T_}([0.7, -0.3]))), zeros(3, 3))
add("c_ux = 0", Array(con_t.cux(x, u, SVector{2,T_}([0.7, -0.3]))), zeros(3, 1))

# ---------------------------------------------------------------------------
# Terminal cost adds  w (B^{T+1} - B_0)^2  with  B^{T+1} = B^N - P_B^N Δt
# ---------------------------------------------------------------------------
Bend = Bv - PB * Δt
add("lN_u = [c^T Δt, 2C_B P_B Δt - 2wΔt(B^{T+1}-B_0), 0]",
    Array(ocp.term_objective.lu(x, u)),
    [cvec[TT] * Δt, 2 * C_B * PB * Δt - 2 * W * Δt * (Bend - B_0), 0.0])
add("lN_x  = 2w (B^{T+1} - B_0)", Array(ocp.term_objective.lx(x, u)),
    [2 * W * (Bend - B_0)])
add("lN_xx = 2w", Array(ocp.term_objective.lxx(x, u)), [2 * W])
add("lN_uu[2,2] = 2 C_B Δt + 2wΔt²",
    Array(ocp.term_objective.luu(x, u))[2, 2], 2 * C_B * Δt + 2 * W * Δt^2)

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
println("="^78)
println("Stage 10 — FilterDDP's objects, pinned to the copper-plate variables")
println("="^78)
@printf("state is B alone (nx = 1); stage %d carries c^t = %.4f, p_L^t = %.4f\n",
        TSTAGE, ct, pLt)
@printf("evaluation point:  B = %.2f, P_Subs = %.2f, P_B = %.2f, s = %.2f\n\n",
        Bv, Psub, PB, sv)

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
cu = Array(con_t.cu(x, u))
nb = nullspace(cu)
@printf("n_u = 3, n_c = 2  ->  null space of c_u has dimension %d\n", size(nb, 2))
d = vec(nb) / vec(nb)[2]
@printf("free direction (normalised on P_B):  ΔP_Subs = %+.3f, ΔP_B = %+.3f, Δs = %+.3f\n",
        d[1], d[2], d[3])
println("  i.e. the ONLY freedom left at a stage is to move power between the")
println("  substation and the battery, with the energy slack following.")
@printf("residual ||c_u d|| = %.2e\n", norm(cu * d))

println("\nWith f_uu = c_uu = l_ux = 0 the stage Hessian is DIAGONAL:")
println("  H = diag(0, 2 C_B Δt, 0) + diag(Sig_L + Sig_U) + f_u' V_BB f_u + reg*I")
println("so the reduced Hessian along the single free direction is the SCALAR")
println("  d'Hd = (Sig_L+Sig_U)_Psub + 2 C_B Δt + Δt^2 V_BB + (Sig_L+Sig_U)_P_B")
println("         + Δt^2 (Sig_L+Sig_U)_s ,   Sig_L = z_l/(u-l),  Sig_U = z_u/(u_hi-u)")
println("\nThe 'dense n_u x n_u KKT solve' is therefore a 1x1 division here, and")
println("the Cholesky that can fail is a positivity test on that one number.")
println("\nWith nx = 1 the value derivatives V_B and V_BB are scalars too, so the")
println("whole backward recursion is scalar arithmetic at every stage.")
