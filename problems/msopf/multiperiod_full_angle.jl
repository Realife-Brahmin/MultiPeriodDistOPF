# multiperiod_full_angle.jl -- the full-angle (rectangular AC) model of full_angle_pf.jl over
# a horizon of T periods, coupled in time by each battery's state of charge.
#
# WHY
# ---
# solve_full_angle is a snapshot: a battery there can charge or discharge at its rating with
# no memory, so "how far can batteries stretch the angle window" comes out optimistic. Here
# every battery carries its energy B from period to period,
#     B[t] = B[t-1] - P_B[t] * dt,      SOC_min * B_R <= B[t] <= SOC_max * B_R,
# so a battery that keeps charging to hold back-feed off eventually fills up.
#
# MODEL
# -----
# Per period, exactly the equations of solve_full_angle (Section II of the MSOPF paper):
# source KVL behind Zs, line KVL with explicit series currents, pi-model half-shunts, nodal
# balance, P_Subs >= 0, and the voltage band on buses not set by a substation. Loads scale
# with the load factor lambda_L[t] (P and Q alike). Each PV unit's real power is fixed,
# p_D[t] = lambda_PV[t] * Pmpp; its reactive power q_D[t] is a decision within the inverter
# rating, |q_D[t]| <= sqrt(S_D^2 - p_D[t]^2) (injecting positive) -- as in the user's DDP and
# tADMM models. `pv_q = :unity` pins q_D = 0 instead. Batteries are one real-power variable
# each (discharging positive), lossless -- the user's convention.
#
# Objective, in dollars, with the conventions of problems/mpopf/spar_ddp/examples/power_system:
#     sum_t c[t] * dt * sum_s P_Subs[s,t]                     energy bought (equal prices)
#   + C_B * dt * sum_{k,t} P_B[k,t]^2                         battery cycling cost, P_B in kW
#   + gamma * sum_k (B[k,T] - B[k,0])^2                       only with terminal = :soft
# Terminal state of charge: `terminal = :cyclic` (the default) is the HARD B[T] == B[0], a
# day that can repeat -- the user's choice for MSOPF, where feasibility is the question.
# `:soft` replaces it with the gamma penalty above, the DDP standing rule; here it lets the
# batteries end the day fuller, so its feasible window is one day's borrowing and grows
# without bound with battery size. `:free` has no terminal term at all.
#
# PROFILES
# --------
# load_profile / price_profile / pv_profile mirror problems/mpopf/tadmm/root_level/config.jl with
# PROFILE_PERIODIC=1, the sampling every current DDP instance uses: phase_k = 2 pi (k-1)/T.
# They resample with T; never slice one.

include(joinpath(@__DIR__, "full_angle_pf.jl"))

# Where a capped solve (TIME_LIMIT / ITERATION_LIMIT) is parked: (; file, stage, case), set by a study.
const PARK = Ref{Any}(nothing)

profile_phase(T) = [2pi * (k - 1) / T for k in 1:T]
"Load factor lambda_L[t] in [0.8, 1.0] (tADMM's LoadShapeLoad, periodic sampling)."
load_profile(T) = T == 1 ? [0.97174] : 0.8 .+ 0.2 .* (sin.(profile_phase(T) .- 0.8) .+ 1) ./ 2
"Energy price c[t] in \$/kWh, in [0.08, 0.20] (tADMM's LoadShapeCost, periodic sampling)."
price_profile(T) = T == 1 ? [0.14] : 0.08 .+ 0.12 .* (sin.(profile_phase(T)) .+ 1) ./ 2
"PV factor lambda_PV[t] in [0, 1]: a half-sine over the middle 50% of the horizon (tADMM's LoadShapePV)."
function pv_profile(T)
    pv = zeros(T)
    if T >= 4
        a, b = max(1, round(Int, 0.25T)), min(T, round(Int, 0.75T))
        n = b - a + 1
        pv[a:b] = n >= 2 ? sin.(range(0, pi, length = n)) : [1.0]
    elseif T == 3
        pv[2] = 1.0
    elseif T == 2
        pv .= 0.5
    else
        pv .= 1.0
    end
    return pv
end

"""
    read_pv(system; kwargs...)

The deck's PV units as (name, bus, kW = Pmpp). Pmpp is the rated maximum power; not
`PVsystems.kW()`, which is the unit's present output and reads 0 under a daily shape (the
bug fixed in problems/mpopf/tadmm/parse_opendss.jl on 2026-10-05). Compiles the deck with the PV
systems enabled; compile it again afterwards for read_network.
"""
function read_pv(system; kwargs...)
    compile_deck(system; kwargs...)
    pvs = []
    for name in ODD.PVsystems.AllNames()
        lowercase(name) == "none" && continue
        ODD.Circuit.SetActiveElement("PVSystem." * name)
        ODD.CktElement.Enabled() || continue
        kW = parse(Float64, ODD.Text.Command("? PVSystem.$name.Pmpp"))
        kVA = parse(Float64, ODD.Text.Command("? PVSystem.$name.kVA"))
        kW > 0 || error("PV unit $name has no rated power (Pmpp = $kW)")
        kVA >= kW || error("PV unit $name: inverter rating $kVA kVA is below Pmpp = $kW kW")
        push!(pvs, (; name = lowercase(name), bus = busname(ODD.CktElement.BusNames()[1]), kW, kVA))
    end
    return pvs
end

"Cycling cost by the user's rule (terminal_soc_penalty.jl): the price-driven swing is alpha x the median rating."
cb_by_rule(net; alpha = 3.0, c_swing = 0.06) = c_swing / (2 * alpha * median_of([bt.kW_rated for bt in net.batteries]))
"Terminal-SOC weight by the user's rule: 2 gamma dB = cbar * kVA_B at dB = 10% of the median energy rating."
gamma_by_rule(net; cbar = 0.14) = 5 * cbar * net.S_base_kVA / median_of([bt.kWh_rated / net.S_base_kVA for bt in net.batteries])
median_of(v) = (s = sort(v); n = length(s); isodd(n) ? s[(n+1)÷2] : (s[n÷2] + s[n÷2+1]) / 2)

"δ for source `name` in period t: a number (constant over the horizon) or a vector of length T."
angle_at(delta_deg, name, t) = (d = delta_deg[name]; d isa AbstractVector ? d[t] : d)

"""
    solve_multiperiod(net, delta_deg; T, load_factor, price, dt_h = 24/T, pvs = [], pv_factor = zeros(T),
                      pv_q = :optimize, battery = :optimize, soc = (0.30, 0.95), soc0 = 0.625,
                      terminal = :cyclic, C_B = 0.0, gamma = 0.0, no_backflow = true,
                      V_band = (0.95, 1.05), ideal_sources = false, start = nothing, max_iter = 3000)

The T-period problem described at the top of this file. `delta_deg[name]` is each source's
angle in degrees, a number or a length-T vector. `battery = :idle` fixes every P_B at 0, which
decouples the periods. `pvs` is the list from `read_pv`; `pv_q = :unity` pins every PV's
reactive power at 0. `start` is a previous result of this function, used as the starting
point (a warm start along an angle sweep). With `T = 1` and `soc = (-Inf, Inf)` this is the
single-period battery OPF with unlimited energy.

`objective = :violation` is the feasibility test. The no-backflow and voltage-band limits get
nonnegative slacks and the objective is their sum (plus a 1e-6 P_B^2 term to pick one
dispatch), so the problem always has a solution and `violation` -- total backflow in per unit
plus total band violation in per unit squared -- is zero exactly when the angles are
feasible. Asking Ipopt to certify infeasibility of the hard-constrained problem instead costs
500+ iterations per point (about 30 s at T = 24) and rests on a local solver's failure status.

`objective = :ray` finds how far the angles can move, in one solve. The angles become
delta_s = delta_deg[s] + r * ray_dir[s] with r >= 0 a decision variable (degrees), held
constant over the horizon, and the objective is to maximize r subject to every hard limit:
the largest r for which a feasible day exists. Sources with ray_dir = 0 stay fixed. The
result carries `r`; `r == ray_limit` means the search limit was reached. A local solver finds
the edge of the feasible set it starts in, so start it from a feasible point (`start`) and
check the edge with `objective = :violation` just inside and just outside.

`objective = :margin` measures how much room no-backflow has at the given angles: it maximizes
m subject to every substation importing at least m in every period (`margin_kW`; batteries
and PV reactive power are free to raise it). The margin is zero at a window edge set by
no-backflow and close to linear in the angle there (its slope is the angle-to-import
sensitivity), so a secant iteration on it finds the edge in a few solves. Band violations are
priced at 100 per unit, not forbidden; check `max_violation`.

`tol` and `constr_viol_tol` are Ipopt's convergence and constraint-violation tolerances (both
1e-9 by default). With a looser `tol` the slacks of the feasibility test end near `tol`, not
zero, so their sum grows with the number of limits; judge feasibility by `max_violation`, the
largest single slack, against a threshold above `tol`.

Returns powers in kW/kvar, energies in kWh, voltages in per unit; `cost` splits the objective.
"""
function solve_multiperiod(net, delta_deg; T, load_factor, price, dt_h = 24 / T, pvs = [], pv_factor = zeros(T),
                           pv_q = :optimize, battery = :optimize, soc = (0.30, 0.95), soc0 = 0.625,
                           terminal = :cyclic, C_B = 0.0, gamma = 0.0, no_backflow = true,
                           V_band = (0.95, 1.05), ideal_sources = false, start = nothing, max_iter = 3000,
                           objective = :cost, ray_dir = nothing, ray_limit = 20.0, margin_band_kw = nothing, margin_reg = 1e-7, ampacity_pu = nothing,
                           source_ampacity_pu = nothing, free_angles = nothing, current_target = nothing, margin_band_hard = false,
                           band_buses = nothing, tol = IPOPT_TOL, constr_viol_tol = tol, solver = :ipopt, solver_attributes = ())
    # `solver = :gurobi` solves the same model to global optimality by spatial branch and bound
    # (objectives :cost, :violation and :margin; :ray has cos and sin of a variable).
    # MULTIPERIOD_LEAN=1 hands the model straight to Ipopt (JuMP's direct mode): no cached copy
    # and no bridge layer, which is most of the memory of a 10k-bus day.
    t_start = time()
    lean = solver === :ipopt && get(ENV, "MULTIPERIOD_LEAN", "0") == "1"
    model = lean ? direct_model(Ipopt.Optimizer()) : new_model(solver)
    lean && set_silent(model)
    if solver === :ipopt
        set_attribute(model, "max_iter", max_iter)
        set_attribute(model, "tol", tol)
        set_attribute(model, "constr_viol_tol", constr_viol_tol)
        # Further Ipopt options as "key=value;key=value", the convention of the DDP scripts: e.g.
        # IPOPT_EXTRA_OPTIONS="linear_solver=ma57;hsllib=<dir>/libma57.dll;linear_system_scaling=none"
        # A warm start (`start`, a nearby solution) is only worth something if Ipopt does not push
        # it back to the interior: start the barrier small and keep the point near its bounds.
        # MULTIPERIOD_WARM=0 restores Ipopt's defaults (mu_init 0.1, bound_push 0.01).
        # A wall-clock cap on this solve (seconds; set by a study from its own reference solve).
        # A capped solve is recorded in PARK[] (if set) for a later uncapped replay.
        cap = tryparse(Float64, get(ENV, "MULTIPERIOD_MAX_WALL", ""))
        cap === nothing || set_attribute(model, "max_wall_time", cap)
        if start !== nothing && get(ENV, "MULTIPERIOD_WARM", "1") == "1"
            set_attribute(model, "mu_init", 1e-4)
            set_attribute(model, "bound_push", 1e-6)
            set_attribute(model, "bound_frac", 1e-6)
        end
        for kv in split(get(ENV, "IPOPT_EXTRA_OPTIONS", ""), ';'; keepempty = false)
            k, v = strip.(split(kv, '='; limit = 2))
            set_attribute(model, String(k), something(tryparse(Int, v), tryparse(Float64, v), String(v)))
        end
    end
    for (k, v) in solver_attributes
        set_attribute(model, k, v)
    end
    S = net.S_base_kVA
    buses, lines = net.buses, [l.name for l in net.lines]
    names, bnames = [s.name for s in net.sources], [bt.name for bt in net.batteries]
    pnames = [p.name for p in pvs]
    Ts = 1:T
    E = Dict((s.name, t) => s.pu * cis(deg2rad(angle_at(delta_deg, s.name, t))) for s in net.sources, t in Ts)

    @variable(model, -1.5 <= e[buses, Ts] <= 1.5)
    @variable(model, -1.5 <= f[buses, Ts] <= 1.5)
    @variable(model, -1e5 <= ir[names, Ts] <= 1e5)
    @variable(model, -1e5 <= ii[names, Ts] <= 1e5)
    @variable(model, -1e5 <= lr[lines, Ts] <= 1e5)
    @variable(model, -1e5 <= li[lines, Ts] <= 1e5)
    # Line ampacity (per unit of the deck's current base, one number for every line or a Dict by line
    # name): |I_l|^2 = lr^2 + li^2 <= I_max^2 in every period. Off by default (`nothing`), which keeps
    # every result computed before 2026-10-21 reproducible.
    amp(n) = ampacity_pu isa AbstractDict ? get(ampacity_pu, n, nothing) : ampacity_pu
    for n in lines, t in Ts
        a = amp(n); a === nothing && continue
        set_lower_bound(lr[n, t], -a); set_upper_bound(lr[n, t], a); set_lower_bound(li[n, t], -a); set_upper_bound(li[n, t], a)
        @constraint(model, lr[n, t]^2 + li[n, t]^2 <= a^2)
    end
    # Source current limit (per unit; a number for every source or a Dict by name): |I_s|^2 <= I_max^2.
    samp(n) = source_ampacity_pu isa AbstractDict ? get(source_ampacity_pu, n, nothing) : source_ampacity_pu
    for n in names, t in Ts
        a = samp(n); a === nothing && continue
        @constraint(model, ir[n, t]^2 + ii[n, t]^2 <= a^2)
    end
    @variable(model, P_subs[names, Ts])
    @variable(model, Q_subs[names, Ts])
    @variable(model, P_B[bnames, Ts])                 # per unit, discharging > 0
    @variable(model, B[bnames, 0:T])                  # per unit hours
    @variable(model, Q_PV[pnames, Ts])                # per unit, injecting > 0

    for p in pvs, t in Ts                             # |q_D| <= sqrt(S_D^2 - p_D^2), p_D fixed
        qmax = pv_q === :optimize ? sqrt(max(p.kVA^2 - (pv_factor[t] * p.kW)^2, 0.0)) / S : 0.0
        if qmax > 0
            set_lower_bound(Q_PV[p.name, t], -qmax); set_upper_bound(Q_PV[p.name, t], qmax)
            set_start_value(Q_PV[p.name, t], start === nothing || !haskey(start.Q_PV_kvar, p.name) ? 0.0 : start.Q_PV_kvar[p.name][t] / S)
        else
            fix(Q_PV[p.name, t], 0.0; force = true)
        end
    end

    for bt in net.batteries
        k, BR = bt.name, bt.kWh_rated / S
        fix(B[k, 0], soc0 * BR; force = true)
        for t in Ts
            isfinite(soc[1]) && set_lower_bound(B[k, t], soc[1] * BR)
            isfinite(soc[2]) && set_upper_bound(B[k, t], soc[2] * BR)
            if battery === :optimize
                set_lower_bound(P_B[k, t], -bt.kW_rated / S)
                set_upper_bound(P_B[k, t], bt.kW_rated / S)
            else
                fix(P_B[k, t], 0.0; force = true)
            end
            @constraint(model, B[k, t] == B[k, t-1] - P_B[k, t] * dt_h)
        end
        terminal === :cyclic && battery === :optimize && @constraint(model, B[k, T] == B[k, 0])
    end

    # Starting point: a previous solution, or the flat start rotated to the sources.
    for t in Ts
        V0 = sum(E[(n, t)] for n in names) / length(names)
        for b in buses
            set_start_value(e[b, t], start === nothing ? real(V0) : real(start.V_pu[b][t]))
            set_start_value(f[b, t], start === nothing ? imag(V0) : imag(start.V_pu[b][t]))
        end
        start === nothing && continue
        for l in lines
            set_start_value(lr[l, t], start.raw.lr[l][t]); set_start_value(li[l, t], start.raw.li[l][t])
        end
        for n in names
            set_start_value(ir[n, t], start.raw.ir[n][t]); set_start_value(ii[n, t], start.raw.ii[n][t])
            set_start_value(P_subs[n, t], start.P_subs_kW[n][t] / S); set_start_value(Q_subs[n, t], start.Q_subs_kvar[n][t] / S)
        end
        if battery === :optimize
            for k in bnames
                set_start_value(P_B[k, t], start.P_B_kW[k][t] / S); set_start_value(B[k, t], start.B_kWh[k][t+1] / S)
            end
        end
    end

    # Fixed injections per bus and period: load (P + jQ) less PV (P only), per unit.
    base_load = Dict(b => 0.0im for b in buses)
    for d in net.loads
        base_load[d.bus] += complex(d.kW, d.kvar) / S
    end
    pv_at = Dict(b => 0.0 for b in buses)
    for p in pvs
        pv_at[p.bus] += p.kW / S
    end
    src_at = Dict(s.bus => s.name for s in net.sources)
    batt_at = Dict(b => [bt.name for bt in net.batteries if bt.bus == b] for b in buses)
    pvs_at = Dict(b => [p.name for p in pvs if p.bus == b] for b in buses)

    # Slacks of the feasibility test: backflow at each source, band violation at each bus.
    elastic = objective === :violation
    # The import margin: the largest m with every substation importing at least m in every
    # period (per unit; negative when some backflow is unavoidable). The band keeps slacks,
    # heavily priced, so the problem has a solution at any angle.
    marg = objective === :margin
    # With `margin_band_kw` (kW of margin per pu of squared voltage) the band is part of the margin
    # instead of an elastic constraint: |V|^2 within [lo^2 + m/k, hi^2 - m/k]. The margin then reaches
    # zero at whichever limit binds first, no-backflow or voltage, and the problem has no penalized
    # slacks (past a voltage edge those made Ipopt take hundreds of iterations).
    vmarg = marg && margin_band_kw !== nothing
    # margin_band_hard = true: the band stays a hard constraint and the margin is on no-backflow alone
    # (the restoration problem: max m = -rho with P_Subs >= m); no slacks, no band room in the margin.
    eband = elastic || (marg && !vmarg && !margin_band_hard)
    # band_buses: where the voltage band applies. nothing = every bus not set by a source (as before); a collection =
    # only those buses (user, 2026-10-24: on decks with secondaries, only the customer buses, the secondary side of each
    # service transformer; the primary then has no band).
    bandbuses = V_band === nothing ? String[] : [b for b in buses if !haskey(src_at, b) && (band_buses === nothing || b in band_buses)]
    bandset = Set(bandbuses)
    @variable(model, s_back[elastic && no_backflow ? names : String[], Ts] >= 0, start = 0.0)
    @variable(model, s_vlo[eband ? bandbuses : String[], Ts] >= 0, start = 0.0)
    @variable(model, s_vhi[eband ? bandbuses : String[], Ts] >= 0, start = 0.0)
    m_imp = marg ? @variable(model, base_name = "m_imp", start = 0.0) : 0.0

    # The ray: one scalar r moves the angles along ray_dir, constant over the horizon.
    objective === :ray && ray_dir === nothing && error("objective = :ray needs ray_dir")
    moving(s) = objective === :ray && get(ray_dir, s.name, 0.0) != 0
    @variable(model, 0 <= r_ray <= ray_limit, start = 0.0)
    objective === :ray || fix(r_ray, 0.0; force = true)
    # `free_angles = Dict(name => (lo, hi))` (degrees): those sources' angles become decision variables in
    # [lo, hi], constant over the horizon, started from delta_deg (used by objective = :maxcurrent).
    freeang = free_angles === nothing ? Dict{String,Tuple{Float64,Float64}}() : free_angles
    @variable(model, theta[collect(keys(freeang))])
    for (n, (lo, hi)) in freeang
        set_lower_bound(theta[n], lo); set_upper_bound(theta[n], hi); set_start_value(theta[n], clamp(angle_at(delta_deg, n, 1), lo, hi))
    end

    for t in Ts
        for s in net.sources
            R, X = ideal_sources ? (0.0, 0.0) : reim(s.z_series * net.Y_base)
            Et = E[(s.name, t)]
            if haskey(freeang, s.name)                # a free angle: E = pu * cis(theta), theta a variable
                th = deg2rad(1.0) * theta[s.name]
                @constraint(model, e[s.bus, t] == s.pu * cos(th) - (R * ir[s.name, t] - X * ii[s.name, t]))
                @constraint(model, f[s.bus, t] == s.pu * sin(th) - (R * ii[s.name, t] + X * ir[s.name, t]))
            elseif moving(s)                          # EMF angle depends on r: E = pu * cis(theta)
                th = deg2rad(1.0) * (angle_at(delta_deg, s.name, t) + ray_dir[s.name] * r_ray)
                @constraint(model, e[s.bus, t] == s.pu * cos(th) - (R * ir[s.name, t] - X * ii[s.name, t]))
                @constraint(model, f[s.bus, t] == s.pu * sin(th) - (R * ii[s.name, t] + X * ir[s.name, t]))
            elseif R == 0 && X == 0
                fix(e[s.bus, t], real(Et); force = true); fix(f[s.bus, t], imag(Et); force = true)
            else
                @constraint(model, e[s.bus, t] == real(Et) - (R * ir[s.name, t] - X * ii[s.name, t]))
                @constraint(model, f[s.bus, t] == imag(Et) - (R * ii[s.name, t] + X * ir[s.name, t]))
            end
            @constraint(model, P_subs[s.name, t] == e[s.bus, t] * ir[s.name, t] + f[s.bus, t] * ii[s.name, t])
            @constraint(model, Q_subs[s.name, t] == f[s.bus, t] * ir[s.name, t] - e[s.bus, t] * ii[s.name, t])
        end
        out_P = Dict(b => zero(QuadExpr) for b in buses)
        out_Q = Dict(b => zero(QuadExpr) for b in buses)
        for l in net.lines
            R, X = reim(-net.Y_base / l.Y[1, 2])
            n = l.name
            @constraint(model, e[l.from, t] - e[l.to, t] == R * lr[n, t] - X * li[n, t])
            @constraint(model, f[l.from, t] - f[l.to, t] == R * li[n, t] + X * lr[n, t])
            Y = l.Y ./ net.Y_base
            for (b, sigma, ysh) in ((l.from, 1, Y[1, 1] + Y[1, 2]), (l.to, -1, Y[2, 2] + Y[2, 1]))
                g, bsh = reim(ysh)
                add_to_expression!(out_P[b], sigma * (e[b, t] * lr[n, t] + f[b, t] * li[n, t]) + g * (e[b, t]^2 + f[b, t]^2))
                add_to_expression!(out_Q[b], sigma * (f[b, t] * lr[n, t] - e[b, t] * li[n, t]) - bsh * (e[b, t]^2 + f[b, t]^2))
            end
        end
        for b in buses
            P_in = (haskey(src_at, b) ? P_subs[src_at[b], t] : 0.0) + sum(P_B[k, t] for k in batt_at[b]; init = 0.0)
            Q_in = (haskey(src_at, b) ? Q_subs[src_at[b], t] : 0.0) + sum(Q_PV[k, t] for k in pvs_at[b]; init = 0.0)
            P_fix = load_factor[t] * real(base_load[b]) - pv_factor[t] * pv_at[b]
            @constraint(model, out_P[b] + P_fix == P_in)
            @constraint(model, out_Q[b] + load_factor[t] * imag(base_load[b]) == Q_in)
            if b in bandset
                vroom = vmarg ? (S / margin_band_kw) * m_imp : 0.0
                @constraint(model, e[b, t]^2 + f[b, t]^2 + (eband ? s_vlo[b, t] : 0.0) >= V_band[1]^2 + vroom)
                @constraint(model, e[b, t]^2 + f[b, t]^2 - (eband ? s_vhi[b, t] : 0.0) <= V_band[2]^2 - vroom)
            end
        end
        no_backflow && @constraint(model, [n in names], P_subs[n, t] + (elastic ? s_back[n, t] : 0.0) >= m_imp)
    end

    # Sums are built in place. `sum(container)` adds term by term and copies the growing
    # expression each time: with half a million slacks (10k buses, T = 24) that copying was
    # 88% of a solve's wall time (profiled 2026-10-13: 171 s per solve at T = 6, Ipopt 13 s).
    energy = sum(price[t] * dt_h * S * P_subs[n, t] for n in names, t in Ts)
    pb2 = zero(QuadExpr)
    for k in bnames, t in Ts
        add_to_expression!(pb2, 1.0, P_B[k, t], P_B[k, t])
    end
    cycling = C_B * dt_h * S^2 * pb2
    term = terminal === :soft ? gamma * sum(((B[k, T] - B[k, 0])^2 for k in bnames); init = zero(QuadExpr)) : zero(QuadExpr)
    viol = zero(AffExpr)
    for c in (s_back, s_vlo, s_vhi), v in c
        add_to_expression!(viol, 1.0, v)
    end
    if elastic
        @objective(model, Min, viol + 1e-6 * pb2)
    elseif objective === :maxcurrent
        # current_target = (:source, name) or (:line, name): maximize |I|^2 summed over the periods (T = 1:
        # that hour's current). All the model's constraints stay hard.
        kind, nm = current_target
        @objective(model, Max, kind === :source ? sum(ir[nm, t]^2 + ii[nm, t]^2 for t in Ts) : sum(lr[nm, t]^2 + li[nm, t]^2 for t in Ts))
    elseif objective === :ray
        @objective(model, Max, r_ray - 1e-7 * pb2)
    elseif marg
        @objective(model, Max, m_imp - 100 * viol - margin_reg * pb2)   # margin_reg = 0: the margin alone
    else
        @objective(model, Min, energy + cycling + term)
    end

    t_built = time()
    optimize!(model)
    status = termination_status(model)
    ok = status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) && primal_status(model) == MOI.FEASIBLE_POINT
    # One line per solve (MULTIPERIOD_SOLVE_LOG=0 silences it): where the time of a study goes.
    if get(ENV, "MULTIPERIOD_SOLVE_LOG", "1") == "1"
        its = try barrier_iterations(model) catch; -1 end
        @printf("      [solve] %-9s T=%-3d %-24s %5d it  build %6.1f s  solve %7.1f s  %s  max|angle| %.4f deg\n", objective, T, status, its,
                t_built - t_start, time() - t_built, start === nothing ? "cold" : "warm", maximum(abs, values(delta_deg)))
        flush(stdout)
    end
    if status in (MOI.TIME_LIMIT, MOI.ITERATION_LIMIT) && PARK[] !== nothing
        p = PARK[]
        line = @sprintf("%s,%s,%s,%d,%s,%s,%d,%.1f,%s", p.stage, p.case, objective, T,
                        join([@sprintf("%s=%.6f", k, v) for (k, v) in sort(collect(delta_deg))], ";"), status,
                        (try barrier_iterations(model) catch; -1 end), time() - t_built, start === nothing ? "cold" : "warm")
        new = !isfile(p.file)
        open(p.file, "a") do io
            new && println(io, "stage,case,objective,T,angles_deg,status,iterations,seconds,start")
            println(io, line)
        end
    end
    # Optional failed-solve summaries for checkpointed campaigns. Never extract or reuse
    # an unaccepted primal snapshot; successful solves and the default API are unchanged.
    if get(ENV, "MULTIPERIOD_FAILED_SNAPSHOT", "1") == "0" &&
       !(status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_LOCALLY_SOLVED))
        failed = (; status, ok=false, T, dt_h, time=solve_time(model),
            iterations=(try barrier_iterations(model) catch; missing end), violation=NaN,
            r=NaN, angles_deg=Dict{String,Float64}(), backflow_kW=NaN, margin_kW=NaN,
            max_violation=NaN, cost=(;energy=NaN,cycling=NaN,terminal=NaN,total=NaN),
            P_subs_kW=Dict{String,Vector{Float64}}(),Q_subs_kvar=Dict{String,Vector{Float64}}(),
            P_B_kW=Dict{String,Vector{Float64}}(),Q_PV_kvar=Dict{String,Vector{Float64}}(),
            B_kWh=Dict{String,Vector{Float64}}(),V_pu=Dict{String,Vector{ComplexF64}}(),
            raw=(;lr=Dict(),li=Dict(),ir=Dict(),ii=Dict()))
        release_solver!(model)
        return failed
    end
    val(x) = has_values(model) ? value(x) : NaN
    series(x, keys, scale) = Dict(k => [val(x[k, t]) * scale for t in Ts] for k in keys)
    result = (; status, ok, T, dt_h, time = solve_time(model), iterations = (try barrier_iterations(model) catch; missing end),
            violation = eband ? val(viol) : 0.0, r = objective === :ray ? val(r_ray) : 0.0,
            angles_deg = Dict(s.name => (haskey(freeang, s.name) ? val(theta[s.name]) : angle_at(delta_deg, s.name, 1)) for s in net.sources),
            backflow_kW = elastic ? val(sum(s_back; init = zero(AffExpr))) * S : 0.0,
            margin_kW = marg ? val(m_imp) * S : NaN,
            # largest single limit violation: backflow in per unit, band in per unit squared
            max_violation = eband && has_values(model) ? maximum(maximum(value.(c); init = 0.0) for c in (s_back, s_vlo, s_vhi)) : (eband ? NaN : 0.0),
            cost = (; energy = val(energy), cycling = val(cycling), terminal = val(term), total = val(energy) + val(cycling) + val(term)),
            P_subs_kW = series(P_subs, names, S), Q_subs_kvar = series(Q_subs, names, S),
            P_B_kW = series(P_B, bnames, S), Q_PV_kvar = series(Q_PV, pnames, S),
            B_kWh = Dict(k => [val(B[k, t]) * S for t in 0:T] for k in bnames),
            V_pu = Dict(b => [complex(val(e[b, t]), val(f[b, t])) for t in Ts] for b in buses),
            raw = (; lr = series(lr, lines, 1.0), li = series(li, lines, 1.0), ir = series(ir, names, 1.0), ii = series(ii, names, 1.0)))
    release_solver!(model)
    return result
end

"""
    release_solver!(model)

Free the Ipopt problem behind `model` now. Its memory is allocated by the solver's C library,
outside Julia's heap, so the garbage collector does not feel it and frees it late: a sweep of
a few hundred T-period solves otherwise grows to several GB and the machine starts swapping
(measured 2026-10-11: 3.8 GB after ~50 s of solves, the process then stuck at 5% CPU). Call
only after every value has been read. `MULTIPERIOD_KEEP_SOLVER=1` skips it, for comparison.
"""
function release_solver!(model)
    get(ENV, "MULTIPERIOD_KEEP_SOLVER", "0") == "1" && return
    try
        inner = unsafe_backend(model).inner
        inner === nothing || finalize(inner)
    catch
    end
    # A 10k-bus day leaves gigabytes of model garbage; collect it before the next solve builds
    # its own (the large10kC run reached 11 GB over two dozen solves without this).
    get(ENV, "MULTIPERIOD_GC", "1") == "1" && GC.gc()
    return
end

"`net` at one period: loads scaled by `lambda`, PV added as negative constant-power loads. For solve_full_angle."
function net_at(net, lambda, pvs, pv_lambda)
    loads = Any[merge(d, (; kW = lambda * d.kW, kvar = lambda * d.kvar)) for d in net.loads]
    for p in pvs
        pv_lambda > 0 && push!(loads, (; name = "pv_" * p.name, bus = p.bus, kW = -pv_lambda * p.kW, kvar = 0.0,
                                       kV = net.kV_base, vminpu = 0.0, vmaxpu = 10.0))
    end
    return merge(net, (; loads))
end

"""
    replay_opendss(net, sol, delta_deg; load_factor, pvs = [], pv_factor = zeros(sol.T))

Hand the dispatch of `sol` to OpenDSS period by period -- loads at `Set LoadMult`, each PV
unit as a constant-power generator at its kW and dispatched kvar, each battery at its
dispatched kW -- and return the
largest disagreement with the JuMP solution over the horizon, plus whether OpenDSS converged
and kept every load constant-PQ throughout. The circuit must be the compiled one `net` came
from, in snapshot mode.
"""
function replay_opendss(net, sol, delta_deg; load_factor, pvs = [], pv_factor = zeros(sol.T))
    existing = Set(lowercase.(ODD.Circuit.AllElementNames()))
    for p in pvs                                          # created once per compiled circuit
        "generator.pvinj_" * p.name in existing && continue
        ODD.Text.Command(@sprintf("New Generator.pvinj_%s phases=1 bus1=%s kV=%.6g kW=1e-6 kvar=0 model=1 Vminpu=0.5 Vmaxpu=1.5", p.name, p.bus, net.kV_base))
    end
    dP = dQ = dV = dB = 0.0
    converged = pq = true
    for t in 1:sol.T
        ODD.Text.Command(@sprintf("Set LoadMult=%.12g", load_factor[t]))
        for p in pvs
            kW, kvar = pv_factor[t] * p.kW, sol.Q_PV_kvar[p.name][t]
            on = kW > 0 || abs(kvar) > 1e-9
            ODD.Text.Command(@sprintf("Edit Generator.pvinj_%s enabled=%s kW=%.12g kvar=%.12g", p.name, on ? "true" : "false", kW, kvar))
        end
        d = Dict(s.name => angle_at(delta_deg, s.name, t) for s in net.sources)
        ds = opendss_point(net, d; P_B = Dict(k => sol.P_B_kW[k][t] for k in keys(sol.P_B_kW)))
        converged &= ds.converged
        pq &= all(dd -> dd.vminpu <= abs(ds.V_pu[dd.bus]) * net.kV_base / dd.kV <= dd.vmaxpu, net.loads)
        dP = max(dP, maximum(abs(sol.P_subs_kW[k][t] - ds.P_subs_kW[k]) for k in keys(ds.P_subs_kW)))
        dQ = max(dQ, maximum(abs(sol.Q_subs_kvar[k][t] - ds.Q_subs_kvar[k]) for k in keys(ds.Q_subs_kvar)))
        dV = max(dV, maximum(abs(sol.V_pu[b][t] - ds.V_pu[b]) for b in net.buses))
        isempty(ds.P_B_kW) || (dB = max(dB, maximum(abs(sol.P_B_kW[k][t] - ds.P_B_kW[k]) for k in keys(ds.P_B_kW))))
    end
    ODD.Text.Command("Set LoadMult=1")
    for p in pvs
        ODD.Text.Command("Edit Generator.pvinj_$(p.name) enabled=false")
    end
    return (; converged, pq_band = pq, dP, dQ, dV, dB)
end

# ---- feasible windows of a constant angle difference -------------------------------------------------

"""
    window_edge(try_point, sgn; step = 0.1, limit = 20.0, bisections = 12)

Edge of the feasible interval containing 0 in direction `sgn` (+1 / -1). `try_point(dd, prev)`
returns a result with `.ok`, given the last feasible result `prev` (or `nothing`) as a warm
start. Steps outward to the first infeasible point, then bisects that step -- see
source_impedance_study.jl for why the search limit is never tested first. Returns
(edge, at_limit, last feasible result).
"""
function window_edge(try_point, sgn; step = 0.1, limit = 20.0, bisections = 12, first = nothing)
    inside, best = 0.0, first
    while abs(inside) < limit
        trial = sgn * min(abs(inside) + step, limit)
        r = try_point(trial, best)
        if !r.ok
            lo, hi = inside, trial
            for _ in 1:bisections
                mid = (lo + hi) / 2
                rm = try_point(mid, best)
                rm.ok ? ((lo, best) = (mid, rm)) : (hi = mid)
            end
            return (lo, false, best)
        end
        inside, best = trial, r
    end
    return (inside, true, best)
end

"Feasible window around 0: (lo, hi, feasible at 0, capped flags, and the solutions at 0 and both edges)."
function feasible_window(try_point; kwargs...)
    r0 = try_point(0.0, nothing)
    r0.ok || return (; feasible = false, lo = NaN, hi = NaN, lo_capped = false, hi_capped = false, at0 = r0, at_lo = nothing, at_hi = nothing)
    lo, lc, rlo = window_edge(try_point, -1; first = r0, kwargs...)
    hi, hc, rhi = window_edge(try_point, +1; first = r0, kwargs...)
    return (; feasible = true, lo, hi, lo_capped = lc, hi_capped = hc, at0 = r0, at_lo = rlo, at_hi = rhi)
end
