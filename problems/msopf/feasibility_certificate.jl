# feasibility_certificate.jl -- does a prescribed angle vector admit a daily AC dispatch?
#
# The clean maximum-margin problem (the model of multiperiod_full_angle.jl, all its constraints
# kept, objective = the margin m alone):
#     max m   s.t.   P_Subs[s,t] >= m                          every source and period
#                    lo^2 + m/k <= |V_b,t|^2 <= hi^2 - m/k     every banded bus and period
#                    the network equations, PV reactive limits, battery power, SOC dynamics,
#                    SOC bounds, B[0], and the hard cyclic terminal B[T] == B[0]
# with k (kW of margin per pu of |V|^2) as in MULTI_VMARGIN. m >= 0 exactly when the original
# no-backflow and band constraints hold, so the sign of the optimum decides feasibility.
#
#   L, a lower bound: a local solve (Ipopt) whose dispatch is re-validated here against the
#      original constraints (`validate_dispatch`); L is the margin of that validated state.
#   U, an upper bound: a sparse SDP relaxation of the same problem (`relax_margin`), with the
#      bound itself computed from the relaxation's dual multipliers (`dual_bound`): any dual
#      point, projected onto the dual cones, gives a valid bound once the box bounds of the
#      variables are used, so U does not rest on the solver's status or primal objective.
#   feasible   if L >= 0        infeasible   if U < -tol        unresolved   otherwise
#
# The relaxation. With the angles fixed, every source EMF E_s,t is a constant, the source and
# line currents are linear in the bus voltages (I = (E - V_s)/z_s, I_l = (V_i - V_j)/z_l), and
# each period is a QCQP in x_t = (e_t, f_t). Every product e_i e_j, e_i f_j, ... that the
# constraints use is replaced by a variable, and for every line (i, j) and period t the 5 x 5
# moment matrix of (1, e_i, f_i, e_j, f_j) must be PSD. For a radial network (a tree of lines,
# the sources being fixed EMFs, not lines) this pattern is chordal, so the line blocks are
# exactly the dense per-period relaxation; with loops they are a weaker but still valid outer
# relaxation. The periods couple only through the linear battery and SOC constraints, kept as
# they are: no matrix spans two periods. Valid bounds used: |e|, |f| <= 1.5 (the box of the
# nonlinear model), hence |W| <= 2.25 for every product variable.
using LinearAlgebra, Printf
import Clarabel

# Banded buses: every bus but the source terminals, or only `kw.band_buses` when given (decks with service
# transformers band their customer buses only, 2026-10-24). An unbanded bus has no voltage limit in the
# model, but the relaxation needs a positive lower bound on |V| to bound injection currents: it gets
# |V| >= V_FLOOR, added to the model being certified (and checked by validate_dispatch). V_FLOOR is far
# below any voltage the feasible dispatches reach (primary buses sit near 0.95-1.05 pu).
const V_FLOOR = 0.5
const Z_IDEAL = 1e-6           # pu: a source impedance below this is treated as ideal in the relaxation
banded_bus(kw, src_at, b) = !haskey(src_at, b) && (get(kw, :band_buses, nothing) === nothing || b in kw.band_buses)

"""
    validate_dispatch(net, delta_deg, kw, x; k)

Re-check a dispatch (`x`: a result of solve_multiperiod, or one recovered from the relaxation)
against the original constraints, recomputing every flow from the bus voltages. Returns the
largest violation by kind (per unit: power balance, source equations, bounds, SOC; band in pu
of |V|) and the margin of this state in kW (min over sources and periods of P_Subs, and of the
band room times k).
"""
function validate_dispatch(net, delta_deg, kw, x; k)
    T = kw.T; dt = get(kw, :dt_h, 24 / T); S = net.S_base_kVA; Yb = net.Y_base
    lo, hi = get(kw, :V_band, (0.95, 1.05)); lam = kw.load_factor
    pvs = get(kw, :pvs, []); pvf = get(kw, :pv_factor, zeros(T)); ideal = get(kw, :ideal_sources, false)
    soc = get(kw, :soc, (0.30, 0.95)); soc0 = get(kw, :soc0, 0.625); terminal = get(kw, :terminal, :cyclic)
    battery = get(kw, :battery, :optimize); pv_q = get(kw, :pv_q, :optimize)
    src_at = Dict(s.bus => s for s in net.sources)
    load = Dict(b => 0.0im for b in net.buses); for d in net.loads; load[d.bus] += complex(d.kW, d.kvar) / S; end
    pvP = Dict(b => 0.0 for b in net.buses); for p in pvs; pvP[p.bus] += p.kW / S; end
    viol = Dict("ampacity" => 0.0, "kvl" => 0.0, "balance" => 0.0, "source" => 0.0, "pv_q" => 0.0, "battery" => 0.0, "soc" => 0.0, "band" => 0.0, "backflow" => 0.0)
    up!(kind, v) = (viol[kind] = max(viol[kind], v))
    minP, minroom = Inf, Inf
    for t in 1:T
        V(b) = x.V_pu[b][t]
        out = Dict(b => 0.0im for b in net.buses)
        for l in net.lines
            Y = l.Y ./ Yb; z = -1 / Y[1, 2]
            # KCL with the state's own series current when it has one (I = dV/z would multiply rounding
            # by 1/|z|, 1e8 on the coupling lines), and KVL checked separately as a voltage residual
            I = hasproperty(x, :raw) ? complex(x.raw.lr[l.name][t], x.raw.li[l.name][t]) : (V(l.from) - V(l.to)) / z
            up!("kvl", abs(V(l.from) - V(l.to) - z * I))
            amax = ampacity_of(kw, l.name); amax === nothing || up!("ampacity", max(0.0, abs(I) - amax))
            out[l.from] += V(l.from) * conj(I) + conj(Y[1, 1] + Y[1, 2]) * abs2(V(l.from))
            out[l.to] += V(l.to) * conj(-I) + conj(Y[2, 2] + Y[2, 1]) * abs2(V(l.to))
        end
        inj = Dict(b => 0.0im for b in net.buses)
        for s in net.sources
            E = s.pu * cis(deg2rad(angle_at(delta_deg, s.name, t)))
            z = ideal ? 0.0im : s.z_series * Yb
            if abs(z) < 1e-12
                up!("source", abs(V(s.bus) - E)); Ssrc = complex(x.P_subs_kW[s.name][t], x.Q_subs_kvar[s.name][t]) / S
            else
                Ssrc = V(s.bus) * conj((E - V(s.bus)) / z)
            end
            inj[s.bus] += Ssrc
            sa = src_amp(kw, s.name, t)
            sa === nothing || abs(z) < 1e-12 || up!("ampacity", max(0.0, abs((E - V(s.bus)) / z) - sa))
            up!("backflow", max(0.0, -real(Ssrc))); minP = min(minP, real(Ssrc) * S)
        end
        for bt in net.batteries
            pb = x.P_B_kW[bt.name][t] / S
            inj[bt.bus] += pb
            up!("battery", max(0.0, abs(pb) - (battery === :optimize ? bt.kW_rated / S : 0.0)))
        end
        for p in pvs
            q = x.Q_PV_kvar[p.name][t] / S
            inj[p.bus] += im * q
            qmax = pv_q === :optimize ? sqrt(max(p.kVA^2 - (pvf[t] * p.kW)^2, 0.0)) / S : 0.0
            up!("pv_q", max(0.0, abs(q) - qmax))
        end
        for b in net.buses
            mis = out[b] + lam[t] * load[b] - pvf[t] * pvP[b] - inj[b]
            up!("balance", max(abs(real(mis)), abs(imag(mis))))
            haskey(src_at, b) && continue
            vm = abs(V(b))
            banded_bus(kw, src_at, b) || (up!("band", max(0.0, V_FLOOR - vm)); continue)
            up!("band", max(0.0, lo - vm, vm - hi))
            minroom = min(minroom, (abs2(V(b)) - lo^2) * k, (hi^2 - abs2(V(b))) * k)
        end
    end
    for bt in net.batteries
        Bk = x.B_kWh[bt.name] ./ S; BR = bt.kWh_rated / S
        up!("soc", abs(Bk[1] - soc0 * BR))
        for t in 1:T
            up!("soc", abs(Bk[t+1] - (Bk[t] - x.P_B_kW[bt.name][t] / S * dt)))
            up!("soc", max(0.0, soc[1] * BR - Bk[t+1], Bk[t+1] - soc[2] * BR))
        end
        terminal === :cyclic && battery === :optimize && up!("soc", abs(Bk[T+1] - Bk[1]))
    end
    eq = maximum(viol[c] for c in ("ampacity", "kvl", "balance", "source", "pv_q", "battery", "soc"))
    return (; viol, max_violation = maximum(values(viol)), equations = eq, min_import_kW = minP, min_band_room_kW = minroom,
            margin_kW = min(minP, minroom))
end

"""
Source-current bound for source `n` at period `t` from `source_amp_pu`: a number for every source, a Dict by source
of numbers, or a Dict by source of per-period vectors (derived bounds, one per hour); nothing if absent.
"""
function src_amp(kw, n, t)
    a = get(kw, :source_ampacity_pu, nothing)
    a === nothing && return nothing
    v = a isa AbstractDict ? get(a, n, nothing) : a
    return v isa AbstractVector ? v[t] : v
end
src_amp_max(kw, n) = (v = [src_amp(kw, n, t) for t in 1:kw.T]; any(isnothing, v) ? nothing : maximum(v))

"Line ampacity in per unit from the problem keywords (`ampacity_pu`: a number for every line or a Dict), or nothing."
ampacity_of(kw, n) = (a = get(kw, :ampacity_pu, nothing); a isa AbstractDict ? get(a, n, nothing) : a)

"""
    variable_bounds(net, kw; k, m_lo_kW, vbox)

Bounds valid for every point of the original model with margin m >= m_lo: bus voltage
magnitudes (the band, widened by m_lo/k, on banded buses; the model's box |e|, |f| <= vbox at
source terminals) and, because the lines form a tree, each series current: it carries the
current of every injection on either side of it, so |I_l| is at most the smaller of the two
sides' sums of injection-current bounds (load, PV, battery and PV reactive ratings over the
lowest voltage, line half-shunts, and at a source terminal |E| + |V| over |z_s|).
"""
function variable_bounds(net, kw; k, m_lo_kW, vbox, band_in_margin = true)
    T = kw.T; S = net.S_base_kVA; Yb = net.Y_base
    lo, hi = get(kw, :V_band, (0.95, 1.05)); lam = kw.load_factor
    pvs = get(kw, :pvs, []); pvf = get(kw, :pv_factor, zeros(T)); battery = get(kw, :battery, :optimize)
    src_at = Dict(s.bus => s for s in net.sources)
    c = S / k
    lo^2 + c * m_lo_kW / S > 0 || error("m_lo too low for a positive voltage bound")
    # band in the margin: |V|^2 in [lo^2 + c m, hi^2 - c m] with m >= m_lo; band hard (restoration): [lo^2, hi^2]
    vmin = Dict(b => haskey(src_at, b) ? 0.0 : !banded_bus(kw, src_at, b) ? V_FLOOR : (band_in_margin ? sqrt(lo^2 + c * m_lo_kW / S) : lo) for b in net.buses)
    vmax = Dict(b => (haskey(src_at, b) || !banded_bus(kw, src_at, b)) ? sqrt(2) * vbox : (band_in_margin ? sqrt(hi^2 - c * m_lo_kW / S) : hi) for b in net.buses)
    Ib = Dict(b => 0.0 for b in net.buses)
    for d in net.loads; Ib[d.bus] += maximum(lam) * abs(complex(d.kW, d.kvar)) / S / vmin[d.bus]; end
    for p in pvs; Ib[p.bus] += max(p.kW * maximum(pvf), p.kVA) / S / vmin[p.bus]; end
    for bt in net.batteries; battery === :optimize && (Ib[bt.bus] += bt.kW_rated / S / vmin[bt.bus]); end
    for s in net.sources
        Emax = s.pu
        sa = src_amp_max(kw, s.name)
        Ib[s.bus] += min((Emax + vmax[s.bus]) / abs(s.z_series * Yb), sa === nothing ? Inf : sa)
    end
    for l in net.lines
        Y = l.Y ./ Yb
        Ib[l.from] += abs(Y[1, 1] + Y[1, 2]) * vmax[l.from]; Ib[l.to] += abs(Y[2, 2] + Y[2, 1]) * vmax[l.to]
    end
    adj = Dict(b => Tuple{String,String}[] for b in net.buses)
    for l in net.lines; push!(adj[l.from], (l.to, l.name)); push!(adj[l.to], (l.from, l.name)); end
    total = sum(values(Ib)); C = Dict{String,Float64}()
    for l in net.lines                                       # the to-side of l: buses reached from l.to without l
        side = Set([l.to]); st = [l.to]
        while !isempty(st)
            u = pop!(st)
            for (v, n) in adj[u]; (n == l.name || v in side) && continue; push!(side, v); push!(st, v); end
        end
        # a line on a loop (closed ties) does not split the network: only the total bound holds for it
        sto = sum(Ib[b] for b in side); C[l.name] = l.from in side ? total : min(sto, total - sto)
        a = ampacity_of(kw, l.name); a === nothing || (C[l.name] = min(C[l.name], a))
    end
    return (; vmin, vmax, C)
end

"""
    relax_margin(net, delta_deg, kw; k, m_lo_kW = -1000, m_hi_kW = 1e4, cuts = Symbol[])

The sparse SDP relaxation of the clean maximum-margin problem at fixed angles (see the top of
this file), built in JuMP for Clarabel as `min -m`, over the points with m >= m_lo (enough to
decide the sign of the optimum). Line currents are variables, as in the nonlinear model: for
each line and period the 7 x 7 moment matrix of (1, e_i, f_i, e_j, f_j, I_re, I_im) is PSD,
and the line's two KVL equations hold for every row of it (valid: for a real point the matrix
is x x', so a' x = 0 gives M a = 0). Eliminating the currents instead puts 1/z of the
near-zero-impedance coupling lines (1e8 pu) into the constraints, which stalls the solver.
"""
function relax_margin(net, delta_deg, kw; k, m_lo_kW = -1000.0, m_hi_kW = 1e4, cuts = Symbol[], vbox = 1.5, solver_attributes = (), band_in_margin = true, src_scale = false)
    T = kw.T; dt = get(kw, :dt_h, 24 / T); S = net.S_base_kVA; Yb = net.Y_base; Ts = 1:T
    lo, hi = get(kw, :V_band, (0.95, 1.05)); lam = kw.load_factor
    pvs = get(kw, :pvs, []); pvf = get(kw, :pv_factor, zeros(T)); ideal = get(kw, :ideal_sources, false)
    soc = get(kw, :soc, (0.30, 0.95)); soc0 = get(kw, :soc0, 0.625); terminal = get(kw, :terminal, :cyclic)
    battery = get(kw, :battery, :optimize); pv_q = get(kw, :pv_q, :optimize)
    buses = net.buses; lines = [l.name for l in net.lines]
    bnames = [bt.name for bt in net.batteries]; pnames = [p.name for p in pvs]
    src_at = Dict(s.bus => s for s in net.sources)
    vb = variable_bounds(net, kw; k, m_lo_kW, vbox, band_in_margin)
    model = Model(Clarabel.Optimizer); set_silent(model)
    for (a, v) in solver_attributes; set_attribute(model, a, v); end
    # Bounds that only the certificate needs (voltage box, derived current bounds, product bounds) go
    # to `box`, not to the solver: in the solver they put 1e5-scale rows beside unit-scale ones.
    box = Dict{VariableRef,Tuple{Float64,Float64}}()
    vx(b) = min(vb.vmax[b], vbox)                             # bound on |e|, |f|
    @variable(model, e[buses, Ts]); @variable(model, f[buses, Ts])
    @variable(model, Wee[buses, Ts]); @variable(model, Wff[buses, Ts]); @variable(model, Wef[buses, Ts])
    for b in buses, t in Ts
        a = vx(b); w = vb.vmax[b]^2
        box[e[b, t]] = (-a, a); box[f[b, t]] = (-a, a); box[Wee[b, t]] = (0.0, a^2); box[Wff[b, t]] = (0.0, a^2); box[Wef[b, t]] = (-a^2, a^2)
    end
    # per line: V_i x V_j products, the current, and current x voltage / current x current products
    nm = (:Lee, :Lff, :Lef, :Lfe, :lr, :li, :eilr, :eili, :filr, :fili, :ejlr, :ejli, :fjlr, :fjli, :lrlr, :lili, :lrli)
    Lv = Dict(s => @variable(model, [lines, Ts], base_name = string(s)) for s in nm)
    for l in net.lines, t in Ts
        n = l.name; ai, aj, C = vx(l.from), vx(l.to), vb.C[l.name]
        bd = Dict(:Lee => ai * aj, :Lff => ai * aj, :Lef => ai * aj, :Lfe => ai * aj, :lr => C, :li => C,
                  :eilr => ai * C, :eili => ai * C, :filr => ai * C, :fili => ai * C, :ejlr => aj * C, :ejli => aj * C,
                  :fjlr => aj * C, :fjli => aj * C, :lrlr => C^2, :lili => C^2, :lrli => C^2)
        for s in nm
            v = Lv[s][n, t]
            box[v] = (s in (:lrlr, :lili) ? 0.0 : -bd[s], bd[s])
            # :ibounds -- the one strengthening tried: the derived (valid) current and product bounds
            # also go to the solver, capping the fictitious losses the relaxation can create
            :ibounds in cuts && (set_lower_bound(v, box[v][1]); set_upper_bound(v, box[v][2]))
        end
    end
    @variable(model, P_B[bnames, Ts]); @variable(model, B[bnames, 0:T]); @variable(model, Q_PV[pnames, Ts])
    @variable(model, m_lo_kW / S <= m <= m_hi_kW / S)
    for p in pvs, t in Ts
        qmax = pv_q === :optimize ? sqrt(max(p.kVA^2 - (pvf[t] * p.kW)^2, 0.0)) / S : 0.0
        qmax > 0 ? (set_lower_bound(Q_PV[p.name, t], -qmax); set_upper_bound(Q_PV[p.name, t], qmax)) : fix(Q_PV[p.name, t], 0.0; force = true)
    end
    for bt in net.batteries
        kk, BR = bt.name, bt.kWh_rated / S
        fix(B[kk, 0], soc0 * BR; force = true)
        for t in Ts
            set_lower_bound(B[kk, t], soc[1] * BR); set_upper_bound(B[kk, t], soc[2] * BR)
            battery === :optimize ? (set_lower_bound(P_B[kk, t], -bt.kW_rated / S); set_upper_bound(P_B[kk, t], bt.kW_rated / S)) : fix(P_B[kk, t], 0.0; force = true)
            @constraint(model, B[kk, t] == B[kk, t-1] - P_B[kk, t] * dt)
        end
        terminal === :cyclic && battery === :optimize && @constraint(model, B[kk, T] == B[kk, 0])
    end
    load = Dict(b => 0.0im for b in buses); for d in net.loads; load[d.bus] += complex(d.kW, d.kvar) / S; end
    pvP = Dict(b => 0.0 for b in buses); for p in pvs; pvP[p.bus] += p.kW / S; end
    battat = Dict(b => [bt.name for bt in net.batteries if bt.bus == b] for b in buses)
    pvat = Dict(b => [p.name for p in pvs if p.bus == b] for b in buses)
    cmul(a::Complex, re, imx) = (real(a) * re - imag(a) * imx, real(a) * imx + imag(a) * re)    # a * (re + j imx)
    psd = []
    for t in Ts
        VV(b) = Wee[b, t] + Wff[b, t]
        outP = Dict(b => zero(AffExpr) for b in buses); outQ = Dict(b => zero(AffExpr) for b in buses)
        for l in net.lines
            i, j, n = l.from, l.to, l.name
            g(s) = Lv[s][n, t]
            R, X = reim(-Yb / l.Y[1, 2])
            M = [1        e[i, t]     f[i, t]     e[j, t]     f[j, t]     g(:lr)   g(:li);
                 e[i, t]  Wee[i, t]   Wef[i, t]   g(:Lee)     g(:Lef)     g(:eilr) g(:eili);
                 f[i, t]  Wef[i, t]   Wff[i, t]   g(:Lfe)     g(:Lff)     g(:filr) g(:fili);
                 e[j, t]  g(:Lee)     g(:Lfe)     Wee[j, t]   Wef[j, t]   g(:ejlr) g(:ejli);
                 f[j, t]  g(:Lef)     g(:Lff)     Wef[j, t]   Wff[j, t]   g(:fjlr) g(:fjli);
                 g(:lr)   g(:eilr)    g(:filr)    g(:ejlr)    g(:fjlr)    g(:lrlr) g(:lrli);
                 g(:li)   g(:eili)    g(:fili)    g(:ejli)    g(:fjli)    g(:lrli) g(:lili)]
            push!(psd, (; t, line = n, c = @constraint(model, Symmetric(M) in PSDCone())))
            amax = ampacity_of(kw, n)                                       # ampacity, lifted: |I|^2 moment <= I_max^2
            amax === nothing || @constraint(model, g(:lrlr) + g(:lili) <= amax^2)
            # KVL  e_i - e_j - R lr + X li = 0,  f_i - f_j - X lr - R li = 0, on every row of M
            ar = [0, 1, 0, -1, 0, -R, X]; ai_ = [0, 0, 1, 0, -1, -X, -R]
            for row in 1:7
                @constraint(model, sum(M[row, c] * ar[c] for c in 1:7 if ar[c] != 0) == 0)
                @constraint(model, sum(M[row, c] * ai_[c] for c in 1:7 if ai_[c] != 0) == 0)
            end
            add_to_expression!(outP[i], g(:eilr) + g(:fili)); add_to_expression!(outQ[i], g(:filr) - g(:eili))
            add_to_expression!(outP[j], -(g(:ejlr) + g(:fjli))); add_to_expression!(outQ[j], -(g(:fjlr) - g(:ejli)))
            Y = l.Y ./ Yb
            for (b, ysh) in ((i, Y[1, 1] + Y[1, 2]), (j, Y[2, 2] + Y[2, 1]))
                add_to_expression!(outP[b], real(ysh) * VV(b)); add_to_expression!(outQ[b], -imag(ysh) * VV(b))
            end
            if :rlt in cuts                                                  # McCormick on the voltage products
                for (p, a, bb) in ((g(:Lee), e[i, t], e[j, t]), (g(:Lff), f[i, t], f[j, t]), (g(:Lef), e[i, t], f[j, t]), (g(:Lfe), f[i, t], e[j, t]))
                    ua, ub = vx(i), vx(j)
                    @constraint(model, p >= -ub * a - ua * bb - ua * ub); @constraint(model, p >= ub * a + ua * bb - ua * ub)
                    @constraint(model, p <= ub * a - ua * bb + ua * ub);  @constraint(model, p <= -ub * a + ua * bb + ua * ub)
                end
            end
        end
        for b in buses
            Pin, Qin = zero(AffExpr), zero(AffExpr)
            for kk in battat[b]; add_to_expression!(Pin, P_B[kk, t]); end
            for kk in pvat[b]; add_to_expression!(Qin, Q_PV[kk, t]); end
            if haskey(src_at, b) && (ideal || abs(src_at[b].z_series * Yb) < Z_IDEAL)
                # A (near-)ideal source: |z_s| = 1.9e-10 pu on the decks with service transformers, whose
                # 1/z_s = 5e9 coefficients stall the solver. Its terminal voltage is E (to within |z_s| |I_s|,
                # below 1e-9 pu here), so V_s and its products are fixed, and its power is whatever the bus balance
                # needs: free variables Ps, Qs with the no-backflow margin, and its current limit as the exact cone
                # |S_s| <= |E| I_max (|V_s| = |E|). The box keeps them within that limit for dual_bound.
                s = src_at[b]; E = s.pu * cis(deg2rad(angle_at(delta_deg, s.name, t)))
                sa = src_amp(kw, s.name, t)
                sa === nothing && error("near-ideal source $(s.name) needs a source current limit (source_ampacity_pu) to bound its power")
                for (v, val) in ((e[b, t], real(E)), (f[b, t], imag(E)), (Wee[b, t], real(E)^2), (Wff[b, t], imag(E)^2), (Wef[b, t], real(E) * imag(E)))
                    fix(v, val; force = true); box[v] = (val, val)
                end
                Smax = abs(E) * sa
                Ps = @variable(model, base_name = "Ps_ideal[$(s.name),$(t)]"); Qs = @variable(model, base_name = "Qs_ideal[$(s.name),$(t)]")
                box[Ps] = (-Smax, Smax); box[Qs] = (-Smax, Smax)
                add_to_expression!(Pin, Ps); add_to_expression!(Qin, Qs)
                @constraint(model, Ps >= m)                                         # no-backflow, with the margin
                @constraint(model, [Smax, Ps, Qs] in SecondOrderCone())
            elseif haskey(src_at, b)
                s = src_at[b]; E = s.pu * cis(deg2rad(angle_at(delta_deg, s.name, t)))
                z = ideal ? 0.0im : s.z_series * Yb
                Pc = e[b, t] * real(E) + f[b, t] * imag(E) - VV(b); Qc = f[b, t] * real(E) - e[b, t] * imag(E)   # V conj(E) - |V|^2
                Ps, Qs = cmul(1 / conj(z), Pc, Qc)
                add_to_expression!(Pin, Ps); add_to_expression!(Qin, Qs)
                @constraint(model, Ps >= m)                                         # no-backflow, with the margin
                sa = src_amp(kw, s.name, t)
                # |I_s|^2 = |E - V_s|^2 / |z_s|^2 = (|E|^2 - 2 Re(conj(E) V_s) + |V_s|^2) / |z_s|^2, lifted: linear in the moments.
                # src_scale = true writes the same constraint times |z_s|^2 (unit-size coefficients instead of 1/|z_s|^2 = 2.8e4).
                if sa !== nothing
                    src_scale ? @constraint(model, abs2(E) - 2 * (real(E) * e[b, t] + imag(E) * f[b, t]) + VV(b) <= sa^2 * abs2(z)) :
                                @constraint(model, (abs2(E) - 2 * (real(E) * e[b, t] + imag(E) * f[b, t]) + VV(b)) / abs2(z) <= sa^2)
                end
            elseif !banded_bus(kw, src_at, b)
                @constraint(model, VV(b) >= V_FLOOR^2)                              # unbanded: the floor only
            else
                if band_in_margin                                                   # the clean margin problem
                    @constraint(model, VV(b) >= lo^2 + (S / k) * m); @constraint(model, VV(b) <= hi^2 - (S / k) * m)
                else                                                                # restoration: band hard, m = -rho on no-backflow only
                    @constraint(model, VV(b) >= lo^2); @constraint(model, VV(b) <= hi^2)
                end
            end
            @constraint(model, outP[b] + lam[t] * real(load[b]) - pvf[t] * pvP[b] == Pin)
            @constraint(model, outQ[b] + lam[t] * imag(load[b]) == Qin)
        end
    end
    # ---- Bus moment blocks (cuts :kcl0 and :jbound), 2026-10-22 ----------------------------------------
    # For a bus b, x_b = (1, e_b, f_b, I_1r, I_1i, ..., I_dr, I_di) over its incident lines; the block
    # reuses the line blocks' variables (V x I, |I|^2 moments) and adds the products of currents of
    # different lines. J_b = sum_l sigma_l I_l + y_sh,b V_b is the current the bus's devices inject
    # (sigma = +1 if b is the line's from-end); y_sh,b sums the half-shunts at b. a_r, a_i: its real
    # and imaginary parts as coefficient vectors on x_b.
    #  :kcl0   at a zero-injection bus (no load, PV or battery, not a source terminal) with a positive
    #          voltage bound: S = V conj(J) = 0 and V != 0 give J = 0, linear; lifted: M a_r = M a_i = 0.
    #  :jbound at a banded bus with devices: |J| = |S| / |V| <= Smax / Vmin, Smax from the fixed net real
    #          power and the worst end of the allowed PV reactive range; lifted: a_r'M a_r + a_i'M a_i <= Jmax^2.
    # Both hold for every original point with m >= m_lo (for a real point M = x x', so they are the
    # original equations and bounds multiplied out). Source terminals are excluded (no band there).
    xcur = Dict{Tuple{String,String,String,Symbol},VariableRef}()      # (bus, line k, line m, :rr/:ri/:ir/:ii)
    busblocks = []
    if :kcl0 in cuts || :jbound in cuts
        inc = Dict(b => [(l, l.from == b ? 1 : -1) for l in net.lines if l.from == b || l.to == b] for b in buses)
        hasdev(b) = load[b] != 0 || pvP[b] != 0 || !isempty(battat[b]) || !isempty(pvat[b])
        for t in Ts
            g(sy, n) = Lv[sy][n, t]
            for b in buses
                haskey(src_at, b) && continue
                zero_inj = !hasdev(b)
                use = (zero_inj && :kcl0 in cuts) || (!zero_inj && :jbound in cuts)
                (use && vb.vmin[b] > 0) || continue
                ls = inc[b]; d = length(ls); n = 3 + 2d
                X = Matrix{AffExpr}(undef, n, n)
                X[1, 1] = AffExpr(1.0); X[1, 2] = X[2, 1] = 1.0 * e[b, t]; X[1, 3] = X[3, 1] = 1.0 * f[b, t]
                X[2, 2] = 1.0 * Wee[b, t]; X[2, 3] = X[3, 2] = 1.0 * Wef[b, t]; X[3, 3] = 1.0 * Wff[b, t]
                for (q, (l, sg)) in enumerate(ls)
                    r, ii = 2 + 2q, 3 + 2q; nm = l.name; fr = sg == 1
                    X[1, r] = X[r, 1] = 1.0 * g(:lr, nm); X[1, ii] = X[ii, 1] = 1.0 * g(:li, nm)
                    X[2, r] = X[r, 2] = 1.0 * g(fr ? :eilr : :ejlr, nm); X[2, ii] = X[ii, 2] = 1.0 * g(fr ? :eili : :ejli, nm)
                    X[3, r] = X[r, 3] = 1.0 * g(fr ? :filr : :fjlr, nm); X[3, ii] = X[ii, 3] = 1.0 * g(fr ? :fili : :fjli, nm)
                    X[r, r] = 1.0 * g(:lrlr, nm); X[ii, ii] = 1.0 * g(:lili, nm); X[r, ii] = X[ii, r] = 1.0 * g(:lrli, nm)
                end
                for (q1, (l1, _)) in enumerate(ls), (q2, (l2, _)) in enumerate(ls)
                    q2 > q1 || continue
                    c1, c2 = vb.C[l1.name], vb.C[l2.name]
                    for (s, (a1, a2)) in ((:rr, (2 + 2q1, 2 + 2q2)), (:ri, (2 + 2q1, 3 + 2q2)), (:ir, (3 + 2q1, 2 + 2q2)), (:ii, (3 + 2q1, 3 + 2q2)))
                        v = @variable(model, base_name = "xc_$(s)[$(b),$(l1.name),$(l2.name),$(t)]")
                        xcur[(string(b, ",", t), l1.name, l2.name, s)] = v; box[v] = (-c1 * c2, c1 * c2)
                        X[a1, a2] = X[a2, a1] = 1.0 * v
                    end
                end
                push!(psd, (; t, line = "bus:" * b, c = @constraint(model, Symmetric(X) in PSDCone())))
                gsh = sum((sg == 1 ? (l.Y[1, 1] + l.Y[1, 2]) : (l.Y[2, 2] + l.Y[2, 1])) / Yb for (l, sg) in ls)
                ar = zeros(n); ai = zeros(n)
                ar[2], ar[3] = real(gsh), -imag(gsh); ai[2], ai[3] = imag(gsh), real(gsh)   # y_sh V = (g e - b f) + j (g f + b e)
                for (q, (l, sg)) in enumerate(ls); ar[2+2q] = sg; ai[3+2q] = sg; end
                if zero_inj
                    for row in 1:n
                        @constraint(model, sum(ar[c] * X[row, c] for c in 1:n if ar[c] != 0) == 0)
                        @constraint(model, sum(ai[c] * X[row, c] for c in 1:n if ai[c] != 0) == 0)
                    end
                    push!(busblocks, (; t, bus = b, kind = :kcl0, n))
                else
                    Pnet = -lam[t] * real(load[b]) + pvf[t] * pvP[b]
                    Pb = sum((bt.kW_rated / S for bt in net.batteries if bt.bus == b && battery === :optimize); init = 0.0)
                    qd = sum((pv_q === :optimize ? sqrt(max(p.kVA^2 - (pvf[t] * p.kW)^2, 0.0)) / S : 0.0 for p in pvs if p.bus == b); init = 0.0)
                    Qnet = -lam[t] * imag(load[b])
                    Smax = sqrt((abs(Pnet) + Pb)^2 + max(abs(Qnet + qd), abs(Qnet - qd))^2)
                    Jmax = Smax / vb.vmin[b]
                    @constraint(model, sum(ar[c1] * ar[c2] * X[c1, c2] + ai[c1] * ai[c2] * X[c1, c2] for c1 in 1:n, c2 in 1:n if ar[c1] * ar[c2] != 0 || ai[c1] * ai[c2] != 0) <= Jmax^2)
                    push!(busblocks, (; t, bus = b, kind = :jbound, n, Jmax))
                end
            end
        end
    end
    @objective(model, Min, -m)
    return (; model, e, f, Wee, Wff, Wef, P_B, B, Q_PV, m, psd, S, bounds = vb, n = 7, box, xcur, busblocks, Lv)
end

"""
    lifted_violation(R, net, x, kw)

The relaxation's validity check: lift a solution `x` of the nonlinear model (with its line
currents, `x.raw`) to the relaxation's variables (every product from the actual values, the
margin from the state) and return the largest constraint violation and the constraint behind it.
Every point of the original model must pass (up to the solver's tolerance).
"""
function lifted_violation(R, net, x, kw; margin_kW)
    T = kw.T; S = R.S; m = R.model; pt = Dict{VariableRef,Float64}()
    for b in net.buses, t in 1:T
        V = x.V_pu[b][t]; ev, fv = real(V), imag(V)
        pt[R.e[b, t]] = ev; pt[R.f[b, t]] = fv; pt[R.Wee[b, t]] = ev^2; pt[R.Wff[b, t]] = fv^2; pt[R.Wef[b, t]] = ev * fv
    end
    for l in net.lines, t in 1:T
        Vi, Vj = x.V_pu[l.from][t], x.V_pu[l.to][t]; ei, fi, ej, fj = real(Vi), imag(Vi), real(Vj), imag(Vj)
        lr, li = x.raw.lr[l.name][t], x.raw.li[l.name][t]
        val = Dict(:Lee => ei * ej, :Lff => fi * fj, :Lef => ei * fj, :Lfe => fi * ej, :lr => lr, :li => li,
                   :eilr => ei * lr, :eili => ei * li, :filr => fi * lr, :fili => fi * li, :ejlr => ej * lr, :ejli => ej * li,
                   :fjlr => fj * lr, :fjli => fj * li, :lrlr => lr^2, :lili => li^2, :lrli => lr * li)
        for (s, v) in val; pt[variable_by_name(m, string(s) * "[" * l.name * "," * string(t) * "]")] = v; end
    end
    for kk in axes(R.P_B, 1), t in 1:T; pt[R.P_B[kk, t]] = x.P_B_kW[kk][t] / S; end
    for kk in axes(R.B, 1), t in 0:T; pt[R.B[kk, t]] = x.B_kWh[kk][t+1] / S; end
    for kk in axes(R.Q_PV, 1), t in 1:T; pt[R.Q_PV[kk, t]] = x.Q_PV_kvar[kk][t] / S; end
    for ((bt, k1, k2, sy), v) in R.xcur                                        # bus blocks: products of two lines' currents
        bb, tt = split(bt, ','); t = parse(Int, tt)
        I1 = complex(x.raw.lr[k1][t], x.raw.li[k1][t]); I2 = complex(x.raw.lr[k2][t], x.raw.li[k2][t])
        pt[v] = sy === :rr ? real(I1) * real(I2) : sy === :ri ? real(I1) * imag(I2) : sy === :ir ? imag(I1) * real(I2) : imag(I1) * imag(I2)
    end
    pt[R.m] = margin_kW / S
    rep = primal_feasibility_report(m, pt; atol = 0.0)
    # the box the certificate uses must contain the lifted point too (a tightened box excluding a real point is invalid)
    bw, bat = 0.0, ""
    for (v, (l, u)) in R.box
        haskey(pt, v) || continue
        d = max(l - pt[v], pt[v] - u, 0.0); d > bw && (bw = d; bat = name(v))
    end
    isempty(rep) && return (; worst = 0.0, at = "", box_worst = bw, box_at = bat)
    (c, v) = argmax(last, collect(rep))
    return (; worst = v, at = string(c), box_worst = bw, box_at = bat)
end

"""
    start_from_voltages(net, delta_deg, kw, V)

A `start` for solve_multiperiod built from bus voltages alone (`V[b][t]`): currents, source
powers and an idle battery from the network equations. Used to start the local solver on other
branches of the AC equations (for example scrambled bus angles).
"""
function start_from_voltages(net, delta_deg, kw, V)
    T = kw.T; S = net.S_base_kVA; Yb = net.Y_base; Ts = 1:T
    lr = Dict{String,Vector{Float64}}(); li = Dict{String,Vector{Float64}}()
    for l in net.lines
        z = -Yb / l.Y[1, 2]; I = [(V[l.from][t] - V[l.to][t]) / z for t in Ts]
        lr[l.name] = real.(I); li[l.name] = imag.(I)
    end
    ir = Dict{String,Vector{Float64}}(); ii = Dict{String,Vector{Float64}}(); P = Dict{String,Vector{Float64}}(); Q = Dict{String,Vector{Float64}}()
    for s in net.sources
        z = s.z_series * Yb
        I = [(s.pu * cis(deg2rad(angle_at(delta_deg, s.name, t))) - V[s.bus][t]) / z for t in Ts]
        Ssrc = [V[s.bus][t] * conj(I[t]) for t in Ts]
        ir[s.name] = real.(I); ii[s.name] = imag.(I); P[s.name] = real.(Ssrc) .* S; Q[s.name] = imag.(Ssrc) .* S
    end
    BR = Dict(bt.name => bt.kWh_rated for bt in net.batteries)
    return (; V_pu = V, raw = (; lr, li, ir, ii), P_subs_kW = P, Q_subs_kvar = Q,
            P_B_kW = Dict(k => zeros(T) for k in keys(BR)), B_kWh = Dict(k => fill(get(kw, :soc0, 0.625) * BR[k], T + 1) for k in keys(BR)),
            Q_PV_kvar = Dict(p.name => zeros(T) for p in get(kw, :pvs, [])))
end

"Vector index of entry (i, j), i <= j, in MOI's upper-triangle column-wise order."
tri(i, j) = (j * (j - 1)) ÷ 2 + i

"""
    dual_bound(model; c = objective_function(model))

A lower bound on  min c'x  over every point that satisfies the model's constraints and its
variable bounds, built from the solver's dual values: each multiplier is projected onto its
dual cone (>= 0, <= 0, free, or -- for a PSD triangle -- the matrix with halved off-diagonals
clipped to PSD), the Lagrangian c'x - sum lambda' g(x) is then a valid under-estimator, and its
minimum over the variable box is attained at a vertex. Valid for any multipliers (up to
floating-point rounding, which is bounded and reported as `rounding`); its tightness depends on
how good the multipliers are. Also reports how far the raw multipliers were from their cones.
"""
function dual_bound(model; c = objective_function(model), box = Dict{VariableRef,Tuple{Float64,Float64}}(), detail = 0)
    c = convert(AffExpr, c)                                    # a bare variable objective, too
    r = Dict{VariableRef,Float64}(); c0 = c.constant
    for (v, a) in c.terms; r[v] = get(r, v, 0.0) + a; end
    absacc = abs(c0); psd_viol, sign_viol = 0.0, 0.0
    csrc = Dict{String,Float64}()                              # where the constant comes from, by constraint kind
    bk = backend(model)
    for (F, Sset) in list_of_constraint_types(model)
        F == VariableRef && continue                           # bounds: the box, used below
        for cref in all_constraints(model, F, Sset)
            obj = constraint_object(cref)
            if obj.set isa MOI.PositiveSemidefiniteConeTriangle
                lamv = MOI.get(bk, MOI.ConstraintDual(), index(cref))
                fv = moi_function(obj)
                n = obj.set.side_dimension
                # MOI's dual for PositiveSemidefiniteConeTriangle is the upper triangle of the dual matrix Z,
                # paired with the constraint matrix by the Frobenius product: off-diagonal entries count twice
                # (checked on toy problems, certificate_validation.jl part 0). Z must be PSD: clip its eigenvalues.
                L = zeros(n, n)
                for j in 1:n, i in 1:j
                    L[i, j] = L[j, i] = lamv[tri(i, j)]
                end
                ev = eigen(Symmetric(L)); psd_viol = max(psd_viol, -minimum(ev.values))
                Lp = ev.vectors * Diagonal(max.(ev.values, 0.0)) * ev.vectors'
                lp = [(i == j ? 1.0 : 2.0) * Lp[i, j] for j in 1:n for i in 1:j]        # weight 2 off the diagonal
                for term in fv.terms
                    v = VariableRef(model, term.scalar_term.variable)
                    r[v] = get(r, v, 0.0) - lp[term.output_index] * term.scalar_term.coefficient
                    absacc += abs(lp[term.output_index] * term.scalar_term.coefficient)
                end
                for (q, cst) in enumerate(fv.constants)
                    c0 -= lp[q] * cst; absacc += abs(lp[q] * cst); csrc["PSD constant entries"] = get(csrc, "PSD constant entries", 0.0) - lp[q] * cst
                end
            elseif obj.set isa MOI.SecondOrderCone
                # f(x) in SOC with multiplier y in SOC (self-dual): y'f(x) >= 0 on the feasible set, so subtracting it keeps
                # the bound valid; y projected onto the cone: (t, u) -> kept if |u| <= t, 0 if |u| <= -t, else ((t+|u|)/2)(1, u/|u|)
                yv = MOI.get(bk, MOI.ConstraintDual(), index(cref)); fv = moi_function(obj)
                tt, uu = yv[1], yv[2:end]; nu = norm(uu)
                yp = nu <= tt ? copy(yv) : nu <= -tt ? zero(yv) : vcat((tt + nu) / 2, (tt + nu) / 2 .* uu ./ nu)
                sign_viol = max(sign_viol, norm(yv - yp))
                for term in fv.terms
                    v = VariableRef(model, term.scalar_term.variable)
                    r[v] = get(r, v, 0.0) - yp[term.output_index] * term.scalar_term.coefficient
                    absacc += abs(yp[term.output_index] * term.scalar_term.coefficient)
                end
                for (q, cst) in enumerate(fv.constants)
                    c0 -= yp[q] * cst; absacc += abs(yp[q] * cst); csrc["SOC constant entries"] = get(csrc, "SOC constant entries", 0.0) - yp[q] * cst
                end
            else
                lam = dual(cref); fn = obj.func; b = MOI.constant(obj.set)
                lp = obj.set isa MOI.GreaterThan ? max(lam, 0.0) : obj.set isa MOI.LessThan ? min(lam, 0.0) : lam
                sign_viol = max(sign_viol, abs(lam - lp))
                for (v, a) in fn.terms
                    r[v] = get(r, v, 0.0) - lp * a; absacc += abs(lp * a)
                end
                c0 -= lp * (fn.constant - b); absacc += abs(lp * (fn.constant - b))
                kname = string(nameof(typeof(obj.set))); csrc[kname] = get(csrc, kname, 0.0) - lp * (fn.constant - b)
            end
        end
    end
    lb, box_term = c0, 0.0; c0_before_box = c0
    contrib = Tuple{Float64,String,Float64,Float64,Float64}[]   # (|term|, variable, residual coefficient, lower, upper)
    for (v, a) in r
        a == 0 && continue
        l, u = is_fixed(v) ? (fix_value(v), fix_value(v)) : (has_lower_bound(v) ? lower_bound(v) : -Inf, has_upper_bound(v) ? upper_bound(v) : Inf)
        haskey(box, v) && (l = max(l, box[v][1]); u = min(u, box[v][2]))
        val = a > 0 ? a * l : a * u
        isfinite(val) || return (; bound = -Inf, unbounded_variable = name(v), psd_viol, sign_viol, box_term, rounding = 0.0)
        lb += val; box_term += abs(val); absacc += abs(val)
        detail > 0 && push!(contrib, (abs(val), name(v), a, l, u))
    end
    top = detail > 0 ? sort(contrib; rev = true)[1:min(detail, length(contrib))] : contrib
    rounding = 4 * eps() * absacc * 10                         # a generous bound on the summation error
    return (; bound = lb - rounding, psd_viol, sign_viol, box_term, rounding, unbounded_variable = nothing, constant = c0_before_box, top, csrc)
end

"Interval product [a, b] x [c, d] and interval square."
iprod(x, y) = (p = (x[1] * y[1], x[1] * y[2], x[2] * y[1], x[2] * y[2]); (minimum(p), maximum(p)))
isq(x) = x[1] <= 0 <= x[2] ? (0.0, max(x[1]^2, x[2]^2)) : (min(x[1]^2, x[2]^2), max(x[1]^2, x[2]^2))
iand(x, y) = (max(x[1], y[1]), min(x[2], y[2]))

"""
    tighten_boxes!(R, net, kw, delta_deg)

Tighten the certificate's variable box with consequences of the restricted model (source and line
current limits in `kw`; valid for every point of that model, not only for observed dispatches):

* a source with |I_s| <= I_s,max: V_s = E_s - Z_s I_s, so |V_s - E_s| <= |Z_s| I_s,max;
* a line with |I_l| <= I_l,max (the tree bound `R.bounds.C`, or the line limit if smaller):
  V_i - V_j = z_l I_l, so |V_i - V_j| <= |z_l| I_l,max;
* hence, by the triangle inequality along the tree path from a source terminal to bus b,
  |V_b - E_s| <= |Z_s| I_s,max + sum over the path of |z_l| I_l,max  -- a disc around E_s; each bus
  keeps the smallest such disc, and its square box bounds e_b and f_b (intersected with the old box);
* every moment variable stands for a product of two of these quantities, so it lies in the interval
  product of their boxes (squares: the interval square), for the bus blocks, the line blocks and the
  cross-current products alike.
Returns the radii by bus. Only `R.box` changes; the solved model does not.
"""
function tighten_boxes!(R, net, kw, delta_deg)
    Yb = net.Y_base; box = R.box; vb = R.bounds
    Is(n, t) = src_amp(kw, n, t)
    Il(l) = vb.C[l.name]                                         # already min(tree bound, line limit)
    adj = Dict(b => Tuple{String,Float64}[] for b in net.buses)
    for l in net.lines
        w = abs(-Yb / l.Y[1, 2]) * Il(l)
        push!(adj[l.from], (l.to, w)); push!(adj[l.to], (l.from, w))
    end
    T = kw.T; radius = nothing
    for t in 1:T
        radius = Dict(b => Inf for b in net.buses); center = Dict{String,String}()
        for s in net.sources
            a = Is(s.name, t); a === nothing && continue
            d = Dict(s.bus => abs(s.z_series * Yb) * a); st = [s.bus]          # a tree: one path to each bus
            while !isempty(st)
                u = pop!(st)
                for (v, w) in adj[u]
                    haskey(d, v) && continue
                    d[v] = d[u] + w; push!(st, v)
                end
            end
            for (b, r) in d
                r < radius[b] && (radius[b] = r; center[b] = s.name)
            end
        end
        ex = Dict{String,Tuple{Float64,Float64}}(); fx = Dict{String,Tuple{Float64,Float64}}()
        for b in net.buses
            e0, f0 = box[R.e[b, t]], box[R.f[b, t]]
            if isfinite(radius[b])
                s = only(x for x in net.sources if x.name == center[b]); E = s.pu * cis(deg2rad(angle_at(delta_deg, s.name, t)))
                e0 = iand(e0, (real(E) - radius[b], real(E) + radius[b])); f0 = iand(f0, (imag(E) - radius[b], imag(E) + radius[b]))
            end
            ex[b], fx[b] = e0, f0
            box[R.e[b, t]] = e0; box[R.f[b, t]] = f0
            box[R.Wee[b, t]] = iand(box[R.Wee[b, t]], isq(e0)); box[R.Wff[b, t]] = iand(box[R.Wff[b, t]], isq(f0))
            box[R.Wef[b, t]] = iand(box[R.Wef[b, t]], iprod(e0, f0))
        end
        g(sy, n) = R.Lv[sy][n, t]
        for l in net.lines
            n = l.name; i, j = l.from, l.to; C = (-Il(l), Il(l))
            for (sy, iv) in ((:Lee, iprod(ex[i], ex[j])), (:Lff, iprod(fx[i], fx[j])), (:Lef, iprod(ex[i], fx[j])), (:Lfe, iprod(fx[i], ex[j])),
                             (:lr, C), (:li, C), (:eilr, iprod(ex[i], C)), (:eili, iprod(ex[i], C)), (:filr, iprod(fx[i], C)), (:fili, iprod(fx[i], C)),
                             (:ejlr, iprod(ex[j], C)), (:ejli, iprod(ex[j], C)), (:fjlr, iprod(fx[j], C)), (:fjli, iprod(fx[j], C)),
                             (:lrlr, isq(C)), (:lili, isq(C)), (:lrli, iprod(C, C)))
                v = g(sy, n); box[v] = iand(box[v], iv)
            end
        end
        for ((bt, k1, k2, sy), v) in R.xcur
            tt = parse(Int, split(bt, ',')[2]); tt == t || continue
            c1 = vb.C[k1]; c2 = vb.C[k2]
            box[v] = iand(box[v], (-c1 * c2, c1 * c2))
        end
    end
    return radius
end

"""
    certify(net, delta_deg, kw; k, cuts = Symbol[])

Solve the relaxation and return U (kW, from the dual bound), the solver's primal and dual
objectives, the residuals, the PSD eigenvalue checks, the rank of the moment blocks, and, when
every block is numerically rank one, the recovered state for `validate_dispatch`.
"""
# Clarabel settings that solve the 24-period relaxation (the defaults stop on a numerical error there)
const CLARABEL_DEFAULTS = [("max_iter", 400), ("static_regularization_constant", 1e-7), ("tol_gap_abs", 1e-6), ("tol_gap_rel", 1e-6), ("tol_feas", 1e-6)]
function certify(net, delta_deg, kw; k, cuts = Symbol[], solver_attributes = CLARABEL_DEFAULTS, m_lo_kW = -1000.0, m_hi_kW = 1e4, band_in_margin = true, tighten = false, detail = 0, src_scale = false)
    t0 = time()
    R = relax_margin(net, delta_deg, kw; k, cuts, solver_attributes, m_lo_kW, m_hi_kW, band_in_margin, src_scale)
    tb = time() - t0
    optimize!(R.model)
    ts = time() - t0 - tb
    st = termination_status(R.model)
    S = R.S
    primal = has_values(R.model) ? -objective_value(R.model) * S : NaN
    dualobj = try -dual_objective_value(R.model) * S catch; NaN end
    box_loose = copy(R.box)
    radius = tighten ? tighten_boxes!(R, net, kw, delta_deg) : nothing
    tv0 = time()
    db = has_duals(R.model) ? dual_bound(R.model; box = R.box, detail) : (; bound = -Inf, psd_viol = NaN, sign_viol = NaN, box_term = NaN, rounding = NaN, unbounded_variable = nothing)
    U = -db.bound * S
    # Infeasible relaxation: validate the Farkas ray independently. With c = 0, dual_bound returns a
    # lower bound on 0 over every point satisfying the constraints and the box; > 0 proves there is none.
    farkas = (; bound = NaN, box_term = NaN, psd_viol = NaN, constant = NaN, top = [], csrc = Dict{String,Float64}(), bound_loose = NaN, box_term_loose = NaN)
    if st in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE) && dual_status(R.model) in (MOI.INFEASIBILITY_CERTIFICATE, MOI.NEARLY_INFEASIBILITY_CERTIFICATE)
        fb = dual_bound(R.model; c = zero(AffExpr), box = R.box, detail)
        fl = dual_bound(R.model; c = zero(AffExpr), box = box_loose)
        farkas = (; bound = fb.bound, box_term = fb.box_term, psd_viol = fb.psd_viol, constant = fb.constant, top = fb.top, csrc = fb.csrc,
                  bound_loose = fl.bound, box_term_loose = fl.box_term)
    end
    verify_s = time() - tv0
    pres = has_values(R.model) ? maximum(values(primal_feasibility_report(R.model; atol = 0.0)); init = 0.0) : NaN
    minev, rank_ratio = Inf, 0.0
    if has_values(R.model)
        for p in R.psd
            X = value.(constraint_object(p.c).func); n = constraint_object(p.c).set.side_dimension
            Mx = zeros(n, n); for j in 1:n, i in 1:j; Mx[i, j] = Mx[j, i] = X[tri(i, j)]; end
            ev = eigvals(Symmetric(Mx)); minev = min(minev, ev[1]); startswith(p.line, "bus:") || (rank_ratio = max(rank_ratio, ev[end-1] / max(ev[end], eps())))
        end
    end
    # where consistency is lost: line losses the relaxation books beyond what its own currents imply,
    # sum_l R_l (lrlr + lili - lr^2 - li^2), and |V|^2 beyond the squared first moments (kW and pu)
    fict, vgap = NaN, NaN
    if has_values(R.model)
        fict = 0.0; vgap = 0.0
        for l in net.lines, t in 1:kw.T
            Rl = real(-net.Y_base / l.Y[1, 2]); g(sy) = value(variable_by_name(R.model, string(sy) * "[" * l.name * "," * string(t) * "]"))
            fict += Rl * (g(:lrlr) + g(:lili) - g(:lr)^2 - g(:li)^2)
        end
        for b in net.buses, t in 1:kw.T
            vgap = max(vgap, value(R.Wee[b, t]) + value(R.Wff[b, t]) - value(R.e[b, t])^2 - value(R.f[b, t])^2)
        end
        fict *= S / kw.T
    end
    rec = nothing
    if has_values(R.model)
        T = kw.T; Ts = 1:T
        V = Dict(b => [complex(value(R.e[b, t]), value(R.f[b, t])) for t in Ts] for b in net.buses)
        lv(sy, n, t) = value(variable_by_name(R.model, string(sy) * "[" * n * "," * string(t) * "]"))
        raw = (; lr = Dict(l.name => [lv(:lr, l.name, t) for t in Ts] for l in net.lines), li = Dict(l.name => [lv(:li, l.name, t) for t in Ts] for l in net.lines))
        rec = (; raw, V_pu = V, P_B_kW = Dict(kk => [value(R.P_B[kk, t]) * S for t in Ts] for kk in axes(R.P_B, 1)),
               Q_PV_kvar = Dict(kk => [value(R.Q_PV[kk, t]) * S for t in Ts] for kk in axes(R.Q_PV, 1)),
               B_kWh = Dict(kk => [value(R.B[kk, t]) * S for t in 0:T] for kk in axes(R.B, 1)),
               P_subs_kW = Dict(s.name => fill(NaN, T) for s in net.sources), Q_subs_kvar = Dict(s.name => fill(NaN, T) for s in net.sources))
    end
    return (; status = st, U, primal, dualobj, gap = primal - dualobj, primal_residual = pres, psd_primal_min_eig = minev,
            psd_dual_violation = db.psd_viol, dual_sign_violation = db.sign_viol, box_term_kW = db.box_term * S,
            rounding_kW = db.rounding * S, rank_ratio, recovered = rec, build_s = tb, solve_s = ts,
            nblocks = length(R.psd), nvars = num_variables(R.model), R, farkas, verify_s, dual_status = dual_status(R.model), db, radius, box_loose, fictitious_loss_kW = fict, vsq_gap = vgap)
end

"""
    export_certificate(model, box, path; c = objective_function(model), note = "")

Write everything an independent check of a solved relaxation needs, so that the bound can be recomputed
without another solve (root_level/verify_certificate.jl): the objective c, every variable's box (the solver
bounds intersected with the certificate box), every scalar constraint (sense, right-hand side, constant,
terms, and the solver's dual value) and every PSD block (side, the dual triangle, and each triangle entry's
constant and terms). Numbers are written in Julia's shortest round-trip form, so they parse back exactly.
With c = 0 (zero(AffExpr)) the file is a Farkas check.
"""
function export_certificate(model, box, path; c = objective_function(model), note = "")
    c = convert(AffExpr, c)
    vars = all_variables(model); vid = Dict(v => i for (i, v) in enumerate(vars))
    open(path, "w") do io
        println(io, "# relaxation certificate data; ", note)
        print(io, "OBJ ", c.constant, " ", length(c.terms)); for (v, a) in c.terms; print(io, " ", vid[v], " ", a); end; println(io)
        for (i, v) in enumerate(vars)
            l, u = is_fixed(v) ? (fix_value(v), fix_value(v)) : (has_lower_bound(v) ? lower_bound(v) : -Inf, has_upper_bound(v) ? upper_bound(v) : Inf)
            haskey(box, v) && (l = max(l, box[v][1]); u = min(u, box[v][2]))
            println(io, "VAR ", i, " ", l, " ", u, " ", replace(name(v), ' ' => '_'))
        end
        bk = backend(model)
        for (F, Sset) in list_of_constraint_types(model)
            F == VariableRef && continue
            for cref in all_constraints(model, F, Sset)
                obj = constraint_object(cref)
                if obj.set isa MOI.PositiveSemidefiniteConeTriangle
                    lamv = MOI.get(bk, MOI.ConstraintDual(), index(cref)); fv = moi_function(obj); n = obj.set.side_dimension
                    print(io, "PSD ", n); for x in lamv; print(io, " ", x); end; println(io)
                    terms = [Tuple{Int,Float64}[] for _ in 1:length(fv.constants)]
                    for term in fv.terms
                        push!(terms[term.output_index], (vid[VariableRef(model, term.scalar_term.variable)], term.scalar_term.coefficient))
                    end
                    for (q, cst) in enumerate(fv.constants)
                        print(io, "E ", q, " ", cst, " ", length(terms[q])); for (v, a) in terms[q]; print(io, " ", v, " ", a); end; println(io)
                    end
                elseif obj.set isa MOI.SecondOrderCone
                    yv = MOI.get(bk, MOI.ConstraintDual(), index(cref)); fv = moi_function(obj)
                    print(io, "SOC ", length(yv)); for x in yv; print(io, " ", x); end; println(io)
                    terms = [Tuple{Int,Float64}[] for _ in 1:length(fv.constants)]
                    for term in fv.terms
                        push!(terms[term.output_index], (vid[VariableRef(model, term.scalar_term.variable)], term.scalar_term.coefficient))
                    end
                    for (q, cst) in enumerate(fv.constants)
                        print(io, "E ", q, " ", cst, " ", length(terms[q])); for (v, a) in terms[q]; print(io, " ", v, " ", a); end; println(io)
                    end
                else
                    fn = obj.func; sense = obj.set isa MOI.GreaterThan ? "GE" : obj.set isa MOI.LessThan ? "LE" : obj.set isa MOI.EqualTo ? "EQ" : error("unsupported set $(typeof(obj.set))")
                    print(io, "LIN ", sense, " ", MOI.constant(obj.set), " ", dual(cref), " ", fn.constant, " ", length(fn.terms))
                    for (v, a) in fn.terms; print(io, " ", vid[v], " ", a); end; println(io)
                end
            end
        end
    end
    return path
end
