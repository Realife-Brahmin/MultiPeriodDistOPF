# full_angle_pf.jl -- full-angle AC power flow for single-phase OpenDSS decks, posed
# as a JuMP optimization (Ipopt or Gurobi) and cross-checked against OpenDSS.
#
# WHY "FULL ANGLE"
# ----------------
# root_level/small2poi_mpopf.jl uses the angle-relaxed branch-flow model: angles never
# enter the optimization and are recovered afterwards (beta = angle(v_i - z* S_ij)), so
# the split between substations is the optimizer's choice. Here every bus voltage is
# complex, V = e + jf, so angles are carried exactly. Fixing each source's magnitude AND
# angle then makes the split a consequence of physics.
#
# With every source angle fixed there is no degree of freedom left: the equality
# constraints are the square power-flow equations. On a network with one load bus (like
# small2poi) those have exactly two roots -- the high-voltage one, which OpenDSS converges
# to, and a low-voltage one. The objective (minimum total substation power, i.e. minimum
# losses) is only there to select the high-voltage root.
#
# SOURCES
# -------
# A Vsource is an EMF E = pu * cis(delta) behind its series impedance Zs, written in
# impedance form, V_bus = E - Zs * I_source. That form is exact for any Zs including zero,
# so the same model covers the deck as written and a "stiff" variant, where Zs is
# OpenDSS's 1e-8 ohm stand-in for an ideal source (the form ieee2522C/large10kC use).
# `ideal_sources = true` forces Zs = 0 exactly.
#
# DATA
# ----
# The network is read from OpenDSS's own primitive admittance matrices (YPrim), so the
# JuMP model and OpenDSS use identical impedances. Two consequences a literal reading of
# the deck would get wrong:
#   * a 1-phase Vsource's series impedance is Zs = (2 Z1 + Z0)/3, not Z1;
#   * every Line carries OpenDSS's default shunt charging (C1 = 3.4 nF per unit length).
# Supported: single-conductor Lines, 1-phase Vsources, constant-PQ Loads (model=1). Any
# other enabled element is an error, never silently dropped.
#
# PER UNIT
# --------
# V_base = the sources' basekv (line-to-neutral for phases=1); S_base = 1000 kVA.
# OpenDSS's own per-unit voltages are never used -- for small2poi_1ph they are sqrt(3)
# off (see check_voltage_bases.jl). Cross-checks compare volts, kW and kvar.

using OpenDSSDirect
using JuMP
using Ipopt
using Gurobi
using Printf

const ODD = OpenDSSDirect
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

"Bus name without its node qualifier: \"1s.1\" -> \"1s\"."
busname(bus) = lowercase(String(split(bus, '.')[1]))

"Node qualifiers of a bus reference: \"1s.1.0\" -> [1, 0]; \"1s\" -> []."
nodes_of(bus) = [parse(Int, n) for n in split(bus, '.')[2:end]]

"""
    compile_deck(system; stiff_sources = false)

Compile `rawData/<system>/Master.dss`. `stiff_sources = true` replaces every Vsource's
impedance with a 1e-8 ohm reactance in memory; the deck on disk is never modified.
OpenDSS's convergence tolerance is tightened from its default 1e-4 so that the
cross-check measures the models, not OpenDSS's stopping rule.
"""
function compile_deck(system; stiff_sources = false, tolerance = 1e-12, max_iterations = 1000)
    master = joinpath(REPO_ROOT, "rawData", system, "Master.dss")
    isfile(master) || error("Master.dss not found: $master")
    ODD.Text.Command("Clear")
    ODD.Text.Command("Redirect \"$master\"")
    if stiff_sources
        for el in ODD.Circuit.AllElementNames()
            startswith(lowercase(el), "vsource.") || continue
            ODD.Circuit.SetActiveElement(el)
            ODD.CktElement.Enabled() || continue
            ODD.Text.Command("Edit $el R1=0 X1=0.00000001 R0=0 X0=0.00000001")
        end
    end
    ODD.Text.Command("Set tolerance=$tolerance")
    ODD.Text.Command("Set maxiterations=$max_iterations")
    # An Edit marks YPrim stale; OpenDSS rebuilds it only when it next builds the system
    # Y matrix. Solve once so read_network() sees the impedances actually in use.
    ODD.Solution.Solve()
    return master
end

"""
    read_network(; S_base_kVA = 1000.0)

Read the compiled circuit into the data `solve_full_angle` needs. Impedances and
admittances stay in ohms and siemens here; `Y_base` converts them to per unit.
"""
function read_network(; S_base_kVA = 1000.0)
    lines, sources, loads = [], [], []
    for el in ODD.Circuit.AllElementNames()
        ODD.Circuit.SetActiveElement(el)
        ODD.CktElement.Enabled() || continue
        class, name = lowercase.(split(el, '.'; limit = 2))
        refs = ODD.CktElement.BusNames()
        all(n -> n in (0, 1), reduce(vcat, nodes_of.(refs))) ||
            error("$el connects to a node other than 1 or ground; only single-phase decks are supported")

        if class == "line"
            ODD.CktElement.NumConductors() == 1 || error("$el: only single-conductor lines are supported")
            push!(lines, (; name, from = busname(refs[1]), to = busname(refs[2]),
                          Y = ODD.CktElement.YPrim()))
        elseif class == "vsource"
            ODD.CktElement.NumPhases() == 1 || error("$el: only 1-phase Vsources are supported")
            z_series = 1 / ODD.CktElement.YPrim()[1, 1]   # terminal 1 = its bus, terminal 2 = ground
            ODD.Vsources.Name(name)
            push!(sources, (; name, bus = busname(refs[1]), kV = ODD.Vsources.BasekV(),
                            pu = ODD.Vsources.PU(), z_series))
        elseif class == "load"
            ODD.Loads.Name(name)
            Int(ODD.Loads.Model()) == 1 || error("$el: only constant-PQ loads (model=1) are supported")
            push!(loads, (; name, bus = busname(refs[1]), kW = ODD.Loads.kW(), kvar = ODD.Loads.kvar(),
                          kV = ODD.Loads.kV(), vminpu = ODD.Loads.Vminpu(), vmaxpu = ODD.Loads.Vmaxpu()))
        else
            error("$el: element class '$class' is not supported by the full-angle model")
        end
    end

    isempty(sources) && error("no enabled Vsource")
    allunique(s.bus for s in sources) || error("two Vsources share a bus; not supported")
    kV_bases = unique(s.kV for s in sources)
    length(kV_bases) == 1 || error("Vsources disagree on basekv: $kV_bases")
    kV_base = only(kV_bases)

    return (; buses = lowercase.(ODD.Circuit.AllBusNames()), lines, sources, loads,
            kV_base, S_base_kVA, Y_base = S_base_kVA * 1e3 / (kV_base * 1e3)^2)
end

# One Gurobi environment for every solve, so the license is read once.
const GRB_ENV = Ref{Gurobi.Env}()
gurobi_env() = isassigned(GRB_ENV) ? GRB_ENV[] : (GRB_ENV[] = Gurobi.Env(; output_flag = 0))

"""
    gurobi_usable()

`(true, "")` if a Gurobi environment starts, else `(false, reason)` -- so a missing or
expired license skips the Gurobi solves with a message instead of aborting a sweep.
"""
function gurobi_usable()
    try
        gurobi_env()
        return (true, "")
    catch err
        return (false, sprint(showerror, err))
    end
end

"Result-shaped placeholder for a solver that did not run: status :SKIPPED, every value NaN."
skipped(net, solver) = (; solver, status = :SKIPPED, ok = false, time = NaN,
                        P_subs_kW = Dict(s.name => NaN for s in net.sources),
                        Q_subs_kvar = Dict(s.name => NaN for s in net.sources),
                        V_pu = Dict(b => complex(NaN, NaN) for b in net.buses))

# Ipopt's tolerance: in per unit on 1000 kVA this network's admittances run to ~5e3, so the
# power-balance residual bottoms out near 1e-12 in floating point. Measured on the small2poi
# sweep: 1e-9 converges cleanly at every point for all three source variants (deck,
# stand-in, ideal); 1e-10 leaves ideal sources "almost solved" at 2 of 33. 1e-9 pu = 1e-6 kW.
const IPOPT_TOL = 1e-9

function new_model(solver)
    if solver == :ipopt
        model = Model(Ipopt.Optimizer)
        set_attribute(model, "tol", IPOPT_TOL)
        set_attribute(model, "constr_viol_tol", IPOPT_TOL)
        set_attribute(model, "max_iter", 500)
    elseif solver == :gurobi
        model = Model(() -> Gurobi.Optimizer(gurobi_env()))
        set_attribute(model, "NonConvex", 2)          # spatial branch-and-bound: global optimum
        set_attribute(model, "FeasibilityTol", 1e-9)  # Gurobi's tightest
        set_attribute(model, "OptimalityTol", 1e-9)
        set_attribute(model, "MIPGap", 1e-9)
        set_attribute(model, "TimeLimit", 60.0)
    else
        error("solver must be :ipopt or :gurobi, got :$solver")
    end
    set_silent(model)
    return model
end

"""
    solve_full_angle(net, delta_deg; solver = :ipopt, ideal_sources = false,
                     sense = MIN_SENSE, V_load_start = nothing)

Full-angle AC power flow with each source's EMF fixed at `pu * cis(delta)`, where
`delta_deg[name]` is that Vsource's angle in degrees, behind the source impedance OpenDSS
reports (`ideal_sources = true`: behind none).

Substation power is what OpenDSS reports for a Vsource: power delivered into the network
at its terminal bus.

The other (low-voltage) power-flow root: Gurobi finds it globally with
`sense = MAX_SENSE`; Ipopt, being local, converges to whichever root is nearest its
start, so pass `V_load_start` (per unit) to start every bus that carries load low.
"""
function solve_full_angle(net, delta_deg; solver = :ipopt, ideal_sources = false,
                          sense = MIN_SENSE, V_load_start = nothing)
    model = new_model(solver)
    buses = net.buses
    names = [s.name for s in net.sources]
    E = Dict(s.name => s.pu * cis(deg2rad(delta_deg[s.name])) for s in net.sources)

    @variable(model, -1.5 <= e[buses] <= 1.5)     # bus voltage V = e + jf, per unit
    @variable(model, -1.5 <= f[buses] <= 1.5)
    @variable(model, -1e5 <= ir[names] <= 1e5)    # source current I = ir + j*ii, per unit
    @variable(model, -1e5 <= ii[names] <= 1e5)

    # Flat start rotated to the sources; buses carrying load optionally started low
    V0 = sum(values(E)) / length(E)
    load_buses = Set(d.bus for d in net.loads)
    for b in buses
        V = (V_load_start !== nothing && b in load_buses) ? V_load_start * cis(angle(V0)) : V0
        set_start_value(e[b], real(V))
        set_start_value(f[b], imag(V))
    end

    # Each source: V_bus = E - Zs * I, i.e. e = Re(E) - (R ir - X ii), f = Im(E) - (R ii + X ir).
    # With Zs = 0 that pins the bus at E, so fix it outright rather than constrain it.
    for s in net.sources
        R, X = ideal_sources ? (0.0, 0.0) : reim(s.z_series * net.Y_base)   # ohm -> per unit
        if R == 0 && X == 0
            fix(e[s.bus], real(E[s.name]); force = true)
            fix(f[s.bus], imag(E[s.name]); force = true)
        else
            @constraint(model, e[s.bus] == real(E[s.name]) - (R * ir[s.name] - X * ii[s.name]))
            @constraint(model, f[s.bus] == imag(E[s.name]) - (R * ii[s.name] + X * ir[s.name]))
        end
    end

    # Power leaving bus i into a line whose primitive row at i is (Yii, Yij):
    #   S = V_i * conj(Yii * V_i + Yij * V_j)
    function flow_out(i, j, Yii, Yij)
        g1, b1 = reim(Yii)
        g2, b2 = reim(Yij)
        vv = e[i]^2 + f[i]^2
        c = e[i] * e[j] + f[i] * f[j]             # Re(V_i * conj(V_j))
        s = f[i] * e[j] - e[i] * f[j]             # Im(V_i * conj(V_j))
        return (P = g1 * vv + g2 * c + b2 * s, Q = -b1 * vv + g2 * s - b2 * c)
    end
    out = Dict(b => [] for b in buses)
    for l in net.lines
        Y = l.Y ./ net.Y_base
        push!(out[l.from], flow_out(l.from, l.to, Y[1, 1], Y[1, 2]))
        push!(out[l.to], flow_out(l.to, l.from, Y[2, 2], Y[2, 1]))
    end

    # Power delivered by each source at its bus: S = V conj(I)
    @variable(model, P_subs[names])
    @variable(model, Q_subs[names])
    for s in net.sources
        b, k = s.bus, s.name
        @constraint(model, P_subs[k] == e[b] * ir[k] + f[b] * ii[k])
        @constraint(model, Q_subs[k] == f[b] * ir[k] - e[b] * ii[k])
    end

    # Power balance at every bus: out through lines + load = delivered by a source there
    load_pu = Dict(b => 0.0im for b in buses)
    for d in net.loads
        load_pu[d.bus] += complex(d.kW, d.kvar) / net.S_base_kVA
    end
    src_at = Dict(s.bus => s.name for s in net.sources)
    for b in buses
        P_in = haskey(src_at, b) ? P_subs[src_at[b]] : 0.0
        Q_in = haskey(src_at, b) ? Q_subs[src_at[b]] : 0.0
        @constraint(model, sum(x.P for x in out[b]) + real(load_pu[b]) == P_in)
        @constraint(model, sum(x.Q for x in out[b]) + imag(load_pu[b]) == Q_in)
    end
    @objective(model, sense, sum(P_subs))

    optimize!(model)
    status = termination_status(model)
    ok = status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) && primal_status(model) == MOI.FEASIBLE_POINT
    val(x) = has_values(model) ? value(x) : NaN
    return (; solver, status, ok, time = solve_time(model),
            P_subs_kW = Dict(k => val(P_subs[k]) * net.S_base_kVA for k in names),
            Q_subs_kvar = Dict(k => val(Q_subs[k]) * net.S_base_kVA for k in names),
            V_pu = Dict(b => complex(val(e[b]), val(f[b])) for b in buses))
end

"""
    opendss_point(net, delta_deg)

Set each Vsource's angle, solve the compiled circuit, and read back what
`solve_full_angle` returns: voltages in per unit of net.kV_base, and substation power in
kW/kvar measured as the power leaving the source's bus through its lines, plus any load
on that bus. (The Vsource's own power reading, kept as `S_vsource_kVA`, is the same
quantity but loses digits when Zs is tiny: it evaluates Ys*(E - V) with |Ys| ~ 1e8 S.)

`pq_band` is false if any load left its constant-PQ band: outside [Vminpu, Vmaxpu]
OpenDSS turns a model=1 load into a constant impedance, and the two models then
legitimately differ.
"""
function opendss_point(net, delta_deg)
    for s in net.sources
        ODD.Vsources.Name(s.name)
        ODD.Vsources.AngleDeg(Float64(delta_deg[s.name]))
    end
    ODD.Solution.Solve()

    V = Dict{String,ComplexF64}()
    for b in net.buses
        ODD.Circuit.SetActiveBus(b)
        V[b] = first(ODD.Bus.Voltages()) / (net.kV_base * 1e3)
    end

    P, Q, S_vsource = Dict{String,Float64}(), Dict{String,Float64}(), Dict{String,ComplexF64}()
    for s in net.sources
        S = 0.0im
        for l in net.lines
            (l.from == s.bus || l.to == s.bus) || continue
            ODD.Circuit.SetActiveElement("Line." * l.name)
            pw = ODD.CktElement.Powers()                 # into the element, per terminal
            S += l.from == s.bus ? pw[1] : pw[2]
        end
        for d in net.loads
            d.bus == s.bus || continue
            ODD.Circuit.SetActiveElement("Load." * d.name)
            S += sum(ODD.CktElement.Powers())
        end
        P[s.name], Q[s.name] = real(S), imag(S)
        ODD.Circuit.SetActiveElement("Vsource." * s.name)
        S_vsource[s.name] = -ODD.CktElement.Powers()[1]  # into the element -> delivered
    end
    pq_band = all(d -> d.vminpu <= abs(V[d.bus]) * net.kV_base / d.kV <= d.vmaxpu, net.loads)
    return (; converged = ODD.Solution.Converged(), iterations = ODD.Solution.Iterations(),
            P_subs_kW = P, Q_subs_kvar = Q, V_pu = V, S_vsource_kVA = S_vsource, pq_band)
end

"Largest disagreement between two solutions: substation P (kW), Q (kvar), bus-voltage phasor (pu)."
function max_diff(a, b)
    return (; dP = maximum(abs(a.P_subs_kW[k] - b.P_subs_kW[k]) for k in keys(a.P_subs_kW)),
            dQ = maximum(abs(a.Q_subs_kvar[k] - b.Q_subs_kvar[k]) for k in keys(a.Q_subs_kvar)),
            dV = maximum(abs(a.V_pu[k] - b.V_pu[k]) for k in keys(a.V_pu)))
end
