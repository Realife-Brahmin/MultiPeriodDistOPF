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
# SOURCES AND LINES
# -----------------
# A Vsource is an EMF E = pu * cis(delta) behind its series impedance Zs, written in
# impedance form, V_bus = E - Zs * I_source. That form is exact for any Zs including zero,
# so the same model covers the deck as written and a "stiff" variant, where Zs is
# OpenDSS's 1e-8 ohm stand-in for an ideal source (the form ieee2522C/large10kC use).
# `ideal_sources = true` forces Zs = 0 exactly. Lines are in impedance form as well, with
# an explicit series current each, so 1e-6 ohm jumpers do not wreck the scaling.
# Rectangular voltages and currents keep every constraint linear or bilinear: the model
# is a QCQP, which Gurobi can solve to global optimality.
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
using LinearAlgebra
using Printf

const ODD = OpenDSSDirect
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

"Bus name without its node qualifier: \"1s.1\" -> \"1s\"."
busname(bus) = lowercase(String(split(bus, '.')[1]))

"Node qualifiers of a bus reference: \"1s.1.0\" -> [1, 0]; \"1s\" -> []."
nodes_of(bus) = [parse(Int, n) for n in split(bus, '.')[2:end]]

"""
    compile_deck(system; extra = String[], sources = nothing, disable = String[],
                 snapshot = false, stiff_sources = false, load_band = nothing)

Compile `rawData/<system>/Master.dss`, then adjust it in memory -- the deck on disk is
never modified:
  * `extra`: DSS commands run right after the deck, e.g. a test battery to add;
  * `sources`: the Vsources to keep (by name); every other Vsource is disabled;
  * `disable`: element classes to switch off wholesale, e.g. ["PVSystem", "Storage"];
  * `snapshot`: snapshot mode, i.e. loads at their base kW with no loadshape applied;
  * `stiff_sources`: replace each kept Vsource's impedance with a 1e-8 ohm reactance;
  * `load_band = (lo, hi)`: set every load's and storage element's [Vminpu, Vmaxpu].
    OpenDSS keeps a model=1 load constant-PQ only inside that band and makes it a
    constant impedance outside, while the JuMP model is constant-PQ everywhere.
OpenDSS's convergence tolerance is tightened from its default 1e-4 so that the
cross-check measures the models, not OpenDSS's stopping rule.

Returns `deck_load_band`, each load's (Vminpu, Vmaxpu) as the deck has it, so callers can
still report against the deck's own band after widening it.
"""
function compile_deck(system; extra = String[], sources = nothing, disable = String[],
                      snapshot = false, stiff_sources = false, load_band = nothing,
                      tolerance = 1e-12, max_iterations = 1000)
    master = joinpath(REPO_ROOT, "rawData", system, "Master.dss")
    isfile(master) || error("Master.dss not found: $master")
    ODD.Text.Command("Clear")
    ODD.Text.Command("Redirect \"$master\"")
    foreach(ODD.Text.Command, extra)
    for class in disable
        ODD.Text.Command("BatchEdit $class..* enabled=false")
    end
    deck_load_band = Dict{String,Tuple{Float64,Float64}}()
    for name in ODD.Loads.AllNames()
        lowercase(name) == "none" && continue
        ODD.Loads.Name(name)
        deck_load_band[lowercase(name)] = (ODD.Loads.Vminpu(), ODD.Loads.Vmaxpu())
    end
    # Below Vlowpu (default 0.5) OpenDSS makes a load constant-Z whatever Vminpu says, so a
    # band reaching below 0.5 pu lowers Vlowpu with it.
    if load_band !== nothing
        lo, hi = load_band
        ODD.Text.Command("BatchEdit Load..* Vminpu=$lo Vmaxpu=$hi Vlowpu=$(min(0.5, lo))")
        any(startswith(lowercase(el), "storage.") for el in ODD.Circuit.AllElementNames()) &&
            ODD.Text.Command("BatchEdit Storage..* Vminpu=$lo Vmaxpu=$hi")
    end
    kept = String[]
    for el in ODD.Circuit.AllElementNames()
        startswith(lowercase(el), "vsource.") || continue
        ODD.Circuit.SetActiveElement(el)
        ODD.CktElement.Enabled() || continue
        name = lowercase(split(el, '.')[2])
        if sources !== nothing && !(name in lowercase.(sources))
            ODD.Text.Command("Edit $el enabled=false")
            continue
        end
        push!(kept, name)
        stiff_sources && ODD.Text.Command("Edit $el R1=0 X1=0.00000001 R0=0 X0=0.00000001")
    end
    sources === nothing || Set(kept) == Set(lowercase.(sources)) ||
        error("asked for Vsources $(sources); enabled in the deck: $kept")
    snapshot && ODD.Text.Command("Set mode=Snapshot")
    ODD.Text.Command("Set tolerance=$tolerance")
    ODD.Text.Command("Set maxiterations=$max_iterations")
    # An Edit marks YPrim stale; OpenDSS rebuilds it only when it next builds the system
    # Y matrix. Solve once so read_network() sees the impedances actually in use.
    ODD.Solution.Solve()
    return (; master, deck_load_band)
end

"""
    read_network(; S_base_kVA = 1000.0)

Read the compiled circuit into the data `solve_full_angle` needs. Impedances and
admittances stay in ohms and siemens here; `Y_base` converts them to per unit.
"""
function read_network(; S_base_kVA = 1000.0)
    lines, sources, loads, batteries = [], [], [], []
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
        elseif class == "storage"
            # A real-power resource here: output P_B within +-kWrated, no reactive power.
            ODD.CktElement.NumPhases() == 1 || error("$el: only 1-phase storage is supported")
            prop(p) = parse(Float64, ODD.Properties.Value(p))
            push!(batteries, (; name, bus = busname(refs[1]), kW_rated = prop("kWrated"),
                              kVA = prop("kVA"), kWh_rated = prop("kWhrated")))
        else
            error("$el: element class '$class' is not supported by the full-angle model")
        end
    end

    isempty(sources) && error("no enabled Vsource")
    allunique(s.bus for s in sources) || error("two Vsources share a bus; not supported")
    kV_bases = unique(s.kV for s in sources)
    length(kV_bases) == 1 || error("Vsources disagree on basekv: $kV_bases")
    kV_base = only(kV_bases)

    return (; buses = lowercase.(ODD.Circuit.AllBusNames()), lines, sources, loads, batteries,
            kV_base, S_base_kVA, Y_base = S_base_kVA * 1e3 / (kV_base * 1e3)^2)
end

"""
    source_separation(net)

Series impedance (ohm) between every pair of source buses, keyed by (name_a, name_b)
with name_a < name_b: the Thevenin impedance between the two buses through the lines
alone (no shunts, loads or source impedance). On a radial network it is simply the
impedance of the path between them -- how electrically far apart two substations sit.
"""
function source_separation(net)
    idx = Dict(b => i for (i, b) in enumerate(net.buses))
    n = length(net.buses)
    L = zeros(ComplexF64, n, n)                  # series-admittance Laplacian
    for l in net.lines
        y, i, j = -l.Y[1, 2], idx[l.from], idx[l.to]
        L[i, i] += y; L[j, j] += y; L[i, j] -= y; L[j, i] -= y
    end
    Z = zeros(ComplexF64, n, n)
    Z[1:n-1, 1:n-1] = inv(L[1:n-1, 1:n-1])      # last bus as the reference node
    zeff(a, b) = Z[a, a] + Z[b, b] - 2 * Z[a, b]
    return Dict((s.name, t.name) => zeff(idx[s.bus], idx[t.bus])
                for s in net.sources, t in net.sources if s.name < t.name)
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
skipped(net, solver) = (; solver, status = :SKIPPED, ok = false, time = NaN, iterations = missing,
                        P_subs_kW = Dict(s.name => NaN for s in net.sources),
                        Q_subs_kvar = Dict(s.name => NaN for s in net.sources),
                        P_B_kW = Dict(bt.name => NaN for bt in net.batteries),
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
                     sense = MIN_SENSE, V_load_start = nothing,
                     battery = :idle, no_backflow = false, V_guard = nothing)

Full-angle AC power flow with each source's EMF fixed at `pu * cis(delta)`, where
`delta_deg[name]` is that Vsource's angle in degrees, behind the source impedance OpenDSS
reports (`ideal_sources = true`: behind none).

Substation power is what OpenDSS reports for a Vsource: power delivered into the network
at its terminal bus.

The other (low-voltage) power-flow root: Gurobi finds it globally with
`sense = MAX_SENSE`; Ipopt, being local, converges to whichever root is nearest its
start, so pass `V_load_start` (per unit) to start every bus that carries load low.

Batteries (`net.batteries`) inject P_B, discharging positive -- the repo's convention --
with no reactive power. `battery` is `:idle` (P_B = 0), a Dict of fixed outputs in kW by
battery name, or `:optimize`: P_B free within +-kWrated, and the objective becomes the
smallest dispatch, min sum(P_B^2). That is the OPF of interest with `no_backflow = true`,
which adds P_Subs >= 0 at every source. `V_guard` (per unit) bounds every bus voltage
from below, which keeps an optimizer off the low-voltage root -- where both substations
import hugely -- as a spurious way to meet P_Subs >= 0.
"""
function solve_full_angle(net, delta_deg; solver = :ipopt, ideal_sources = false,
                          sense = MIN_SENSE, V_load_start = nothing,
                          battery = :idle, no_backflow = false, V_guard = nothing)
    model = new_model(solver)
    buses = net.buses
    names = [s.name for s in net.sources]
    E = Dict(s.name => s.pu * cis(deg2rad(delta_deg[s.name])) for s in net.sources)

    bnames = [bt.name for bt in net.batteries]
    @variable(model, P_B[bnames])                 # battery output, discharging > 0, per unit
    for bt in net.batteries
        if battery === :optimize
            set_lower_bound(P_B[bt.name], -bt.kW_rated / net.S_base_kVA)
            set_upper_bound(P_B[bt.name], bt.kW_rated / net.S_base_kVA)
        else
            kW = battery === :idle ? 0.0 : battery[bt.name]
            fix(P_B[bt.name], kW / net.S_base_kVA; force = true)
        end
    end

    lines = [l.name for l in net.lines]
    @variable(model, -1.5 <= e[buses] <= 1.5)     # bus voltage V = e + jf, per unit
    @variable(model, -1.5 <= f[buses] <= 1.5)
    @variable(model, -1e5 <= ir[names] <= 1e5)    # source current into its bus, I = ir + j*ii
    @variable(model, -1e5 <= ii[names] <= 1e5)
    @variable(model, -1e5 <= lr[lines] <= 1e5)    # line series current, from -> to, I = lr + j*li
    @variable(model, -1e5 <= li[lines] <= 1e5)

    # Flat start rotated to the sources (every current then starts at zero). For the
    # low-voltage root, start the load buses low, the source buses at their EMF, and each
    # current at what KVL and KCL give for that start -- the heavy-current regime the low
    # root lives in, rather than zero current, which leads back to the high root.
    V0 = sum(values(E)) / length(E)
    Vs = Dict(b => V0 for b in buses)
    if V_load_start !== nothing
        for d in net.loads
            Vs[d.bus] = V_load_start * cis(angle(V0))
        end
        for s in net.sources
            Vs[s.bus] = E[s.name]
        end
        Is = Dict(l.name => (Vs[l.from] - Vs[l.to]) * l.Y[1, 2] / -net.Y_base for l in net.lines)
        for l in net.lines
            set_start_value(lr[l.name], real(Is[l.name]))
            set_start_value(li[l.name], imag(Is[l.name]))
        end
        for s in net.sources
            I = sum((l.from == s.bus ? 1 : -1) * Is[l.name] for l in net.lines if s.bus in (l.from, l.to))
            set_start_value(ir[s.name], real(I))
            set_start_value(ii[s.name], imag(I))
        end
    end
    for b in buses
        set_start_value(e[b], real(Vs[b]))
        set_start_value(f[b], imag(Vs[b]))
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

    # Each line: pi model read off its YPrim -- series Z = -1/Y12, shunt Y11 + Y12 and
    # Y22 + Y21 at the two ends. The series branch is in impedance form too,
    # V_from - V_to = Z * I, so near-zero jumpers (1e-6 ohm in ieee123_5poi_1ph, ~7e7 pu
    # as an admittance) stay well scaled.
    for l in net.lines
        R, X = reim(-net.Y_base / l.Y[1, 2])      # series impedance, per unit
        @constraint(model, e[l.from] - e[l.to] == R * lr[l.name] - X * li[l.name])
        @constraint(model, f[l.from] - f[l.to] == R * li[l.name] + X * lr[l.name])
    end

    # Power leaving bus b into line l: S = V_b * conj(sigma * I_l + y_sh * V_b), where
    # sigma = +1 at the from-end and -1 at the to-end, and y_sh is that end's shunt.
    function line_out(b, l, sigma, y_sh)
        g, bsh = reim(y_sh)
        vv = e[b]^2 + f[b]^2
        return (P = sigma * (e[b] * lr[l] + f[b] * li[l]) + g * vv,
                Q = sigma * (f[b] * lr[l] - e[b] * li[l]) - bsh * vv)
    end
    out = Dict(b => [] for b in buses)
    for l in net.lines
        Y = l.Y ./ net.Y_base
        push!(out[l.from], line_out(l.from, l.name, 1, Y[1, 1] + Y[1, 2]))
        push!(out[l.to], line_out(l.to, l.name, -1, Y[2, 2] + Y[2, 1]))
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
    # plus the output of any battery there
    load_pu = Dict(b => 0.0im for b in buses)
    for d in net.loads
        load_pu[d.bus] += complex(d.kW, d.kvar) / net.S_base_kVA
    end
    src_at = Dict(s.bus => s.name for s in net.sources)
    batt_at = Dict(b => [bt.name for bt in net.batteries if bt.bus == b] for b in buses)
    for b in buses
        P_in = (haskey(src_at, b) ? P_subs[src_at[b]] : 0.0) + sum(P_B[k] for k in batt_at[b]; init = 0.0)
        Q_in = haskey(src_at, b) ? Q_subs[src_at[b]] : 0.0
        @constraint(model, sum(x.P for x in out[b]) + real(load_pu[b]) == P_in)
        @constraint(model, sum(x.Q for x in out[b]) + imag(load_pu[b]) == Q_in)
    end

    no_backflow && @constraint(model, [k in names], P_subs[k] >= 0)
    V_guard === nothing || @constraint(model, [b in buses], e[b]^2 + f[b]^2 >= V_guard^2)
    if battery === :optimize
        @objective(model, Min, sum((P_B[k]^2 for k in bnames); init = zero(QuadExpr)))
    else
        @objective(model, sense, sum(P_subs))
    end

    optimize!(model)
    status = termination_status(model)
    ok = status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) && primal_status(model) == MOI.FEASIBLE_POINT
    val(x) = has_values(model) ? value(x) : NaN
    iterations = try barrier_iterations(model) catch; missing end
    return (; solver, status, ok, time = solve_time(model), iterations,
            P_subs_kW = Dict(k => val(P_subs[k]) * net.S_base_kVA for k in names),
            Q_subs_kvar = Dict(k => val(Q_subs[k]) * net.S_base_kVA for k in names),
            P_B_kW = Dict(k => val(P_B[k]) * net.S_base_kVA for k in bnames),
            V_pu = Dict(b => complex(val(e[b]), val(f[b])) for b in buses))
end

"""
    opendss_point(net, delta_deg; P_B = nothing)

Set each Vsource's angle and each battery's output (`P_B`, kW by name, discharging
positive; `nothing` = all zero), solve the compiled circuit, and read back what
`solve_full_angle` returns: voltages in per unit of net.kV_base, and substation power in
kW/kvar measured as the power leaving the source's bus through its lines, plus any load
on that bus, less any battery output there. (The Vsource's own power reading, kept as
`S_vsource_kVA`, is the same quantity but loses digits when Zs is tiny: it evaluates
Ys*(E - V) with |Ys| ~ 1e8 S.) `P_B_kW` is what each battery actually delivered.

A battery is dispatched with `Edit Storage.<name> kW=<P_B>`, which OpenDSS turns into
charging (P_B < 0) or discharging at exactly that terminal power; at P_B = 0 it idles
and draws %IdlingkW, which must therefore be 0 for the two models to match.

`pq_band` is false if any load left its constant-PQ band: outside [Vminpu, Vmaxpu]
OpenDSS turns a model=1 load into a constant impedance, and the two models then
legitimately differ.
"""
function opendss_point(net, delta_deg; P_B = nothing)
    for s in net.sources
        ODD.Vsources.Name(s.name)
        ODD.Vsources.AngleDeg(Float64(delta_deg[s.name]))
    end
    for bt in net.batteries
        ODD.Text.Command(@sprintf("Edit Storage.%s kW=%.12g", bt.name, P_B === nothing ? 0.0 : P_B[bt.name]))
    end
    ODD.Solution.Solve()

    V = Dict{String,ComplexF64}()
    for b in net.buses
        ODD.Circuit.SetActiveBus(b)
        V[b] = first(ODD.Bus.Voltages()) / (net.kV_base * 1e3)
    end
    S_B = Dict{String,ComplexF64}()
    for bt in net.batteries
        ODD.Circuit.SetActiveElement("Storage." * bt.name)
        S_B[bt.name] = -sum(ODD.CktElement.Powers())     # into the element -> delivered
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
        S -= sum((S_B[bt.name] for bt in net.batteries if bt.bus == s.bus); init = 0.0im)
        P[s.name], Q[s.name] = real(S), imag(S)
        ODD.Circuit.SetActiveElement("Vsource." * s.name)
        S_vsource[s.name] = -ODD.CktElement.Powers()[1]  # into the element -> delivered
    end
    pq_band = all(d -> d.vminpu <= abs(V[d.bus]) * net.kV_base / d.kV <= d.vmaxpu, net.loads)
    return (; converged = ODD.Solution.Converged(), iterations = ODD.Solution.Iterations(),
            P_subs_kW = P, Q_subs_kvar = Q, V_pu = V, S_vsource_kVA = S_vsource,
            P_B_kW = Dict(k => real(v) for (k, v) in S_B), pq_band)
end

"Largest disagreement between two solutions: substation P (kW), Q (kvar), bus-voltage phasor (pu)."
function max_diff(a, b)
    return (; dP = maximum(abs(a.P_subs_kW[k] - b.P_subs_kW[k]) for k in keys(a.P_subs_kW)),
            dQ = maximum(abs(a.Q_subs_kvar[k] - b.Q_subs_kvar[k]) for k in keys(a.Q_subs_kvar)),
            dV = maximum(abs(a.V_pu[k] - b.V_pu[k]) for k in keys(a.V_pu)))
end
