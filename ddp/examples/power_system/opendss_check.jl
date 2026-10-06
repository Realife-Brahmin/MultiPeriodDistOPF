# OpenDSS feasibility check of an MPOPF dispatch.
#
# The optimizers work on the relaxed branch-flow model and, for FilterDDP, may
# stop at a near-optimal point. This script hands only the DISPATCH (battery
# powers, PV reactive powers) to OpenDSS, solves the actual power flow of each
# period on the system's own DSS files, and reports how far the optimizer's
# network state is from it:
#   - voltages: largest difference, and buses outside their limits in OpenDSS;
#   - substation import and losses: optimizer against OpenDSS;
#   - the objective re-evaluated with OpenDSS's import.
#
# What comes from where, so the check stays independent of the OPF parser:
#   network, source           the DSS files (BranchData, Vsource);
#   loads                     the DSS Load elements, times the load shape;
#   PV real power             the DSS PVSystem ratings, times the PV shape;
#   device buses              the DSS PVSystem and Storage elements;
#   shapes, bases, limits     the exported instance (shapes are not in the DSS
#                             files: the instance resamples them with T);
#   battery powers, PV vars   the solution.
# PV systems and batteries are disabled and re-injected as constant-power
# elements at the dispatched values, and loads are held at constant power
# outside 0.95-1.05 pu as well, since that is the model the OPF solves. The
# source is set to 1.05 pu, the value the OPF fixes (the ieee123 DSS file says
# 1.03). Battery efficiency is not modelled by the OPF and plays no role here.
# Devices on the substation bus are skipped: the OPF's balance rows ignore them.
#
#   REDUCED_PROFILE=periodic REDUCED_CB=system TERMINAL_SOC_SOFT=1 \
#   julia --startup-file=no --project=envs/tadmm \
#         ddp/examples/power_system/opendss_check.jl <system> <T> <solution.jls> [label]

using LinearAlgebra
using OpenDSSDirect
using Printf
using Serialization

const ODD = OpenDSSDirect
const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "terminal_soc_penalty.jl"))
include(joinpath(@__DIR__, "control_layout.jl"))

system = ARGS[1]
T = parse(Int, ARGS[2])
solfile = ARGS[3]
label = length(ARGS) >= 4 ? ARGS[4] : basename(solfile)
ptag = haskey(ENV, "REDUCED_PROFILE") ? "_" * ENV["REDUCED_PROFILE"] : ""
data = deserialize(joinpath(REPO, "ddp", "results", "network_filterddp", "network_data_$(system)_T$(T)$(ptag).jls"))
haskey(ENV, "REDUCED_CB") && (data[:C_B] = battery_cb(system, ENV["REDUCED_CB"]))
gammaT = terminal_soc_soft() ? gamma_terminal(system) : 0.0
sol = deserialize(solfile)
u = sol[:u]
idx, nu = control_layout(data)
(length(u) == T && length(u[1]) == nu) || error("solution does not match $(system) T=$T")

buses, lines, bset, dset = data[:Nset], data[:Lset], data[:Bset], data[:Dset]
root, pbase, dt, kvb = data[:substationBus], data[:kVA_B], data[:delta_t_h], data[:kV_B]
buspos = Dict(j => k for (k, j) in enumerate(buses))
batpos = Dict(j => k for (k, j) in enumerate(bset))
derpos = Dict(j => k for (k, j) in enumerate(dset))
busnum(s) = parse(Int, split(strip(s), ".")[1])
ask(q) = strip(ODD.Text.Command("? " * q))

# ---- the circuit, from the system's own files ------------------------------
ODD.Text.Command("Clear")
ODD.Text.Command("Redirect \"$(joinpath(REPO, "rawData", system, "Master.dss"))\"")
source_pu_file = parse(Float64, ask("Vsource.source.pu"))
ODD.Text.Command("Edit Vsource.source pu=1.05")
ODD.Text.Command("Set mode=snapshot")
ODD.Text.Command("Set controlmode=off")
ODD.Text.Command("Set tolerance=1e-9")
ODD.Text.Command("Set maxiterations=200")
# OpenDSS gives every Line a default shunt capacitance (C1 = 3.4, C0 = 1.6 nF
# per unit length) unless the file sets it; the system files do not, and the
# OPF's branch-flow model has no shunt. OPENDSS_ZERO_LINE_C=1 removes it, to
# separate that modelling difference from everything else.
zero_c = get(ENV, "OPENDSS_ZERO_LINE_C", "0") != "0"
zero_c && ODD.Text.Command("BatchEdit Line..* c1=0 c0=0")
line_c1 = parse(Float64, ask("Line.$(first(ODD.Lines.AllNames())).c1"))

loads = Tuple{String,Float64,Float64}[]
for name in ODD.Loads.AllNames()
    ODD.Loads.Name(name)
    push!(loads, (name, ODD.Loads.kW(), ODD.Loads.kvar()))
end
ODD.Text.Command("BatchEdit Load..* Vminpu=0.3 Vmaxpu=1.7")

pvs = Tuple{Int,Float64,String}[]              # bus, Pmpp (kW), kV
for name in ODD.PVsystems.AllNames()
    name == "NONE" && continue
    push!(pvs, (busnum(ask("PVSystem.$name.bus1")), parse(Float64, ask("PVSystem.$name.Pmpp")), ask("PVSystem.$name.kv")))
end
bats = Tuple{Int,String}[]                      # bus, kV
for name in ODD.Storages.AllNames()
    name == "NONE" && continue
    push!(bats, (busnum(ask("Storage.$name.bus1")), ask("Storage.$name.kv")))
end
isempty(pvs) || ODD.Text.Command("BatchEdit PVSystem..* enabled=no")
isempty(bats) || ODD.Text.Command("BatchEdit Storage..* enabled=no")
skipped = count(p -> p[1] == root, pvs) + count(b -> b[1] == root, bats)
filter!(p -> p[1] != root, pvs); filter!(b -> b[1] != root, bats)
(sort(first.(pvs)) == sort(dset) || sort(first.(pvs)) == sort(filter(!=(root), dset))) ||
    error("PV buses in the DSS files differ from the instance")
(sort(first.(bats)) == sort(bset) || sort(first.(bats)) == sort(filter(!=(root), bset))) ||
    error("battery buses in the DSS files differ from the instance")
for (j, _, kv) in pvs
    ODD.Text.Command("New Load.injpv$(j) Bus1=$(j).1 Phases=1 Model=1 kV=$(kv) kW=0 kvar=0 Vminpu=0.3 Vmaxpu=1.7")
end
for (j, kv) in bats
    ODD.Text.Command("New Load.injb$(j) Bus1=$(j).1 Phases=1 Model=1 kV=$(kv) kW=0 kvar=0 Vminpu=0.3 Vmaxpu=1.7")
end

# ---- one power flow per period ----------------------------------------------
shapeL, shapePV, price = data[:LoadShapeLoad], data[:LoadShapePV], data[:LoadShapeCost]
qmaxpu(j, t) = sqrt(max(0.0, data[:S_D_R][j]^2 - data[:p_D_pu][j,t]^2))
vlo(j) = data[:Vminpu][j]; vhi(j) = data[:Vmaxpu][j]
VTOL = 1e-4                                    # pu; a limit is "violated" beyond this

ps_opf = zeros(T); ps_dss = zeros(T); loss_opf = zeros(T); loss_dss = zeros(T); fict = zeros(T)
dvmax = zeros(T); vmin_dss = zeros(T); vmax_dss = zeros(T); vmin_opf = zeros(T)
nlow = zeros(Int, T); nhigh = zeros(Int, T); worstlow = zeros(T); worsthigh = zeros(T)
pv_mismatch = 0.0; converged = true
for t in 1:T
    for (name, kw, kvar) in loads
        ODD.Text.Command("Edit Load.$name kW=$(kw * shapeL[t]) kvar=$(kvar * shapeL[t])")
    end
    for (j, pmpp, _) in pvs
        p = pmpp * shapePV[t]
        global pv_mismatch = max(pv_mismatch, abs(p - data[:p_D_pu][j,t] * pbase))
        q = qmaxpu(j, t) * u[t][idx.qnorm[derpos[j]]] * pbase
        ODD.Text.Command("Edit Load.injpv$(j) kW=$(-p) kvar=$(-q)")
    end
    for (j, _) in bats
        ODD.Text.Command("Edit Load.injb$(j) kW=$(-u[t][idx.pb[batpos[j]]] * pbase) kvar=0")
    end
    ODD.Solution.Solve()
    global converged &= ODD.Solution.Converged()

    vmag = ODD.Circuit.AllBusVMag() ./ (kvb * 1000)
    nodes = ODD.Circuit.AllNodeNames()
    vd = Dict(busnum(n) => vmag[i] for (i, n) in enumerate(nodes))
    vo = Dict(j => sqrt(u[t][idx.v[buspos[j]]]) for j in buses)
    dvmax[t] = maximum(abs(vd[j] - vo[j]) for j in buses)
    vmin_dss[t] = minimum(vd[j] for j in buses); vmax_dss[t] = maximum(vd[j] for j in buses if j != root)
    vmin_opf[t] = minimum(values(vo))
    for j in buses
        j == root && continue
        lo = vlo(j) - vd[j]; hi = vd[j] - vhi(j)
        lo > VTOL && (nlow[t] += 1; worstlow[t] = max(worstlow[t], lo))
        hi > VTOL && (nhigh[t] += 1; worsthigh[t] = max(worsthigh[t], hi))
    end
    ps_dss[t] = -real(ODD.Circuit.TotalPower()[1])          # kW into the feeder
    loss_dss[t] = real(ODD.Circuit.Losses()[1]) / 1000      # kW
    ps_opf[t] = u[t][idx.ps] * pbase
    for (k, e) in enumerate(lines)
        r = data[:rdict_pu][e]
        loss_opf[t] += r * u[t][idx.ell[k]] * pbase
        # loss the relaxation carries beyond the physical current
        fict[t] += r * (u[t][idx.ell[k]] - (u[t][idx.P[k]]^2 + u[t][idx.Q[k]]^2) / u[t][idx.v[buspos[e[1]]]]) * pbase
    end
end

# ---- objective with the optimizer's import and with OpenDSS's --------------
battery = sum(data[:C_B] * pbase^2 * dt * sum(u[t][k]^2 for k in idx.pb) for t in 1:T)
terminal = 0.0
if gammaT > 0
    for (b, j) in enumerate(bset)
        energy = data[:B0_pu][j] - dt * sum(u[t][idx.pb[b]] for t in 1:T)
        global terminal += gammaT * (energy - data[:B0_pu][j])^2
    end
end
J_opf = sum(price[t] * dt * ps_opf[t] for t in 1:T) + battery + terminal
J_dss = sum(price[t] * dt * ps_dss[t] for t in 1:T) + battery + terminal

@printf("OPENDSS_CHECK %s | %s T=%d | converged=%s | source %.2f pu in the DSS file, 1.05 used | devices on the substation bus skipped: %d | PV kW, DSS against instance: %.1e | line C1 %.3g nF\n",
        label, system, T, converged, source_pu_file, skipped, pv_mismatch, line_c1)
@printf("  %3s %10s %10s %8s %9s %9s %9s %9s %8s %8s %6s %6s\n", "t", "P_opf kW", "P_dss kW", "diff %", "loss_opf", "loss_dss", "fict kW", "dV max", "Vmin dss", "Vmax dss", "low", "high")
for t in 1:T
    @printf("  %3d %10.2f %10.2f %8.4f %9.3f %9.3f %9.3f %9.2e %8.4f %8.4f %6d %6d\n", t, ps_opf[t], ps_dss[t],
            100 * (ps_opf[t] - ps_dss[t]) / abs(ps_dss[t]), loss_opf[t], loss_dss[t], fict[t], dvmax[t], vmin_dss[t], vmax_dss[t], nlow[t], nhigh[t])
end
nb = (length(buses) - 1) * T
@printf("OPENDSS_SUMMARY %s | import: optimizer %.2f kWh, OpenDSS %.2f kWh (%+.4f%%) | losses: optimizer %.2f, OpenDSS %.2f kWh, of which not physical %.2f kWh\n",
        label, dt * sum(ps_opf), dt * sum(ps_dss), 100 * (sum(ps_opf) - sum(ps_dss)) / sum(ps_dss), dt * sum(loss_opf), dt * sum(loss_dss), dt * sum(fict))
@printf("OPENDSS_SUMMARY %s | voltage: largest difference %.2e pu | OpenDSS range %.4f-%.4f pu (optimizer lowest %.4f) | below limit: %d of %d bus-periods, worst %.2e pu | above limit: %d, worst %.2e pu (tolerance %.0e)\n",
        label, maximum(dvmax), minimum(vmin_dss), maximum(vmax_dss), minimum(vmin_opf), sum(nlow), nb, maximum(worstlow), sum(nhigh), maximum(worsthigh), VTOL)
@printf("OPENDSS_SUMMARY %s | objective: optimizer %.6f, with OpenDSS import %.6f (%+.4f%%)%s\n", label, J_opf, J_dss, 100 * (J_dss - J_opf) / J_opf,
        haskey(sol, :objective) ? @sprintf(" | solver-reported %.6f", sol[:objective]) : "")
