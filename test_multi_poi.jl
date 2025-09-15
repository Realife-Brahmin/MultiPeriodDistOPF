# multi_poi.jl
# Onetime script - uncomment if needed, comment once finished
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "envs", "multi_poi"))
# Pkg.add("Crayons")
# Pkg.add("JuMP")
# Pkg.add("Ipopt")
# Pkg.add("OpenDSSDirect")
# Pkg.instantiate()
# Pkg.precompile()

# ----------------------------------------------------------------------

data = Dict()
V_1_pu = 1.07
delta_1_deg = 0.0
V_2_pu = 1.07
delta_2_deg = 0.0
alpha_share = 0.75

r1_ohm = 0.025; r2_ohm=0.025; x1_ohm=0.005; x2_ohm=0.005;
P_L_kW = 50000
Q_L_kW = P_L_kW*0.75
kVA_B = 1000
kV_B = 11.547 # 20kV line-to-line, 11.547kV line-to-neutral
C_1_dollar_per_kWh = 0.50
C_2_dollar_per_kWh = 0.50

data = Dict(
    :V_1_pu => V_1_pu,
    :delta_1_deg => delta_1_deg,
    :V_2_pu => V_2_pu,
    :delta_2_deg => delta_2_deg,
    :alpha_share => alpha_share,
    :r1_ohm => r1_ohm,
    :x1_ohm => x1_ohm,
    :r2_ohm => r2_ohm,
    :x2_ohm => x2_ohm,
    :P_L_kW => P_L_kW,
    :Q_L_kW => Q_L_kW,
    :kVA_B => kVA_B,
    :kV_B => kV_B,
    :C_1_dollar_per_kWh => C_1_dollar_per_kWh,
    :C_2_dollar_per_kWh => C_2_dollar_per_kWh
)

function process_data!(data)
    # Base values
    S_base = data[:kVA_B]
    V_base = data[:kV_B]

    # Per-unit conversions
    P_L_pu = data[:P_L_kW] / S_base
    Q_L_pu = data[:Q_L_kW] / S_base
    r1_pu = data[:r1_ohm] / ((V_base^2) / S_base)
    x1_pu = data[:x1_ohm] / ((V_base^2) / S_base)
    r2_pu = data[:r2_ohm] / ((V_base^2) / S_base)
    x2_pu = data[:x2_ohm] / ((V_base^2) / S_base)

    # Cost per unit energy (pu)
    C_1_dollar_pu = data[:C_1_dollar_per_kWh] / S_base
    C_2_dollar_pu = data[:C_2_dollar_per_kWh] / S_base

    # Angles in radians
    delta_1_rad = data[:delta_1_deg] * pi / 180
    delta_2_rad = data[:delta_2_deg] * pi / 180

    # Voltage limits
    Vminpu = 0.90
    Vmaxpu = 1.10

    # Update data dict with all required fields
    data[:P_L_pu] = P_L_pu
    data[:Q_L_pu] = Q_L_pu
    data[:r1_pu] = r1_pu
    data[:x1_pu] = x1_pu
    data[:r2_pu] = r2_pu
    data[:x2_pu] = x2_pu
    data[:C_1_dollar_pu] = C_1_dollar_pu
    data[:C_2_dollar_pu] = C_2_dollar_pu
    data[:delta_1_rad] = delta_1_rad
    data[:delta_2_rad] = delta_2_rad
    data[:Vminpu] = Vminpu
    data[:Vmaxpu] = Vmaxpu

end

process_data!(data)

using JuMP
using Ipopt

function solve_two_poi_opf(data)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    # Variables
    @variable(model, P_1j >= 0)
    @variable(model, Q_1j)
    @variable(model, l_1j >= 0)
    @variable(model, P_2j >= 0)
    @variable(model, Q_2j)
    @variable(model, l_2j >= 0)
    @variable(model, v_j >= data[:Vminpu]^2)
    @constraint(model, v_j <= data[:Vmaxpu]^2)

    # Extract parameters
    r1j = data[:r1_pu]
    x1j = data[:x1_pu]
    r2j = data[:r2_pu]
    x2j = data[:x2_pu]
    v_1 = data[:V_1_pu]^2
    v_2 = data[:V_2_pu]^2
    alpha = data[:alpha_share]
    P_L_j = data[:P_L_pu]
    Q_L_j = data[:Q_L_pu]
    C_1 = data[:C_1_dollar_pu]
    C_2 = data[:C_2_dollar_pu]

    # Constraints
    # 1. NRealPB
    @constraint(model, (P_1j - r1j * l_1j) + (P_2j - r2j * l_2j) == P_L_j)
    # 2. NReacPB
    @constraint(model, (Q_1j - x1j * l_1j) + (Q_2j - x2j * l_2j) == Q_L_j)
    # 3. KVL1
    @constraint(model, v_j == v_1 - 2 * (r1j * P_1j + x1j * Q_1j) + (r1j^2 + x1j^2) * l_1j)
    # 4. KVL2
    @constraint(model, v_j == v_2 - 2 * (r2j * P_2j + x2j * Q_2j) + (r2j^2 + x2j^2) * l_2j)
    # 5. BCPF1
    @constraint(model, P_1j^2 + Q_1j^2 >= v_1 * l_1j)
    # 6. BCPF2
    @constraint(model, P_2j^2 + Q_2j^2 >= v_2 * l_2j)
    # 7. Power sharing
    @constraint(model, P_2j == alpha * P_1j)
    # 8. Voltage limits (already in variable bounds)
    # 9. P_1j, P_2j >= 0 (already in variable bounds)

    # Objective
    @objective(model, Min, C_1 * P_1j + C_2 * P_2j)

    # Solve
    optimize!(model)

    # Collect results
    modelDict = Dict(
        :P_1j => value(P_1j),
        :Q_1j => value(Q_1j),
        :l_1j => value(l_1j),
        :P_2j => value(P_2j),
        :Q_2j => value(Q_2j),
        :l_2j => value(l_2j),
        :v_j => value(v_j)
    )
    return modelDict
end

# Example usage:
modelDict = solve_two_poi_opf(data)

import OpenDSSDirect as dss
using Printf

systemName = "small2poi_1ph"
rawDataFolder = "rawData/"
systemFile = rawDataFolder * systemName * "/Master.dss"

dss.Text.Command("Clear")
dss.Text.Command("Redirect " * systemFile)

deltas = [-10; -5; collect(-2:0.5:2); 5; 10]
deltas = vcat(
    [-10, -5],
    collect(-2:0.5:-0.5),
    collect(-0.4:0.1:0.4),
    collect(0.5:0.5:2),
    [5, 10]
)
PSubs1 = []
QSubs1 = []
PSubs2 = []
QSubs2 = []

for delta in deltas
    dss.Text.Command("Edit Vsource.grid2 angle=" * string(-delta))
    dss.Text.Command("Solve")

    dss.Circuit.SetActiveElement("Vsource.grid1")
    powers1 = dss.CktElement.Powers()
    push!(PSubs1, -real(powers1[1]))
    push!(QSubs1, -imag(powers1[1]))

    dss.Circuit.SetActiveElement("Vsource.grid2")
    powers2 = dss.CktElement.Powers()
    push!(PSubs2, -real(powers2[1]))
    push!(QSubs2, -imag(powers2[1]))
end




# Generalized retrieval of substation voltages (pu) and BasekV
Vsource_names = sort(setdiff(dss.Vsources.AllNames(), ["source"]))
Vsource_pus = Dict{String, Float64}()
Vsource_basekv = Dict{String, Float64}()
for vname in Vsource_names
    dss.Vsources.Name(vname)
    Vsource_pus[vname] = dss.Vsources.PU()
    Vsource_basekv[vname] = dss.Vsources.BasekV()
end

# Generalized retrieval of total loading (sum of all loads)
load_names = sort(dss.Loads.AllNames())
global total_P = 0.0
global total_Q = 0.0
for lname in load_names
    dss.Loads.Name(lname)
    global total_P += dss.Loads.kW()
    global total_Q += dss.Loads.kvar()
end
P_L_kW = total_P
Q_L_kVAr = total_Q

header = @sprintf("%-6s | %-12s | %-12s | %-12s | %-12s", "delta", "PSubs1 [kW]", "PSubs2 [kW]", "QSubs1 [kVAr]", "QSubs2 [kVAr]")
separator = "-"^length(header)
println(separator)

# Print to console
println("Total loading: $(round(P_L_kW, digits=2)) kW, $(round(Q_L_kVAr, digits=2)) kVAr")

println("Substation voltages (pu): " *
    join(["$(name)=$( @sprintf("%.2f", Vsource_pus[name]) )" for name in Vsource_names], ", "))
println("Substation BasekV: " *
    join(["$(name)=$( @sprintf("%.4f", Vsource_basekv[name]) )" for name in Vsource_names], ", "))

println(separator)

using Crayons
println(header)
println(separator)
for i in eachindex(deltas)
    p1 = PSubs1[i]
    p2 = PSubs2[i]
    p1_str = p1 > 0 ? Crayon(foreground = :green)(@sprintf("%-12.2f", p1)) : Crayon(foreground = :red)(@sprintf("%-12.2f", p1))
    p2_str = p2 > 0 ? Crayon(foreground = :green)(@sprintf("%-12.2f", p2)) : Crayon(foreground = :red)(@sprintf("%-12.2f", p2))
    local delta_str = (p1 < 0 || p2 < 0) ? Crayon(foreground = :red)(@sprintf("%-6.1f", deltas[i])) : @sprintf("%-6.1f", deltas[i])
    @printf("%s | %s | %s | %-12.2f | %-12.2f\n",
        delta_str,
        p1_str,
        p2_str,
        QSubs1[i],
        QSubs2[i]
    )
end
# Add OpenDSS optimal row
opt_p1 = modelDict[:P_1j] * kVA_B
opt_p2 = modelDict[:P_2j] * kVA_B
opt_q1 = modelDict[:Q_1j] * kVA_B
opt_q2 = modelDict[:Q_2j] * kVA_B
opt_delta_str = Crayon(foreground = :yellow)(@sprintf("%-6s", "'Optimal'"))
opt_p1_str = Crayon(foreground = :yellow)(@sprintf("%-12.2f", opt_p1))
opt_p2_str = Crayon(foreground = :yellow)(@sprintf("%-12.2f", opt_p2))
@printf("%s | %s | %s | %-12.2f | %-12.2f\n",
    opt_delta_str,
    opt_p1_str,
    opt_p2_str,
    opt_q1,
    opt_q2
)
println(separator)

# Write to txt file
using Printf
import Dates
output_dir = "processedData/" * systemName * "/"
isdir(output_dir) || mkpath(output_dir)
output_file = output_dir * "substation-power-distribution-load_$(Int(round(P_L_kW)))_kW.txt"
open(output_file, "w") do io
    println(io, "Total loading: $(round(P_L_kW, digits=2)) kW, $(round(Q_L_kVAr, digits=2)) kVAr")
    println(io, "Substation voltages (pu): " *
        join(["$(name)=$( @sprintf("%.2f", Vsource_pus[name]) )" for name in Vsource_names], ", "))
    println(io, "Substation BasekV: " *
        join(["$(name)=$( @sprintf("%.4f", Vsource_basekv[name]) )" for name in Vsource_names], ", "))
    println(io, header)
    println(io, separator)
    for i in eachindex(deltas)
        @printf(io, "%-6.1f | %-12.2f | %-12.2f | %-12.2f | %-12.2f\n",
            deltas[i],
            PSubs1[i],
            PSubs2[i],
            QSubs1[i],
            QSubs2[i]
        )
    end
    # Add OpenDSS optimal row
    @printf(io, "%-6s | %-12.2f | %-12.2f | %-12.2f | %-12.2f\n",
        "'Optimal'",
        opt_p1,
        opt_p2,
        opt_q1,
        opt_q2
    )
end

    # --- Additional simulation: gen2 as fixed generator, Vsource.grid2 disabled ---

    # Set up OpenDSS for gen2 simulation
    dss.Text.Command("Edit Vsource.grid2 enabled=no")
    dss.Text.Command("Edit Generator.gen2 enabled=yes")
    dss.Text.Command(@sprintf("Edit Generator.gen2 enabled=yes kW=%.6f kvar=%.6f", modelDict[:P_2j]*kVA_B, modelDict[:Q_2j]*kVA_B))
    # dss.Text.Command("Edit Generator.gen2 enabled=yes kW=$(modelDict[:P_2j] * kVA_B)")
    dss.Text.Command("Solve")

    # Retrieve delta at bus 2s (angle of bus 2s w.r.t. reference)
    println("\n--- Simulation with gen2 as fixed generator (P,Q from OPF) ---")
    println(@sprintf("gen2 set to: P = %.2f kW, Q = %.2f kVAr", modelDict[:P_2j]*kVA_B, modelDict[:Q_2j]*kVA_B))
    # println(@sprintf("gen2 set to: P = %.2f kW", modelDict[:P_2j] * kVA_B))


    # --- Retrieve angle at bus 1s and 2s using CktElement.VoltagesMagAng() ---
    dss.Vsources.Name("grid1")
    dss.Circuit.SetActiveElement("Vsource.grid1")
    vmang_grid1 = dss.CktElement.VoltagesMagAng()
    angle_1s = vmang_grid1[2,1]  # 2nd row, 1st col: angle at terminal 1 of grid1

    dss.Generators.Name("gen2")
    dss.Circuit.SetActiveElement("Generator.gen2")
    vmang_gen2 = dss.CktElement.VoltagesMagAng()
    angle_2s = vmang_gen2[2,1]   # 2nd row, 1st col: angle at terminal 1 of gen2

    delta_12 = angle_1s - angle_2s

    delta_crayon = Crayon(foreground = :blue, bold = true)
    delta_str = delta_crayon(@sprintf("%.4f deg", delta_12))
    println(delta_str)

    # Also write to txt file
    open(output_file, "a") do io
        println(io, "$(round(delta_12, digits=2)) deg")
    end