"""
Simple converter: just use main backbone linedata.txt
Ignore Area Data directories and CB decomposition labels
"""

using DelimitedFiles

# Define paths
script_dir = @__DIR__
flexparams_dir = joinpath(dirname(dirname(script_dir)), "FlexibleParams-DistOPF")
output_dir = joinpath(script_dir, "large10k_1ph")

# Create output directory
mkpath(output_dir)

println("Reading data from: $flexparams_dir")
println("Writing OpenDSS files to: $output_dir")

# System parameters
base_kva = 1000.0
base_kv_ll = 12.47  # kV line-to-line
base_kv_ln = base_kv_ll / sqrt(3)  # kV line-to-neutral

# Default impedances
r_default = 0.07  # pu
x_default = 0.01  # pu

# Load parameters
p_load_mw = 0.1  # MW per bus
q_load_mvar = 0.01  # MVAr per bus

# DG parameters (50% penetration, every 10th bus)
p_dg_mw = 0.07  # MW
s_dg_mva = 0.07 * 1.2  # MVA

println("\n=== System Parameters ===")
println("Base: $(base_kva) kVA, $(base_kv_ll) kV L-L, $(round(base_kv_ln, digits=4)) kV L-N")
println("Default impedance: r=$(r_default) pu, x=$(x_default) pu")
println("Default load: P=$(p_load_mw) MW, Q=$(q_load_mvar) MVAr")

# ===== Read main backbone feeder ONLY =====
println("\n=== Reading Main Backbone Feeder ===")
main_feeder_file = joinpath(flexparams_dir, "linedata.txt")
main_feeder = readdlm(main_feeder_file, '\t', Int)
println("Main feeder: $(size(main_feeder, 1)) lines")

# Find all unique buses
all_buses = Set{Int}()
for i in 1:size(main_feeder, 1)
    push!(all_buses, main_feeder[i, 1])
    push!(all_buses, main_feeder[i, 2])
end
all_buses_sorted = sort(collect(all_buses))
println("Total buses: $(length(all_buses_sorted)) ($(minimum(all_buses_sorted)) to $(maximum(all_buses_sorted)))")

# ===== Generate BranchData.dss =====
println("\n=== Generating BranchData.dss ===")
open(joinpath(output_dir, "BranchData.dss"), "w") do f
    for i in 1:size(main_feeder, 1)
        fbus = main_feeder[i, 1]
        tbus = main_feeder[i, 2]

        write(f, "New Line.L$i Phases=1 Bus1=$fbus Bus2=$tbus ")
        write(f, "r1=$r_default x1=$x_default\n")
    end
end
println("Generated $(size(main_feeder, 1)) lines")

# ===== Generate Loads.dss =====
println("\n=== Generating Loads.dss ===")
load_counter = Ref(0)

open(joinpath(output_dir, "Loads.dss"), "w") do f
    write(f, "! Loads: $(p_load_mw) MW, $(q_load_mvar) MVAr per bus\n\n")

    for bus in all_buses_sorted
        # Skip source bus 1
        if bus == 1
            continue
        end

        load_counter[] += 1
        write(f, "New load.S$(bus)P Bus=$bus.1 Phases=1 Model=1 ")
        write(f, "kV=$(round(base_kv_ln, digits=4)) ")
        write(f, "kW=$(p_load_mw * 1000) ")
        write(f, "kVAr=$(q_load_mvar * 1000) ")
        write(f, "Vminpu=0.95 Vmaxpu=1.05 Daily=LoadShapeLoadDefault\n")
    end
end
println("Generated $(load_counter[]) loads")

# ===== Generate PVSystem.dss =====
println("\n=== Generating PVSystem.dss ===")
pv_counter = Ref(0)

open(joinpath(output_dir, "PVSystem.dss"), "w") do f
    write(f, "! DG systems: 50% penetration at every 10th bus\n\n")

    for bus in 2:10:maximum(all_buses_sorted)
        if bus in all_buses_sorted
            pv_counter[] += 1
            kva = s_dg_mva * 1000
            pmpp = p_dg_mw * 1000

            write(f, "New PVsystem.PV$bus irrad=1.0 Phases=1 Bus1=$bus.1 ")
            write(f, "kV=$(round(base_kv_ln, digits=4)) ")
            write(f, "kVA=$kva Pmpp=$pmpp ")
            write(f, "%cutin=0.001 %cutout=0.001 Daily=LoadShapePVDefault\n")
        end
    end
end
println("Generated $(pv_counter[]) PV systems")

# ===== Generate Storage.dss =====
open(joinpath(output_dir, "Storage.dss"), "w") do f
    write(f, "! Storage systems can be added here\n")
end

# ===== Generate Master.dss =====
println("\n=== Generating Master.dss ===")
open(joinpath(output_dir, "Master.dss"), "w") do f
    write(f, """Clear

New Circuit.large10k_1ph

Edit "Vsource.source" bus1=1 pu=1.05 R0=0 X0=0.00000001 phases=1 basekv=$(round(base_kv_ln, digits=4)) R1=0 X1=0.00000001

Redirect BranchData.dss

Redirect ../LoadShapeDefault.dss
Redirect ../LoadShapePVDefault.dss
Redirect Loads.dss
Redirect PVSystem.dss
Redirect Storage.dss

Set VoltageBases = [$base_kv_ll]
Set mode = Daily
Set stepsize = 1h
Set number = 1

CalcVoltageBases

""")
end

println("\n" * "="^60)
println("✓ Conversion Complete!")
println("="^60)
println("System statistics:")
println("  - Total buses:       $(length(all_buses_sorted))")
println("  - Total lines:       $(size(main_feeder, 1))")
println("  - Loads:             $(load_counter[])")
println("  - PV systems:        $(pv_counter[])")
println("="^60)
