"""
Convert FlexibleParams-DistOPF 10k bus system to OpenDSS format.

This system consists of:
- Main backbone feeder (linedata.txt): 100 buses
- 101 areas (Area Data/Area1-Area101/linedata.txt): ~100 buses each
- Cross-boundary connections (CB_full.txt): connections from backbone to areas

Total: ~10,000+ buses
"""

using DelimitedFiles

# Define paths
script_dir = @__DIR__
flexparams_dir = joinpath(dirname(dirname(script_dir)), "FlexibleParams-DistOPF")
output_dir = joinpath(script_dir, "flexparams_10k_1ph")

# Create output directory
mkpath(output_dir)

println("Reading data from: $flexparams_dir")
println("Writing OpenDSS files to: $output_dir")

# System parameters (from NL_OPF_dist.m)
base_kva = 1000.0
base_kv_ll = 12.47  # kV line-to-line
base_kv_ln = base_kv_ll / sqrt(3)  # kV line-to-neutral (single phase)

# Default impedances (from MATLAB code)
r_default = 0.07  # pu
x_default = 0.01  # pu
r_boundary = 0.0001  # pu (for boundary connections)
x_boundary = 0.0001  # pu

# Load parameters (from MATLAB code)
p_load_mw = 0.1  # MW per bus
q_load_mvar = 0.01  # MVAr per bus

# DG parameters (50% penetration, every 10th bus)
p_dg_mw = 0.07  # MW
s_dg_mva = 0.07 * 1.2  # MVA

println("\n=== System Parameters ===")
println("Base: $(base_kva) kVA, $(base_kv_ll) kV L-L, $(round(base_kv_ln, digits=4)) kV L-N")
println("Default impedance: r=$(r_default) pu, x=$(x_default) pu")
println("Default load: P=$(p_load_mw) MW, Q=$(q_load_mvar) MVAr")
println("DG rating: P=$(p_dg_mw) MW, S=$(s_dg_mva) MVA")

# ===== Read main backbone feeder =====
println("\n=== Reading Main Backbone Feeder ===")
main_feeder_file = joinpath(flexparams_dir, "linedata.txt")
main_feeder = readdlm(main_feeder_file, '\t', Int)
println("Main feeder: $(size(main_feeder, 1)) lines, buses 1-$(maximum(main_feeder))")

# ===== Read cross-boundary connections =====
println("\n=== Reading Cross-Boundary Connections ===")
cb_file = joinpath(flexparams_dir, "CB_full.txt")
cb_data = readdlm(cb_file, '\t', Int)
println("Cross-boundary connections: $(size(cb_data, 1)) lines")
println("Format: [backbone_bus, area_local_bus, area_id]")

# ===== Read all area data =====
println("\n=== Reading Area Data ===")
area_data_dir = joinpath(flexparams_dir, "Area Data")
num_areas = 101

area_lines = Dict{Int, Array{Int,2}}()
for area_id in 1:num_areas
    area_dir = joinpath(area_data_dir, "Area$area_id")
    linedata_file = joinpath(area_dir, "linedata.txt")

    if isfile(linedata_file)
        lines = readdlm(linedata_file, '\t', Int)
        area_lines[area_id] = lines
    end
end

total_area_lines = sum(size(lines, 1) for lines in values(area_lines))
println("Loaded $(length(area_lines)) areas with $(total_area_lines) total lines")

# ===== Bus numbering strategy =====
# - Main feeder: 1, 2, 3, ..., 100
# - Area X, local bus Y → global bus: X*1000 + Y
#   e.g., Area 101, bus 1 → 101001
#        Area 101, bus 2 → 101002
println("\n=== Bus Numbering Strategy ===")
println("Main feeder: 1-100")
println("Area X, local bus Y → global bus X*1000 + Y")
println("Example: Area 101, bus 5 → 101005")

# ===== Generate BranchData.dss =====
println("\n=== Generating BranchData.dss ===")
line_counter = Ref(0)

open(joinpath(output_dir, "BranchData.dss"), "w") do f
    # 1. Main backbone feeder lines
    for i in 1:size(main_feeder, 1)
        fbus = main_feeder[i, 1]
        tbus = main_feeder[i, 2]
        line_counter[] += 1

        write(f, "New Line.L$(line_counter[]) Phases=1 Bus1=$fbus Bus2=$tbus ")
        write(f, "r1=$r_default x1=$x_default\n")
    end

    # 2. Cross-boundary connections (from main feeder to areas)
    for i in 1:size(cb_data, 1)
        backbone_bus = cb_data[i, 1]
        area_local_bus = cb_data[i, 2]
        area_id = cb_data[i, 3]

        # Global bus number for area bus
        area_global_bus = area_id * 1000 + area_local_bus

        line_counter[] += 1
        write(f, "New Line.L$(line_counter[]) Phases=1 Bus1=$backbone_bus Bus2=$area_global_bus ")
        write(f, "r1=$r_boundary x1=$x_boundary\n")
    end

    # 3. Lines within each area
    for area_id in sort(collect(keys(area_lines)))
        lines = area_lines[area_id]

        for i in 1:size(lines, 1)
            local_fbus = lines[i, 1]
            local_tbus = lines[i, 2]

            # Convert to global bus numbers
            global_fbus = area_id * 1000 + local_fbus
            global_tbus = area_id * 1000 + local_tbus

            line_counter[] += 1
            write(f, "New Line.L$(line_counter[]) Phases=1 Bus1=$global_fbus Bus2=$global_tbus ")
            write(f, "r1=$r_default x1=$x_default\n")
        end
    end
end

println("Generated $(line_counter[]) lines total")

# ===== Find all unique buses =====
println("\n=== Identifying All Buses ===")
all_buses = Set{Int}()

# From main feeder
for i in 1:size(main_feeder, 1)
    push!(all_buses, main_feeder[i, 1])
    push!(all_buses, main_feeder[i, 2])
end

# From cross-boundary connections
for i in 1:size(cb_data, 1)
    backbone_bus = cb_data[i, 1]
    area_local_bus = cb_data[i, 2]
    area_id = cb_data[i, 3]
    area_global_bus = area_id * 1000 + area_local_bus

    push!(all_buses, backbone_bus)
    push!(all_buses, area_global_bus)
end

# From area lines
for (area_id, lines) in area_lines
    for i in 1:size(lines, 1)
        push!(all_buses, area_id * 1000 + lines[i, 1])
        push!(all_buses, area_id * 1000 + lines[i, 2])
    end
end

all_buses_sorted = sort(collect(all_buses))
num_buses = length(all_buses_sorted)
println("Total unique buses: $(num_buses)")

# ===== Generate Loads.dss =====
println("\n=== Generating Loads.dss ===")
load_counter = Ref(0)

open(joinpath(output_dir, "Loads.dss"), "w") do f
    write(f, "! Loads based on MATLAB NL_OPF_dist.m template\n")
    write(f, "! Default: $(p_load_mw) MW, $(q_load_mvar) MVAr per bus\n\n")

    for bus in all_buses_sorted
        # Skip source bus 1
        if bus == 1
            continue
        end

        # Add load at this bus
        load_counter[] += 1
        write(f, "New load.S$(bus)P Bus=$bus.1 Phases=1 Model=1 ")
        write(f, "kV=$(round(base_kv_ln, digits=4)) ")
        write(f, "kW=$(p_load_mw * 1000) ")
        write(f, "kVAr=$(q_load_mvar * 1000) ")
        write(f, "Vminpu=0.95 Vmaxpu=1.05 Daily=LoadShapeLoadDefault\n")
    end
end

println("Generated $(load_counter[]) loads")

# ===== Generate PVSystem.dss (DG at every 10th bus in each area) =====
println("\n=== Generating PVSystem.dss ===")
pv_counter = Ref(0)

open(joinpath(output_dir, "PVSystem.dss"), "w") do f
    write(f, "! DG systems: 50% penetration at every 10th bus\n")
    write(f, "! Rating: $(p_dg_mw) MW, $(s_dg_mva) MVA\n\n")

    # DG in main feeder (every 10th bus: 2, 12, 22, ..., 92)
    for bus in 2:10:100
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

    # DG in each area (local buses 2, 12, 22, ..., 92)
    for area_id in 1:num_areas
        for local_bus in 2:10:100
            global_bus = area_id * 1000 + local_bus

            if global_bus in all_buses_sorted
                pv_counter[] += 1
                kva = s_dg_mva * 1000
                pmpp = p_dg_mw * 1000

                write(f, "New PVsystem.PV$global_bus irrad=1.0 Phases=1 Bus1=$global_bus.1 ")
                write(f, "kV=$(round(base_kv_ln, digits=4)) ")
                write(f, "kVA=$kva Pmpp=$pmpp ")
                write(f, "%cutin=0.001 %cutout=0.001 Daily=LoadShapePVDefault\n")
            end
        end
    end
end

println("Generated $(pv_counter[]) PV systems")

# ===== Generate Storage.dss (can add batteries later) =====
println("\n=== Generating Storage.dss ===")
open(joinpath(output_dir, "Storage.dss"), "w") do f
    write(f, "! Storage systems can be added here\n")
    write(f, "! Example: Add battery at bus 50\n")
    write(f, "! New Storage.Battery50 phases=1 Bus1=50 kV=$(round(base_kv_ln, digits=4)) kVA=500 kWrated=400 kWhrated=1600 %stored=62.5 %reserve=30 %EffCharge=95 %EffDischarge=95 %idlingkW=0 Vminpu=0.95 Vmaxpu=1.05 DispMode=External\n")
end

# ===== Generate Master.dss =====
println("\n=== Generating Master.dss ===")
open(joinpath(output_dir, "Master.dss"), "w") do f
    write(f, """Clear

New Circuit.flexparams_10k_1ph

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

! Let DSS estimate the voltage bases
CalcVoltageBases

""")
end

# ===== Summary =====
println("\n" * "="^60)
println("✓ Conversion Complete!")
println("="^60)
println("Output directory: $output_dir")
println("\nGenerated files:")
println("  - Master.dss")
println("  - BranchData.dss     ($(line_counter[]) lines)")
println("  - Loads.dss          ($(load_counter[]) loads)")
println("  - PVSystem.dss       ($(pv_counter[]) PV systems)")
println("  - Storage.dss        (template)")
println("\nSystem statistics:")
println("  - Total buses:       $num_buses")
println("  - Total lines:       $(line_counter[])")
println("  - Main feeder:       $(size(main_feeder, 1)) lines (buses 1-100)")
println("  - Cross-boundary:    $(size(cb_data, 1)) connections")
println("  - Area lines:        $total_area_lines lines")
println("  - Loads:             $(load_counter[])")
println("  - PV systems:        $(pv_counter[])")
println("\nTo use with tADMM:")
println("  1. Update tadmm_socp.jl: systemName = \"flexparams_10k_1ph\"")
println("  2. Set T = 6 or T = 12 for initial testing")
println("  3. Run: JULIA_NUM_THREADS=16 julia tadmm_socp.jl")
println("="^60)
