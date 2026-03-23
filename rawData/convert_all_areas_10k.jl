"""
Convert FlexibleParams 10k system: main backbone + ALL 101 area subdirectories
"""

using DelimitedFiles

# Paths
script_dir = @__DIR__
flexparams_dir = joinpath(dirname(dirname(script_dir)), "FlexibleParams-DistOPF")
output_dir = joinpath(script_dir, "large10k_1ph")

mkpath(output_dir)

println("Reading from: $flexparams_dir")
println("Writing to: $output_dir\n")

# System parameters
base_kva = 1000.0
base_kv_ll = 12.47
base_kv_ln = base_kv_ll / sqrt(3)
r_default = 0.07  # pu
x_default = 0.01  # pu
r_boundary = 0.0001  # pu (CB connections)
x_boundary = 0.0001  # pu
p_load_mw = 0.1
q_load_mvar = 0.01
p_dg_mw = 0.07
s_dg_mva = 0.084

# Read main backbone
println("=== Main Backbone ===")
main_feeder = readdlm(joinpath(flexparams_dir, "linedata.txt"), '\t', Int)
println("Loaded: $(size(main_feeder, 1)) lines")

# Read ALL area data
println("\n=== Loading All Areas ===")
area_data_dir = joinpath(flexparams_dir, "Area Data")
area_lines = Dict{Int, Array{Int,2}}()

for area_id in 1:101
    linedata_file = joinpath(area_data_dir, "Area$area_id", "linedata.txt")
    if isfile(linedata_file)
        lines = readdlm(linedata_file, '\t', Int)
        area_lines[area_id] = lines
        if area_id <= 5 || area_id == 101
            println("  Area$area_id: $(size(lines, 1)) lines")
        end
    end
end
total_area_lines = sum(size(lines, 1) for lines in values(area_lines))
println("  ... (showing first 5 and last)")
println("Total: $(length(area_lines)) areas, $total_area_lines lines")

# Read cross-boundary connections
println("\n=== Cross-Boundary Connections ===")
cb_data = readdlm(joinpath(flexparams_dir, "CB_full.txt"), '\t', Int)
println("Loaded: $(size(cb_data, 1)) CB connections")

# Bus numbering: Area X, local bus Y → global bus X*1000 + Y
println("\n=== Generating BranchData.dss ===")
line_counter = Ref(0)

open(joinpath(output_dir, "BranchData.dss"), "w") do f
    # 1. Main backbone
    for i in 1:size(main_feeder, 1)
        fbus = main_feeder[i, 1]
        tbus = main_feeder[i, 2]
        line_counter[] += 1
        write(f, "New Line.L$(line_counter[]) Phases=1 Bus1=$fbus Bus2=$tbus r1=$r_default x1=$x_default\n")
    end

    # 2. Cross-boundary connections
    for i in 1:size(cb_data, 1)
        backbone_bus = cb_data[i, 1]
        area_local_bus = cb_data[i, 2]
        cb_area_id = cb_data[i, 3]  # This is 101-120 in CB_full.txt

        # Map CB area_id to actual Area directory number
        # CB_full.txt uses area_id 101-120
        # Area directories are Area1-Area101
        # Mapping: CB area_id 101 → Area1, 102 → Area2, etc.
        area_dir_num = cb_area_id - 100  # 101 → 1, 102 → 2, ..., 120 → 20

        area_global_bus = area_dir_num * 1000 + area_local_bus

        line_counter[] += 1
        write(f, "New Line.L$(line_counter[]) Phases=1 Bus1=$backbone_bus Bus2=$area_global_bus r1=$r_boundary x1=$x_boundary\n")
    end

    # 3. Internal lines within each area
    for area_id in sort(collect(keys(area_lines)))
        lines = area_lines[area_id]
        for i in 1:size(lines, 1)
            local_fbus = lines[i, 1]
            local_tbus = lines[i, 2]
            global_fbus = area_id * 1000 + local_fbus
            global_tbus = area_id * 1000 + local_tbus

            line_counter[] += 1
            write(f, "New Line.L$(line_counter[]) Phases=1 Bus1=$global_fbus Bus2=$global_tbus r1=$r_default x1=$x_default\n")
        end
    end
end
println("Generated $(line_counter[]) total lines")

# Find all unique buses
println("\n=== Identifying All Buses ===")
all_buses = Set{Int}()

for i in 1:size(main_feeder, 1)
    push!(all_buses, main_feeder[i, 1])
    push!(all_buses, main_feeder[i, 2])
end

for i in 1:size(cb_data, 1)
    push!(all_buses, cb_data[i, 1])
    push!(all_buses, cb_data[i, 3] * 1000 + cb_data[i, 2])
end

for (area_id, lines) in area_lines
    for i in 1:size(lines, 1)
        push!(all_buses, area_id * 1000 + lines[i, 1])
        push!(all_buses, area_id * 1000 + lines[i, 2])
    end
end

all_buses_sorted = sort(collect(all_buses))
println("Total unique buses: $(length(all_buses_sorted))")

# Generate Loads
println("\n=== Generating Loads.dss ===")
load_counter = Ref(0)

open(joinpath(output_dir, "Loads.dss"), "w") do f
    for bus in all_buses_sorted
        if bus == 1  # Skip source
            continue
        end

        load_counter[] += 1
        write(f, "New load.S$(bus)P Bus=$bus.1 Phases=1 Model=1 ")
        write(f, "kV=$(round(base_kv_ln, digits=4)) ")
        write(f, "kW=$(p_load_mw * 1000) kVAr=$(q_load_mvar * 1000) ")
        write(f, "Vminpu=0.95 Vmaxpu=1.05 Daily=LoadShapeLoadDefault\n")
    end
end
println("Generated $(load_counter[]) loads")

# Generate PV (every 10th bus)
println("\n=== Generating PVSystem.dss ===")
pv_counter = Ref(0)

open(joinpath(output_dir, "PVSystem.dss"), "w") do f
    # Main feeder PV
    for bus in 2:10:100
        if bus in all_buses_sorted
            pv_counter[] += 1
            write(f, "New PVsystem.PV$bus irrad=1.0 Phases=1 Bus1=$bus.1 ")
            write(f, "kV=$(round(base_kv_ln, digits=4)) kVA=$(s_dg_mva * 1000) Pmpp=$(p_dg_mw * 1000) ")
            write(f, "%cutin=0.001 %cutout=0.001 Daily=LoadShapePVDefault\n")
        end
    end

    # Area PV (every 10th bus in each area)
    for area_id in 1:101
        for local_bus in 2:10:120
            global_bus = area_id * 1000 + local_bus
            if global_bus in all_buses_sorted
                pv_counter[] += 1
                write(f, "New PVsystem.PV$global_bus irrad=1.0 Phases=1 Bus1=$global_bus.1 ")
                write(f, "kV=$(round(base_kv_ln, digits=4)) kVA=$(s_dg_mva * 1000) Pmpp=$(p_dg_mw * 1000) ")
                write(f, "%cutin=0.001 %cutout=0.001 Daily=LoadShapePVDefault\n")
            end
        end
    end
end
println("Generated $(pv_counter[]) PV systems")

# Storage template
open(joinpath(output_dir, "Storage.dss"), "w") do f
    write(f, "! Storage can be added here\n")
end

# Master file
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

println("\n" * "="^70)
println("✓ CONVERSION COMPLETE")
println("="^70)
println("Buses:     $(length(all_buses_sorted))")
println("Lines:     $(line_counter[])")
println("  - Backbone:     $(size(main_feeder, 1))")
println("  - Cross-bound:  $(size(cb_data, 1))")
println("  - Area internal: $total_area_lines")
println("Loads:     $(load_counter[])")
println("PV:        $(pv_counter[])")
println("="^70)
