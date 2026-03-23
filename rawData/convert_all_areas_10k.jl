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

# Read cross-boundary connections first
println("\n=== Cross-Boundary Connections ===")
cb_data = readdlm(joinpath(flexparams_dir, "CB_full.txt"), '\t', Int)
println("Loaded: $(size(cb_data, 1)) CB connections")

# Determine which areas are actually connected via CB
println("\n=== Identifying Connected Areas from CB ===")
cb_area_ids = unique(cb_data[:, 3])
connected_area_dirs = sort([cb_id - 100 for cb_id in cb_area_ids])  # 101→1, 102→2, etc.
println("CB mentions $(length(cb_area_ids)) unique area IDs: $(sort(cb_area_ids))")
println("Mapping to area directories: $(connected_area_dirs)")

# Read ONLY the connected area data
println("\n=== Loading Connected Areas ===")
area_data_dir = joinpath(flexparams_dir, "Area Data")
area_lines = Dict{Int, Array{Int,2}}()

for area_id in connected_area_dirs
    linedata_file = joinpath(area_data_dir, "Area$area_id", "linedata.txt")
    if isfile(linedata_file)
        lines = readdlm(linedata_file, '\t', Int)
        area_lines[area_id] = lines
        if area_id <= 5 || area_id == maximum(connected_area_dirs)
            println("  Area$area_id: $(size(lines, 1)) lines")
        end
    else
        println("  WARNING: Area$area_id not found!")
    end
end
total_area_lines = sum(size(lines, 1) for lines in values(area_lines))
println("  ... (showing first 5 and last)")
println("Total: $(length(area_lines)) areas, $total_area_lines lines")

# Bus numbering: 6-digit format AAABBB (Area X, local bus Y → X*1000 + Y)

# Build cross-boundary area mapping (before file writing)
println("\n=== Building CB Area Mapping ===")
cb_areas_seen = Set{Int}()
cb_area_to_backbone = Dict{Int, Int}()

for i in 1:size(cb_data, 1)
    backbone_bus = cb_data[i, 1]
    cb_area_id = cb_data[i, 3]

    # Only create ONE connection per area (to avoid mesh)
    if cb_area_id ∉ cb_areas_seen
        area_dir_num = cb_area_id - 100  # 101 → 1, 102 → 2, etc.
        cb_area_to_backbone[area_dir_num] = backbone_bus
        push!(cb_areas_seen, cb_area_id)
    end
end
println("Mapped $(length(cb_area_to_backbone)) areas to backbone connections")

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

    # 2. Cross-boundary connections to area ROOT buses
    for (area_dir_num, backbone_bus) in sort(collect(cb_area_to_backbone))
        area_root_bus = area_dir_num * 1000 + 1  # Area root is local bus 1

        line_counter[] += 1
        write(f, "New Line.L$(line_counter[]) Phases=1 Bus1=$backbone_bus Bus2=$area_root_bus r1=$r_boundary x1=$x_boundary\n")
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

# Backbone buses
for i in 1:size(main_feeder, 1)
    push!(all_buses, main_feeder[i, 1])
    push!(all_buses, main_feeder[i, 2])
end

# Area root buses (connected via CB)
for (area_dir_num, backbone_bus) in cb_area_to_backbone
    push!(all_buses, backbone_bus)  # Backbone connection point
    push!(all_buses, area_dir_num * 1000 + 1)  # Area root (local bus 1)
end

# Area internal buses
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

    # Area PV (every 10th bus in connected areas only)
    for area_id in connected_area_dirs
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
println("  - Bus numbering: 6-digit AAABBB format")
println("  - Backbone: 1-$(maximum(main_feeder))")
println("  - Areas: $(minimum(connected_area_dirs))001-$(maximum(connected_area_dirs))XXX")
println("Lines:     $(line_counter[])")
println("  - Backbone:     $(size(main_feeder, 1))")
println("  - Cross-bound:  $(length(cb_area_to_backbone))")
println("  - Area internal: $total_area_lines")
println("Loads:     $(load_counter[])")
println("PV:        $(pv_counter[])")
println("="^70)
