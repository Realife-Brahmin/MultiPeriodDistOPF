# Zone mappings
zone_subs = Dict(1 => "1s", 2 => "2s", 3 => "3s", 4 => "4s", 5 => "5s")
zone_names = Dict(1 => "SW-Residential", 2 => "NW-Commercial", 3 => "N-Industrial", 4 => "NE-Office", 5 => "SE-Charging")

pv_zones = Dict(
    1 => [3, 9, 18, 26, 35],
    2 => [41, 48, 53],
    3 => [61, 67, 73, 79],
    4 => [87, 94, 101, 108, 116],
    5 => []
)

batt_zones = Dict(
    1 => [3, 7, 11, 18, 22, 30, 34, 37],
    2 => [41, 47, 50, 53],
    3 => [58, 62, 67, 71, 75, 79],
    4 => [84, 87, 92, 97, 101, 106, 111, 116],
    5 => []
)

# Parse PV file
pv_capacities = Dict{Int, Float64}()
for line in eachline("rawData/ieee123_5poi_1ph/PVSystem.dss")
    if occursin("New PVsystem", line)
        # Extract bus number and Pmpp
        m_bus = match(r"Bus1=(\d+)", line)
        m_pmpp = match(r"Pmpp=([\d.]+)", line)
        if !isnothing(m_bus) && !isnothing(m_pmpp)
            bus = parse(Int, m_bus.captures[1])
            pmpp = parse(Float64, m_pmpp.captures[1])
            pv_capacities[bus] = pmpp
        end
    end
end

# Parse Storage file
batt_capacities = Dict{Int, Float64}()
for line in eachline("rawData/ieee123_5poi_1ph/Storage.dss")
    if occursin("New Storage", line)
        # Extract bus number and kWrated
        m_bus = match(r"Bus1=(\d+)", line)
        m_kw = match(r"kWrated=([\d.]+)", line)
        if !isnothing(m_bus) && !isnothing(m_kw)
            bus = parse(Int, m_bus.captures[1])
            kw = parse(Float64, m_kw.captures[1])
            batt_capacities[bus] = kw
        end
    end
end

# Calculate zone totals
zone_pv_kw = Dict{Int, Float64}()
zone_batt_kw = Dict{Int, Float64}()

for zone in 1:5
    pv_total = sum(get(pv_capacities, bus, 0.0) for bus in pv_zones[zone]; init=0.0)
    zone_pv_kw[zone] = pv_total
    
    batt_total = sum(get(batt_capacities, bus, 0.0) for bus in batt_zones[zone]; init=0.0)
    zone_batt_kw[zone] = batt_total
end

# System totals
total_pv = sum(values(pv_capacities))
total_batt = sum(values(batt_capacities))
total_pv_units = length(pv_capacities)
total_batt_units = length(batt_capacities)

println("\n=== System Totals ===")
println("Total PV: $(round(total_pv, digits=1)) kW ($total_pv_units units)")
println("Total Batteries: $(round(total_batt, digits=1)) kW ($total_batt_units units)")

println("\n=== Zone Resource Distribution ===")
println("Zone | Subs | PV (kW) | Batteries (kW)")
println("-" ^ 50)
for zone in 1:5
    pv_count = length(pv_zones[zone])
    batt_count = length(batt_zones[zone])
    println("Zone $zone | $(zone_subs[zone]) | $(round(zone_pv_kw[zone], digits=1)) kW ($pv_count units) | $(round(zone_batt_kw[zone], digits=1)) kW ($batt_count units)")
end

# Generate LaTeX table
function generate_latex_table()
    latex_str = """
\\begin{table}[t]
\\centering
\\caption{Resource Distribution Across Network Zones}
\\label{tab:zone_resources}
\\begin{tabular}{lcccc}
\\toprule
Zone & Substation & PV Capacity & Battery Units & Battery Power \\\\
     &            & (kW)        &               & (kW) \\\\
\\midrule
"""
    
    for zone in 1:5
        pv_count = length(pv_zones[zone])
        batt_count = length(batt_zones[zone])
        pv_kw = round(Int, zone_pv_kw[zone])
        batt_kw = round(Int, zone_batt_kw[zone])
        latex_str *= "Zone $zone & $(zone_subs[zone]) & $pv_kw & $batt_count & $batt_kw \\\\\n"
    end
    
    latex_str *= """\\midrule
Total & -- & $(round(Int, total_pv)) & $(total_batt_units) & $(round(Int, total_batt)) \\\\
\\bottomrule
\\end{tabular}
\\end{table}
"""
    return latex_str
end

latex_table = generate_latex_table()

# Save LaTeX table
mkpath("envs/multi_poi")
open("envs/multi_poi/zone_resource_table.tex", "w") do f
    write(f, latex_table)
end

println("\n✓ LaTeX table saved to: envs/multi_poi/zone_resource_table.tex")

# Generate text summary
txt_output = """
IEEE 123-Bus 5-POI System - Zone Resource Summary
==================================================

Zone 1 (Substation 1s - SW-Residential):
  - PV: $(round(zone_pv_kw[1], digits=1)) kW ($(length(pv_zones[1])) units at buses: $(join(pv_zones[1], ", ")))
  - Batteries: $(round(zone_batt_kw[1], digits=1)) kW ($(length(batt_zones[1])) units at buses: $(join(batt_zones[1], ", ")))

Zone 2 (Substation 2s - NW-Commercial):
  - PV: $(round(zone_pv_kw[2], digits=1)) kW ($(length(pv_zones[2])) units at buses: $(join(pv_zones[2], ", ")))
  - Batteries: $(round(zone_batt_kw[2], digits=1)) kW ($(length(batt_zones[2])) units at buses: $(join(batt_zones[2], ", ")))

Zone 3 (Substation 3s - N-Industrial):
  - PV: $(round(zone_pv_kw[3], digits=1)) kW ($(length(pv_zones[3])) units at buses: $(join(pv_zones[3], ", ")))
  - Batteries: $(round(zone_batt_kw[3], digits=1)) kW ($(length(batt_zones[3])) units at buses: $(join(batt_zones[3], ", ")))

Zone 4 (Substation 4s - NE-Office):
  - PV: $(round(zone_pv_kw[4], digits=1)) kW ($(length(pv_zones[4])) units at buses: $(join(pv_zones[4], ", ")))
  - Batteries: $(round(zone_batt_kw[4], digits=1)) kW ($(length(batt_zones[4])) units at buses: $(join(batt_zones[4], ", ")))

Zone 5 (Substation 5s - SE-Charging):
  - PV: $(round(zone_pv_kw[5], digits=1)) kW ($(length(pv_zones[5])) units)
  - Batteries: $(round(zone_batt_kw[5], digits=1)) kW ($(length(batt_zones[5])) units)

System Totals:
  - Total PV: $(round(total_pv, digits=1)) kW ($(total_pv_units) units)
  - Total Batteries: $(round(total_batt, digits=1)) kW ($(total_batt_units) units)
"""

open("envs/multi_poi/zone_resource_summary.txt", "w") do f
    write(f, txt_output)
end

println("✓ Text summary saved to: envs/multi_poi/zone_resource_summary.txt")
println("\nDone!")
