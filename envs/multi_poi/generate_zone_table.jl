using OpenDSSDirect

# Load the system
OpenDSSDirect.Text.Command("clear")
OpenDSSDirect.Text.Command("Redirect rawData/ieee123_5poi_1ph/Master.dss")
OpenDSSDirect.Solution.Solve()

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

# Calculate totals
global total_loads = 0.0
global total_pv = 0.0
global total_batt = 0.0

for ln in OpenDSSDirect.Loads.AllNames()
    OpenDSSDirect.Loads.Name(ln)
    global total_loads += OpenDSSDirect.Loads.kW()
end

for pv in OpenDSSDirect.PVsystems.AllNames()
    OpenDSSDirect.PVsystems.Name(pv)
    global total_pv += OpenDSSDirect.PVsystems.Pmpp()
end

for b in OpenDSSDirect.Storages.AllNames()
    OpenDSSDirect.Storages.Name(b)
    global total_batt += OpenDSSDirect.Storages.kW()
end

num_loads = length(OpenDSSDirect.Loads.AllNames())

println("\n=== System Totals ===")
println("Total Load Buses: $num_loads")
println("Total Load: $(round(total_loads, digits=1)) kW")
println("Total PV: $(round(total_pv, digits=1)) kW")
println("Total Batteries: $(round(total_batt, digits=1)) kW")

# Get PV and Battery capacities by zone
zone_pv_kw = Dict{Int, Float64}()
zone_batt_kw = Dict{Int, Float64}()

for zone in 1:5
    local pv_kw = 0.0
    for bus in pv_zones[zone]
        pv_name = "pv$bus"
        try
            OpenDSSDirect.PVsystems.Name(pv_name)
            pv_kw += OpenDSSDirect.PVsystems.Pmpp()
        catch
        end
    end
    zone_pv_kw[zone] = pv_kw
    
    local batt_kw = 0.0
    for bus in batt_zones[zone]
        batt_name = "batt$bus"
        try
            OpenDSSDirect.Storages.Name(batt_name)
            batt_kw += OpenDSSDirect.Storages.kW()
        catch
        end
    end
    zone_batt_kw[zone] = batt_kw
end

println("\n=== Zone Resource Distribution ===")
println("Zone | Subs | PV (kW) | Batteries (kW)")
println("-" ^ 50)
for zone in 1:5
    pv_count = length(pv_zones[zone])
    batt_count = length(batt_zones[zone])
    println("Zone $zone | $(zone_subs[zone]) | $(round(zone_pv_kw[zone], digits=1)) kW ($pv_count units) | $(round(zone_batt_kw[zone], digits=1)) kW ($batt_count units)")
end

# Generate LaTeX table
latex_output = """
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
    latex_output *= "Zone $zone & $(zone_subs[zone]) & $pv_kw & $batt_count & $batt_kw \\\\\n"
end

latex_output *= """\\midrule
Total & -- & $(round(Int, total_pv)) & $(length(vcat(values(batt_zones)...))) & $(round(Int, total_batt)) \\\\
\\bottomrule
\\end{tabular}
\\end{table}
"""

# Save LaTeX table
open("envs/multi_poi/zone_resource_table.tex", "w") do f
    write(f, latex_output)
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
  - Total Load Buses: $num_loads
  - Total Load: $(round(total_loads, digits=1)) kW
  - Total PV: $(round(total_pv, digits=1)) kW ($(sum(length(pv_zones[z]) for z in 1:5)) units)
  - Total Batteries: $(round(total_batt, digits=1)) kW ($(sum(length(batt_zones[z]) for z in 1:5)) units)
"""

open("envs/multi_poi/zone_resource_summary.txt", "w") do f
    write(f, txt_output)
end

println("✓ Text summary saved to: envs/multi_poi/zone_resource_summary.txt")
println("\nDone!")
