# ZonalPlotter.jl - Zonal aggregate plotting for Multi-POI MPOPF
# Plots PV and battery actions aggregated by zone (substation)

module ZonalPlotter

using Plots
using Printf
using LaTeXStrings
using Statistics

export get_zone_assignments, plot_zonal_pv_dispatch, plot_zonal_battery_actions, print_zone_statistics

"""
    get_zone_assignments(data)

Build dictionaries mapping zones to their PV and battery buses.

# Returns
- Dict with keys :pv_by_zone, :batt_by_zone
"""
function get_zone_assignments(data)
    # PV buses
    pv_by_zone = Dict{String, Vector{Int}}(
        "1s" => [3, 9, 18, 26, 35],
        "2s" => [41, 48, 53],
        "3s" => [61, 67, 73, 79],
        "4s" => [87, 94, 101, 108, 116],
        "5s" => Int[]
    )
    
    # Battery buses
    batt_by_zone = Dict{String, Vector{Int}}(
        "1s" => [3, 7, 11, 18, 22, 30, 34, 37],
        "2s" => [41, 47, 50, 53],
        "3s" => [58, 62, 67, 71, 75, 79],
        "4s" => [84, 87, 92, 97, 101, 106, 111, 116],
        "5s" => Int[]
    )
    
    return Dict(:pv_by_zone => pv_by_zone, :batt_by_zone => batt_by_zone)
end

"""
    plot_zonal_pv_dispatch(result, data, slack_sub; showPlots=true, savePlots=false, filename="zonal_pv_dispatch.png")

Plot aggregated PV real power dispatch for each zone.
Uses fixed p_D from solar profiles (not optimized).

# Arguments
- result: OPF solution result
- data: System data dictionary
- slack_sub: Slack substation identifier
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot
- filename: Filename for saving
"""
function plot_zonal_pv_dispatch(result, data, slack_sub; showPlots::Bool=true, savePlots::Bool=false, filename::String="zonal_pv_dispatch.png")
    
    T = data[:T]
    Tset = 1:T
    kVA_B = data[:kVA_B]
    p_D_pu = data[:p_D_pu]
    Dset = data[:Dset]
    
    zone_assignments = get_zone_assignments(data)
    pv_by_zone = zone_assignments[:pv_by_zone]
    
    # Aggregate PV by zone
    p_pv_zone_kW = Dict{String, Vector{Float64}}()
    
    for (zone, pv_buses) in pv_by_zone
        p_zone = zeros(T)
        for bus in pv_buses
            if bus in Dset
                p_zone .+= p_D_pu[bus, :] .* kVA_B  # Convert to kW
            end
        end
        p_pv_zone_kW[zone] = p_zone
    end
    
    # Colors matching your scheme
    zone_colors = [
        RGB(0x00/255, 0x72/255, 0xB2/255),    # Zone 1: deep blue
        RGB(0xE6/255, 0x9F/255, 0x00/255),    # Zone 2: amber
        RGB(0x00/255, 0x9E/255, 0x73/255),    # Zone 3: teal
        RGB(0xCC/255, 0x79/255, 0xA7/255),    # Zone 4: magenta
        RGB(0xD5/255, 0x5E/255, 0x00/255)     # Zone 5: rust red
    ]
    
    marker_shapes = [:circle, :square, :diamond, :utriangle, :star5]
    
    # Create plot
    gr()
    theme(:dao)
    
    p = plot(
        dpi=400,
        xlabel="Time Period (t)",
        ylabel="PV Real Power [kW]",
        legend=:topright,
        legend_columns=2,
        legendfontsize=8,
        legend_background_color=RGBA(1,1,1,0.7),
        framestyle=:box,
        grid=true,
        gridlinewidth=0.5,
        gridalpha=0.25,
        gridstyle=:solid,
        titlefont=font(11, "Computer Modern"),
        guidefont=font(10, "Computer Modern"),
        tickfont=font(9, "Computer Modern"),
        left_margin=3Plots.mm,
        right_margin=3Plots.mm,
        title="Zonal PV Dispatch (Slack: $slack_sub)"
    )
    
    # Plot each zone
    zones = ["1s", "2s", "3s", "4s", "5s"]
    for (i, zone) in enumerate(zones)
        n_pv = length(pv_by_zone[zone])
        plot!(
            p, Tset, p_pv_zone_kW[zone],
            color=zone_colors[i],
            lw=2.2,
            markershape=marker_shapes[i],
            markersize=5,
            markerstrokewidth=1.5,
            markerstrokecolor=:black,
            markeralpha=0.9,
            label="Zone $zone ($n_pv PV)"
        )
    end
    
    if showPlots
        display(p)
    end
    
    if savePlots
        savefig(p, filename)
    end
    
    return p
end

"""
    plot_zonal_battery_actions(result, data, slack_sub; showPlots=true, savePlots=false, filename="zonal_battery_actions.png")

Plot aggregated battery power and energy for each zone.
Two subplots:
- Top: Charging (green) and discharging (wine red) power
- Bottom: State of charge (SOC)

# Arguments
- result: OPF solution result
- data: System data dictionary  
- slack_sub: Slack substation identifier
- showPlots: Whether to display the plot
- savePlots: Whether to save the plot
- filename: Filename for saving
"""
function plot_zonal_battery_actions(result, data, slack_sub; showPlots::Bool=true, savePlots::Bool=false, filename::String="zonal_battery_actions.png")
    
    T = data[:T]
    Tset = 0:T  # Include t=0 for initial SOC
    kVA_B = data[:kVA_B]
    E_BASE = kVA_B * 1.0
    Bset = data[:Bset]
    B0_pu = data[:B0_pu]
    B_R = data[:B_R]
    
    zone_assignments = get_zone_assignments(data)
    batt_by_zone = zone_assignments[:batt_by_zone]
    
    # Aggregate battery power and energy by zone
    P_batt_zone_kW = Dict{String, Vector{Float64}}()
    B_zone_kWh = Dict{String, Vector{Float64}}()
    B_rated_zone_kWh = Dict{String, Float64}()
    
    for (zone, batt_buses) in batt_by_zone
        P_zone = zeros(T)
        B_zone = zeros(T+1)  # Include t=0
        B_rated = 0.0
        
        for bus in batt_buses
            if bus in Bset
                # Power (t=1 to T)
                for t in 1:T
                    P_zone[t] += result[:P_B][bus, t] * kVA_B
                end
                
                # Energy (t=0 to T)
                B_zone[1] += B0_pu[bus] * E_BASE  # Initial SOC
                for t in 1:T
                    B_zone[t+1] += result[:B][bus, t] * E_BASE
                end
                
                # Rated capacity
                B_rated += B_R[bus]
            end
        end
        
        P_batt_zone_kW[zone] = P_zone
        B_zone_kWh[zone] = B_zone
        B_rated_zone_kWh[zone] = B_rated
    end
    
    # Colors
    zone_colors = [
        RGB(0x00/255, 0x72/255, 0xB2/255),    # Zone 1: deep blue
        RGB(0xE6/255, 0x9F/255, 0x00/255),    # Zone 2: amber
        RGB(0x00/255, 0x9E/255, 0x73/255),    # Zone 3: teal
        RGB(0xCC/255, 0x79/255, 0xA7/255),    # Zone 4: magenta
        RGB(0xD5/255, 0x5E/255, 0x00/255)     # Zone 5: rust red
    ]
    
    marker_shapes = [:circle, :square, :diamond, :utriangle, :star5]
    
    # Create plots
    gr()
    theme(:dao)
    
    # Top panel: Charging/Discharging power
    p_power = plot(
        dpi=400,
        ylabel="Battery Power [kW]",
        legend=:topright,
        legend_columns=2,
        legendfontsize=8,
        legend_background_color=RGBA(1,1,1,0.7),
        framestyle=:box,
        grid=true,
        gridlinewidth=0.5,
        gridalpha=0.25,
        gridstyle=:solid,
        titlefont=font(11, "Computer Modern"),
        guidefont=font(10, "Computer Modern"),
        tickfont=font(9, "Computer Modern"),
        left_margin=3Plots.mm,
        right_margin=3Plots.mm,
        top_margin=2Plots.mm,
        bottom_margin=0Plots.mm,
        xformatter=_->"",  # Hide x-labels
        title="Zonal Battery Actions (Slack: $slack_sub)"
    )
    
    # Plot power for each zone
    zones = ["1s", "2s", "3s", "4s", "5s"]
    for (i, zone) in enumerate(zones)
        n_batt = length(batt_by_zone[zone])
        if n_batt > 0
            plot!(
                p_power, 1:T, P_batt_zone_kW[zone],
                color=zone_colors[i],
                lw=2.2,
                markershape=marker_shapes[i],
                markersize=5,
                markerstrokewidth=1.5,
                markerstrokecolor=:black,
                markeralpha=0.9,
                label="Zone $zone ($n_batt Batt)"
            )
        end
    end
    
    hline!(p_power, [0], color=:black, lw=2, linestyle=:dash, label=false)
    
    # Bottom panel: State of charge
    p_soc = plot(
        dpi=400,
        xlabel="Time Period (t)",
        ylabel="Battery Energy [kWh]",
        legend=:topright,
        legend_columns=2,
        legendfontsize=8,
        legend_background_color=RGBA(1,1,1,0.7),
        framestyle=:box,
        grid=true,
        gridlinewidth=0.5,
        gridalpha=0.25,
        gridstyle=:solid,
        guidefont=font(10, "Computer Modern"),
        tickfont=font(9, "Computer Modern"),
        left_margin=3Plots.mm,
        right_margin=3Plots.mm,
        top_margin=0Plots.mm,
        bottom_margin=4Plots.mm
    )
    
    # Plot SOC for each zone
    for (i, zone) in enumerate(zones)
        n_batt = length(batt_by_zone[zone])
        if n_batt > 0
            plot!(
                p_soc, Tset, B_zone_kWh[zone],
                color=zone_colors[i],
                lw=2.2,
                markershape=marker_shapes[i],
                markersize=5,
                markerstrokewidth=1.5,
                markerstrokecolor=:black,
                markeralpha=0.9,
                label="Zone $zone ($n_batt Batt)"
            )
        end
    end
    
    # Combine panels
    p_combined = plot(p_power, p_soc, layout=(2,1), size=(7.1*100, 4.6*100))
    
    if showPlots
        display(p_combined)
    end
    
    if savePlots
        savefig(p_combined, filename)
    end
    
    return p_combined
end

"""
    print_zone_statistics(data)

Print summary statistics of resource distribution by zone.
"""
function print_zone_statistics(data)
    zone_assignments = get_zone_assignments(data)
    pv_by_zone = zone_assignments[:pv_by_zone]
    batt_by_zone = zone_assignments[:batt_by_zone]
    
    println("\n" * "="^80)
    println("ZONE RESOURCE DISTRIBUTION")
    println("="^80)
    
    zones = ["1s", "2s", "3s", "4s", "5s"]
    zone_names = ["SW", "NW", "N", "NE", "SE"]
    
    println("\n  Zone  | Location | PV Units | Battery Units | Load Pattern")
    println("  " * "-"^70)
    
    for (i, zone) in enumerate(zones)
        n_pv = length(pv_by_zone[zone])
        n_batt = length(batt_by_zone[zone])
        
        load_pattern = if zone == "1s"
            "Residential (sinusoid)"
        elseif zone == "2s"
            "Commercial (sinusoid)"
        elseif zone == "3s"
            "Industrial (3-shift)"
        elseif zone == "4s"
            "Office (on/off)"
        else
            "Charging station (ramp)"
        end
        
        @printf("  %-5s | %-8s | %8d | %13d | %s\n", zone, zone_names[i], n_pv, n_batt, load_pattern)
    end
    
    println("  " * "-"^70)
    println("  Total |          | %8d | %13d |", sum(length(pv_by_zone[z]) for z in zones), 
                                                 sum(length(batt_by_zone[z]) for z in zones))
    println("="^80)
end

end  # module
