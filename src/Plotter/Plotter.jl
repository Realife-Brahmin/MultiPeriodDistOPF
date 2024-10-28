module Plotter

using LaTeXStrings
using Plots
using Parameters: @unpack
import JuMP: value  # Import JuMP's value function to extract values of decision variables
import Base.Filesystem: mkpath, isdir  # To create directories

plot_font = "Computer Modern"
default(fontfamily=plot_font,
    linewidth=2, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1.3)

export plot_battery_actions, plot_substation_power

function plot_battery_actions(model, data;
    showPlots::Bool=false,
    savePlots::Bool=true,
    macroItrNum::Int=1)

    # Extract necessary parameters from the `data` dictionary
    @unpack Tset, Bset, kVA_B, B_R_pu, P_B_R, Bref_pu, systemName, numAreas, T, DER_percent, Batt_percent, alpha = data;

    Tset = sort(collect(Tset))

    # @unpack P_c, P_d, B = model; # Why does this not work?

    P_c = model[:P_c]
    P_d = model[:P_d]
    B = model[:B]

    # Set the base directory for saving plots
    base_dir = joinpath("processedData", systemName, "numAreas_$(numAreas)",
        "batteryActionPlots", "Horizon_$(T)",
        "pv_$(DER_percent)_batt_$(Batt_percent)",
        "macroItr_$(macroItrNum)")

    # Create the base directory if it does not exist
    if savePlots && !isdir(base_dir)
        println("Creating directory: $base_dir")
        mkpath(base_dir)  # This will create the directory and its parents if needed
    end

    # Loop through all battery buses and create plots
    for j in Bset
        time_intervals = collect(Tset)

        # Collect charging, discharging power, and state of charge values
        charging_power_kW = [value(P_c[j, t]) * kVA_B for t in Tset]
        discharging_power_kW = [value(P_d[j, t]) * kVA_B for t in Tset]

        # State of charge percentage, starting with Bref_pu for t=0 and continuing with B[t]
        soc = [Bref_pu[j] / B_R_pu[j] * 100; [value(B[j, t]) / B_R_pu[j] * 100 for t in Tset]]

        # Use P_B_R[j] to set the y-limits for charging/discharging
        ylimit = (-P_B_R[j], P_B_R[j])

        # Create a plot for charging and discharging
        charging_discharge_plot = bar(time_intervals, charging_power_kW, label="Charging", color=:green, legend=:bottomleft, xlabel="Time Interval Number", ylabel="[kW]", ylim=ylimit, xticks=1:T)
        bar!(time_intervals, -discharging_power_kW, label="Discharging", color=:darkred)

        # Add a horizontal line at y=0 to clearly separate charging and discharging
        hline!(charging_discharge_plot, [0], color=:black, lw=2, label=false)  # Add prominent black line at y=0

        # Add a title for each battery bus
        bus_label = "Bus $j"
        title!(charging_discharge_plot, "Battery at $bus_label\nCharging and Discharging")

        # Create a plot for SOC, starting from t=0 with Bref_pu and continuing with the values in SOC
        soc_plot = bar(0:T, soc, label="Battery State of Charge", color=:purple, legend=:bottomleft, xlabel="Time Interval Number", ylabel="[%]", ylim=(0, 100), xticks=0:T)
        title!(soc_plot, "SOC")

        # Combine the two plots in a layout
        plot_combined = plot(charging_discharge_plot, soc_plot, layout=@layout([a; b]))

        # Show the plot if `showPlots` is true
        if showPlots
            display(plot_combined)
        end

        # Save the plot if `savePlots` is true
        if savePlots
            # Create the filename with the required format
            filename = joinpath(base_dir, "Battery_$(j)_alpha_$(alpha).png")
            println("Saving plot to: $filename")
            png(plot_combined, filename)
        end
    end
end

# Function to plot Substation Power over time
function plot_substation_power(data)
    @unpack Tset, PSubs_vs_t_1toT_kW, T = data  # Assuming T is the last time period

    yvalues = PSubs_vs_t_1toT_kW;

    gr()
    plot(
        Tset, PSubs_vs_t_1toT_kW,
        label=L"(P^t_{Subs})",
        xlabel="Time Period " * L"t",
        ylabel=L"P_{Subs} \, [kW]",
        title="Substation Power " * L"(P_{Subs})" * " across the Horizon",
        legend=:topright,
        gridstyle=:solid,
        gridlinewidth=1.0,
        gridalpha=0.2,
        minorgrid=true,
        minorgridstyle=:solid,
        minorgridalpha=0.05,
        lw=2,
        marker=:circle,
        markersize=4,
        xlims=(0, T+1),
        xticks=1:1:T,  # Set xticks for every hour from 1 to T
        ylims=(minimum(yvalues)*0.95, maximum(yvalues)*1.05), 
        titlefont=font(12, "Computer Modern"),
        guidefont=font(15, "Computer Modern"),
        tickfontfamily="Computer Modern"
    )
end


end  # module Plotter
