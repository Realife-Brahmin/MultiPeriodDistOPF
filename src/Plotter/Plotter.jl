module Plotter

using Plots
import JuMP: value  # Import JuMP's value function to extract values of decision variables
import Base.Filesystem: mkdir, isdir  # To create directories

export plot_battery_actions

function plot_battery_actions(model, data;
    showPlots::Bool=false,
    savePlots::Bool=true,
    macroItrNum::Int=1)

    # Extract necessary parameters from the `data` dictionary
    Tset = data[:Tset]
    Bset = data[:Bset]
    P_c = model[:P_c]
    P_d = model[:P_d]
    B = model[:B]
    B_R_pu = data[:B_R_pu]  # Rated storage capacity (for SOC % conversion)

    systemName = data[:systemName]
    numAreas = data[:numAreas]
    T = data[:T]
    DER_percent = data[:DER_percent]
    Batt_percent = data[:Batt_percent]
    alpha = data[:alpha]

    # Set the base directory for saving plots
    base_dir = joinpath("processedData", systemName, "numAreas_$(numAreas)",
        "batteryActionPlots", "Horizon_$(T)",
        "pv_$(DER_percent)_batt_$(Batt_percent)",
        "macroItr_$(macroItrNum)")

    # # Create the base directory if it does not exist
    # if savePlots && !isdir(base_dir)
    #     println("Creating directory: $base_dir")
    #     mkdir(base_dir; recursive=true)
    # end

    # Create the base directory if it does not exist
    if savePlots && !isdir(base_dir)
        println("Creating directory: $base_dir")
        mkpath(base_dir)  # This will create the directory and its parents if needed
    end

    # Loop through all battery buses and create plots
    for j in Bset
        time_intervals = collect(Tset)

        # Collect charging, discharging power, and state of charge values
        charging_power = [value(P_c[j, t]) for t in Tset]
        discharging_power = [value(P_d[j, t]) for t in Tset]

        # State of charge percentage, converted by dividing by the rated capacity in pu
        soc = [value(B[j, t]) / B_R_pu[j] * 100 for t in Tset]

        # Create a plot for charging and discharging
        charging_discharge_plot = bar(time_intervals, charging_power, label="Charging", color=:green, legend=:top, xlabel="Time Interval Number", ylabel="[kW]", ylim=(-5, 5))
        bar!(time_intervals, -discharging_power, label="Discharging", color=:purple)

        # Add a title for each battery bus
        bus_label = "Bus $j"
        title!(charging_discharge_plot, "Battery at $bus_label\nCharging and Discharging")

        # Create a plot for SOC
        soc_plot = bar(time_intervals, soc, label="Battery State of Charge", color=:purple, legend=:bottom, xlabel="Time Interval Number", ylabel="[%]", ylim=(0, 100))
        title!(soc_plot, "SOC at $bus_label")

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

end  # module Plotter
