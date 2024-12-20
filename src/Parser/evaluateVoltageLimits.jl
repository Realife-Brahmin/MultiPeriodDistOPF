# evaluateVoltageLimits.jl
module evaluateVoltageLimits

export evaluate_voltage_limits

#region evaluate_voltage_limits
"""
    evaluate_voltage_limits(load_data, pv_data, battery_data)

Evaluate the voltage limits for a given set of load, PV, and battery data.

This function calculates the most restrictive voltage limits for each bus in the combined set of load, PV, and battery components. 
It determines the highest lower limit and the lowest upper limit for each bus.

# Arguments
- `load_data::Dict`: A dictionary containing load data, including voltage limits.
- `pv_data::Dict`: A dictionary containing PV data, including voltage limits.
- `battery_data::Dict`: A dictionary containing battery data, including voltage limits.

# Returns
- `component_data::Dict`: A dictionary containing the combined component set and the calculated voltage limits.
"""
function evaluate_voltage_limits(load_data, pv_data, battery_data)
    # Create Compset as the union of the sets from the three components
    Compset = union(load_data[:NLset], pv_data[:Dset], battery_data[:Bset])

    # Initialize dictionaries to store Vminpu_Comp and Vmaxpu_Comp
    Vminpu_Comp = Dict()
    Vmaxpu_Comp = Dict()

    # Iterate over each bus in Compset to calculate voltage limits
    for j in Compset
        # Get Vminpu and Vmaxpu for load, PV, and battery components if they exist
        Vminpu_L = get(load_data[:Vminpu_L], j, -Inf)   # Default to Inf if not in set
        Vmaxpu_L = get(load_data[:Vmaxpu_L], j, Inf)  # Default to -Inf if not in set

        Vminpu_D = get(pv_data[:Vminpu_D], j, -Inf)
        Vmaxpu_D = get(pv_data[:Vmaxpu_D], j, Inf)

        Vminpu_B = get(battery_data[:Vminpu_B], j, -Inf)
        Vmaxpu_B = get(battery_data[:Vmaxpu_B], j, Inf)

        # Calculate the most restrictive voltage limits for this bus
        Vminpu_Comp[j] = maximum([Vminpu_L, Vminpu_D, Vminpu_B])  # Choose the highest lower limit
        Vmaxpu_Comp[j] = minimum([Vmaxpu_L, Vmaxpu_D, Vmaxpu_B])  # Choose the lowest upper limit
    end

    # Create the output dictionary
    component_data = Dict(
        :Compset => Compset,
        :Vminpu_Comp => Vminpu_Comp,
        :Vmaxpu_Comp => Vmaxpu_Comp
    )

    return component_data
end
#endregion

end # module evaluateVoltageLimits