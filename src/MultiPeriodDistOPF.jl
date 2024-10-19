module MultiPeriodDistOPF

include("./parseOpenDSSFiles.jl")
using .parseOpenDSSFiles

include("./Plotter/Plotter.jl")
using .Plotter: plot_battery_actions  # Import the plotter module

# Re-export the function from parseOpenDSSFiles
export myprintln, parse_all_data, parse_system_simulation_data, parse_branch_data, parse_load_data, parse_pv_data, parse_battery_data, evaluate_voltage_limits, generateBinaryLoadShape, plot_battery_actions

end
