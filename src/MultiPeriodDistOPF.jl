module MultiPeriodDistOPF

include("./computeOutputs.jl")
using .computeOutputs: compute_output_values

include("./functionRetriever.jl")
# using .functionRetriever: get_scd
using .functionRetriever

include("./helperFunctions.jl")
using .helperFunctions: myprintln

include("./Parser/parseOpenDSSFiles.jl")
using .parseOpenDSSFiles

# Import the plotter module
include("./Plotter/Plotter.jl")
using .Plotter: plot_battery_actions, plot_line_losses, plot_substation_power, plot_substation_power_cost

include("./exporter.jl")
using .Exporter: export_decision_variables, export_simulation_key_results_txt

# Re-export the function from parseOpenDSSFiles
export compute_output_values, export_decision_variables, export_simulation_key_results_txt, get_scd, myprintln, parse_all_data, parse_system_simulation_data, parse_branch_data, parse_load_data, parse_pv_data, parse_battery_data, evaluate_voltage_limits, generateBinaryLoadShape, plot_battery_actions, plot_line_losses, plot_substation_power, plot_substation_power_cost

end # module
