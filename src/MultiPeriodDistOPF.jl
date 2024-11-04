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
using .Plotter

include("./exporter.jl")
using .Exporter

include("./openDSSValidator.jl")
using .openDSSValidator

# Re-export the function from parseOpenDSSFiles
# export compute_output_values, export_decision_variables, export_optimization_model, export_simulation_key_results_txt, get_scd, myprintln, parse_all_data, parse_system_simulation_data, parse_branch_data, parse_load_data, parse_pv_data, parse_battery_data, evaluate_voltage_limits, generateBinaryLoadShape, plot_battery_actions, plot_input_forecast_curves, plot_line_losses, plot_substation_power, plot_substation_power_cost, validate_opf_against_opendss

export compute_output_values,
    evaluate_voltage_limits,
    export_decision_variables,
    export_optimization_model,
    export_simulation_key_results_txt,
    generateBinaryLoadShape,
    get_scd,
    get_source_bus,
    get_substation_lines,
    myprintln,
    parse_all_data,
    parse_battery_data,
    parse_branch_data,
    parse_load_data,
    parse_pv_data,
    parse_system_simulation_data,
    plot_battery_actions,
    plot_input_forecast_curves,
    plot_line_losses,
    plot_substation_power,
    plot_substation_power_cost,
    validate_opf_against_opendss

end # module
