# solveMPOPF.jl
using Revise
using MultiPeriodDistOPF
using Parameters: @unpack

Revise.revise()

systemName = "ads10_1ph"
# systemName = "ieee123_1ph"
T0 = 24
# factor = 1/8
factor = 1
T = Int(T0*factor) 
numAreas = 1
temporal_decmp = false
# objfun0 = "powerflow"
# objfun0 = "lineLossMin"
objfun0 = "subsPowerCostMin"
# objfun2 = "none"
objfun2 = "scd"
inputForecastDescription = "bilevelCosts"
alpha = 1e-3
PSubsMax_kW = Inf # Inf means no limit
solver = "Ipopt"
# solver = "EAGO"
# solver = "Gurobi"
# solver = "Juniper"
# solver = "MadNLP"

# Parse all data
data = parse_all_data(systemName, T, numAreas=numAreas, alpha=alpha, objfun0=objfun0, objfun2=objfun2,temporal_decmp=temporal_decmp, PSubsMax_kW=PSubsMax_kW, inputForecastDescription=inputForecastDescription, solver=solver)

# model = optimize_MPOPF_1ph_NL(data)
model = optimize_MPOPF_1ph_NL_TemporallyBruteforced(data)

# postsim computation, plotting, logging
begin

    verbose = false
    # verbose = true

    export_optimization_model(model, data, verbose=verbose)

    data = compute_output_values(model, data, verbose=verbose)

    export_decision_variables(model, data, verbose=verbose)

    # Todo: Maybe separately save the simulation times? It is annoying to have file content differences every single run (for same exact sim)
    export_simulation_key_results_txt(model, data, verbose=verbose)

    savePlots = false
    savePlots = true

    plot_battery_actions(model, data, showPlots=false, savePlots=savePlots, verbose=verbose)

    plot_input_forecast_curves(data, filenameSuffix=inputForecastDescription, showPlots=false, verbose=verbose)
    plot_substation_power(data, savePlots=savePlots, verbose=verbose)
    plot_substation_power_cost(data, savePlots=savePlots, verbose=verbose)
    plot_line_losses(data, savePlots=savePlots, verbose=verbose)

    vald = validate_opf_against_opendss(model, data, verbose=verbose)

    export_validation_key_results(vald, data, verbose=verbose)

    export_validation_decision_variables(vald, data, verbose=verbose)
end

