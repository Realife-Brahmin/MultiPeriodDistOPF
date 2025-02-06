# solveMPOPF.jl
using Revise
using MultiPeriodDistOPF
using Parameters

Revise.revise()

systemName = "ads10_1ph"
# systemName = "ieee123_1ph"
# T0 = 7
T0 = 24
# T0 = 11
# factor = 1/2
factor = 4
T = Int(T0*factor) 
numAreas = 1
temporal_decmp = false
temporal_decmp = true
maxiter_ddp = 33
savePlots = false
savePlots = true
# objfun0 = "powerflow"
# objfun0 = "lineLossMin"
objfun0 = "subsPowerCostMin"
objfun2 = "scd"
inputForecastDescription = "bilevelCosts"
relax_terminal_soc_constraint = false
relax_terminal_soc_constraint = true
tSOC_hard = false
# tSOC_hard = true
PSubsMax_kW = Inf # Inf means no limit
solver = "Ipopt"
# solver = "Gurobi"
# solver = "Juniper"

# Parse all data
data = parse_all_data(systemName, T, temporal_decmp=temporal_decmp, relax_terminal_soc_constraint=relax_terminal_soc_constraint)

if !temporal_decmp
    modelDict = optimize_MPOPF_1ph_NL_TemporallyBruteforced(data)
    @unpack model, modelVals, data = modelDict
    # Print mu values
    # Print mu values
    print_mu(modelDict)
    # @unpack mu = modelDict
    # @unpack Tset, Bset = data
    # print_mu(mu, Tset, Bset)
elseif temporal_decmp
    modelDict = optimize_MPOPF_1ph_NL_DDP(data, maxiter=maxiter_ddp) # modelDict is basically ddpModel
    @unpack modelVals, data = modelDict
    # Print mu values
    print_mu(modelDict)
else
    error("temporal_decmp must be either true or false")
end

if temporal_decmp
    plot_substation_power_cost_allT_vs_k(modelDict, savePlots=savePlots)
end

# postsim computation, plotting, logging
begin

    verbose = false
    # verbose = true

    modelDict = compute_output_values(modelDict, verbose=verbose)
    # @unpack data = modelDict

    export_decision_variables(modelDict, verbose=verbose)

    # Todo: Maybe separately save the simulation times? It is annoying to have file content differences every single run (for same exact sim)
    export_simulation_key_results_txt(modelDict, verbose=verbose)

    plot_battery_actions(modelDict, showPlots=false, savePlots=savePlots, verbose=verbose)

    @unpack data = modelDict;
    plot_input_forecast_curves(data, filenameSuffix=inputForecastDescription, showPlots=false, verbose=verbose)
    plot_substation_power(modelDict, savePlots=savePlots, verbose=verbose)
    plot_substation_power_cost(modelDict, savePlots=savePlots, verbose=verbose)
    plot_line_losses(modelDict, savePlots=savePlots, verbose=verbose)

    modelDict = validate_opf_against_opendss(modelDict, verbose=verbose)

    export_validation_key_results(modelDict, verbose=verbose)

    export_validation_decision_variables(modelDict, verbose=verbose)
end

