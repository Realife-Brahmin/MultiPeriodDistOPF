module MultiPeriodDistOPF

# Export submodules
export computeOutputs
export DDP
export functionRetriever
export helperFunctions
export ModelBuilder
export Playbook_of_MPOPF
export parseOpenDSSFiles
export Plotter
export Exporter
export openDSSValidator

using JuMP
using Parameters: @unpack, @pack!

# Include and use submodules
include("./computeOutputs.jl")
include("./DDP/DDP.jl")
include("./DDP/DDPLinear.jl")
include("./functionRetriever.jl")
include("./helperFunctions.jl")
include("./ModelBuilder/ModelBuilder.jl")
include("./ModelCopier/ModelCopier.jl")
include("./playbook_of_mpopf.jl")
include("./Parser/parseOpenDSSFiles.jl")
include("./Plotter/Plotter.jl")
include("./exporter.jl")
include("./openDSSValidator.jl")

# Use your own submodules
using .computeOutputs
using .DDP
using .DDPLinear
using .functionRetriever
using .helperFunctions
using .ModelBuilder
using .Playbook_of_MPOPF
using .parseOpenDSSFiles
using .Plotter
using .Exporter
using .openDSSValidator

if isdefined(Main, :Revise)
    @info "Tracking submodules with Revise"
    Revise.track(@__MODULE__)
end

# Explicitly import and re-export macros from external packages (correct way)
export @unpack, @pack!

# Explicitly export your submodule functions you want accessible externally
export optimize_MPOPF_1ph_NL_TemporallyBruteforced,
    optimize_MPOPF_1ph_L,
    optimize_MPOPF_1ph_NL_DDP,
    optimize_MPOPF_1ph_L_DDP,
    parse_all_data,
    compute_output_values,
    export_decision_variables,
    export_simulation_key_results_txt,
    export_validation_decision_variables,
    export_validation_key_results,
    plot_battery_actions,
    plot_input_forecast_curves,
    plot_line_losses,
    plot_substation_power,
    plot_substation_power_cost,
    plot_substation_power_cost_allT_vs_k,
    validate_opf_against_opendss

end # module MultiPeriodDistOPF
