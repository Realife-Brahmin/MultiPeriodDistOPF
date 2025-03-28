module Exporter

export export_decision_variables, 
    export_optimization_model,
    export_simulation_key_results_txt,
    export_validation_decision_variables,
    export_validation_key_results

using CSV
using DelimitedFiles  # To write CSV files
using JuMP
using XLSX
using Parameters
include("./helperFunctions.jl")
import .helperFunctions as HF

#region export_optimization_model
"""
    export_optimization_model(modelDict; verbose::Bool=false)

Export the optimization model's decision variables and results.

This function extracts the optimization model from `modelDict` and exports its decision variables and results. 
If temporal decomposition is enabled, it handles the forward step and appends the appropriate suffix to the filenames.
"""
function export_optimization_model(modelDict;
    verbose::Bool=false)
    @unpack data = modelDict
    @unpack temporal_decmp = data
    ddp_appendix = ""
    if !temporal_decmp
        @unpack model = modelDict
    elseif temporal_decmp
        ddpModel = modelDict
        @unpack models_ddp_vs_t_vs_k, t_ddp, k_ddp = ddpModel
        model_t0 = models_ddp_vs_t_vs_k[t_ddp, k_ddp]
        model = model_t0
        ddp_appendix = "_forward_step_t_$(t_ddp)_k_$(k_ddp)"
    else
        @error "temporal_decmp must be either true or false"
    end
    # Define the path and filename based on the specified structure
    @unpack T, systemName, numAreas, gedAppendix, machine_ID, objfunAppendix, simNatureAppendix, linearizedModelAppendix = data
    base_dir = joinpath("processedData", systemName, gedAppendix, "Horizon_$(T)", "numAreas_$(numAreas)", "models")

    # Create the directory if it doesn't exist
    if !isdir(base_dir)
        println("Creating directory: $base_dir")
        mkpath(base_dir)
    end

    # Define the filename with the appropriate structure
    filename = joinpath(base_dir, "optimizationModel_for_$(objfunAppendix)_via_$(simNatureAppendix)_$(linearizedModelAppendix)" * ddp_appendix * ".txt")

    # Check if the file already exists, and delete it if so
    if isfile(filename)
        rm(filename)
    end

    # Open a new file and write the model contents to it
    open(filename, "w") do f

        # Write model constraints
        println(f, "\nConstraints:")
        for constr in JuMP.all_constraints(model, include_variable_in_set_constraints=true)
            println(f, constr)
        end

        # Write model objective
        println(f, "\nObjective:")
        println(f, JuMP.objective_function(model))
    end

    # myprintln(verbose, "Model successfully written to $filename")
end
#endregion

#region export_decision_variables
"""
    export_decision_variables(modelDict; verbose::Bool=false)

Export the decision variables of the optimization model to CSV files.

This function extracts the decision variables from `modelDict` and writes them to CSV files for further analysis.
"""
function export_decision_variables(modelDict;
    filename::String="decision_variables.csv",
    verbose::Bool=false)
    
    @unpack modelVals, data = modelDict;
    # Define the path and filename based on the specified structure
    @unpack T, systemName, numAreas, gedAppendix, machine_ID, objfunAppendix, simNatureAppendix, solver, linearizedModelAppendix = data
    base_dir = joinpath("processedData", systemName, gedAppendix, "Horizon_$(T)", "numAreas_$(numAreas)", "decvars")

    if !isdir(base_dir)
        println("Creating directory: $base_dir")
        mkpath(base_dir)
    end

    ext = ".csv"
    filename = joinpath(base_dir, "decisionVariables_for_$(objfunAppendix)_via_$(simNatureAppendix)_with_$(linearizedModelAppendix)"*ext)
    
    # Log current working directory
    HF.myprintln(verbose, "Current working directory: $(pwd())")

    # Ensure the filename has a valid path or is saved in the current directory
    HF.myprintln(verbose, "Saving to filename: $filename")

    # Unpack the necessary data from `data`
    Tset = sort(collect(data[:Tset]))  # Time steps
    LoadShapeLoad = data[:LoadShapeLoad]  # Load shape values
    LoadShapePV = data[:LoadShapePV]  # Irradiance values
    LoadShapeCost = data[:LoadShapeCost] * 100  # Convert from $/kWh to cents/kWh

    Lset = sort(collect(data[:Lset]))  # Line set
    Nset = sort(collect(data[:Nset]))  # Node set
    Bset = sort(collect(data[:Bset]))  # Battery set
    Dset = sort(collect(data[:Dset]))  # PV set

    P_Subs = modelVals[:P_Subs]
    P = modelVals[:P]
    Q = modelVals[:Q]
    l = modelVals[:l]
    v = modelVals[:v]
    q_D = modelVals[:q_D]
    q_B = modelVals[:q_B]
    P_c = modelVals[:P_c]
    P_d = modelVals[:P_d]
    B = modelVals[:B]

    # Prepare data to write in CSV format as a matrix of strings
    data_matrix = []

    # First row: Time step "t"
    push!(data_matrix, ["t"; Tset...])

    # Second row: LoadShapeLoad values (lambda)
    push!(data_matrix, ["lambda"; LoadShapeLoad...])

    # Third row: Irradiance (LoadShapePV)
    push!(data_matrix, ["Irrad"; LoadShapePV...])

    # Fourth row: cents/kWh (LoadShapeCost)
    push!(data_matrix, ["cents/kWh"; LoadShapeCost...])

    # Next row: P_Subs for all time intervals
    push!(data_matrix, ["P_Subs"; [P_Subs[t] for t in Tset]...])

    # Next rows: P_ij for all (i, j) pairs in Lset
    for (i, j) in Lset
        push!(data_matrix, ["P_ij_$(i)_$(j)"; [P[(i, j), t] for t in Tset]...])
    end

    # Next rows: Q_ij for all (i, j) pairs in Lset
    for (i, j) in Lset
        push!(data_matrix, ["Q_ij_$(i)_$(j)"; [Q[(i, j), t] for t in Tset]...])
    end

    # Next rows: l_ij for all (i, j) pairs in Lset
    for (i, j) in Lset
        push!(data_matrix, ["l_ij_$(i)_$(j)"; [l[(i, j), t] for t in Tset]...])
    end

    # Next rows: v_j for all nodes in Nset
    for j in Nset
        push!(data_matrix, ["v_j_$(j)"; [v[j, t] for t in Tset]...])
    end

    # Next rows: q_D_j for all PV buses in Dset
    for j in Dset
        push!(data_matrix, ["q_D_j_$(j)"; [q_D[j, t] for t in Tset]...])
    end

    # Next rows: q_B_j for all battery buses in Bset
    for j in Bset
        push!(data_matrix, ["q_B_j_$(j)"; [q_B[j, t] for t in Tset]...])
    end

    # Next rows: P_c_j for all battery buses in Bset
    for j in Bset
        push!(data_matrix, ["P_c_j_$(j)"; [P_c[j, t] for t in Tset]...])
    end

    # Next rows: P_d_j for all battery buses in Bset
    for j in Bset
        push!(data_matrix, ["P_d_j_$(j)"; [P_d[j, t] for t in Tset]...])
    end

    # Next rows: B_j for all battery buses in Bset
    for j in Bset
        push!(data_matrix, ["B_j_$(j)"; [B[j, t] for t in Tset]...])
    end

    # Write the matrix to the CSV file
    try
        HF.myprintln(verbose, "Opening file: $filename")
        writedlm(filename, data_matrix, ',')
        HF.myprintln(verbose, "Data written successfully.")
    catch e
        HF.myprintln(verbose, "Error occurred: $e")
    end

    # Check if file was created
    if isfile(filename)
        HF.myprintln(verbose, "File exists: $filename")
    else
        HF.myprintln(verbose, "File does not exist: $filename")
    end

    HF.myprintln(verbose, "Decision variables exported to $filename")
end
#endregion

#region export_simulation_key_results_txt
"""
    export_simulation_key_results_txt(modelDict; verbose::Bool=false)

Export key results of the simulation to a text file.

This function extracts key results from `modelDict` and writes them to a text file for documentation and review.
"""
function export_simulation_key_results_txt(modelDict; filename::String="simulation_results.txt", verbose::Bool=false)

    @unpack data = modelDict;
    # Define the path and filename based on the specified structure
    @unpack T, systemName, numAreas, gedAppendix, machine_ID, objfunConciseDescription, simNatureAppendix, solver, linearizedModelAppendix = data
    base_dir = joinpath("processedData", systemName, gedAppendix, "Horizon_$(T)", "numAreas_$(numAreas)")

    if !isdir(base_dir)
        HF.myprintln(verbose, "Creating directory: $base_dir")
        mkpath(base_dir)
    end

    filename = joinpath(base_dir, "simResults_for_$(objfunConciseDescription)_via_$(simNatureAppendix)_with_$(linearizedModelAppendix).txt")

    # Extract system information and parameters from `data`
    macroItrsCompleted = get(data, :macroItrsCompleted, 0)  # Default to 1 if not set
    solution_time = get(data, :solution_time, -1)  # Placeholder, real solution time

    # Open the file and write each section
    open(filename, "w") do f
        # Initialize output item counter
        item_counter = 1

        # Header Section
        println(f, "---------------------------------------------")
        println(f, "$(item_counter). Machine ID: $(data[:machine_ID])")
        item_counter += 1
        println(f, "$(item_counter). Solver Used: $(data[:solver])")
        item_counter += 1
        println(f, "$(item_counter). System Name: $(data[:systemName])")
        item_counter += 1
        println(f, "$(item_counter). Horizon Duration: $(data[:T])")
        item_counter += 1
        println(f, "$(item_counter). Nature of Optimization Simulation: $(data[:simNatureString])")
        item_counter += 1  # Placeholder
        println(f, "$(item_counter). Optimization Model: $(data[:linearizedModelString])")
        item_counter += 1  # Placeholder
        println(f, "$(item_counter). Objective: $(data[:objfunString])")
        item_counter += 1  # Placeholder
        println(f, "$(item_counter). GED Configuration: $(data[:gedAppendix])")
        item_counter += 1
        println(f, "$(item_counter). Maximum Substation Power Allowed: $(data[:PSubsMax_kW]) kW")
        item_counter += 1
        println(f, "---------------------------------------------")

        # Horizon Results
        println(f, "Full $(data[:T]) Hour Horizon")

        # Example metrics using the iterator
        println(f, "$(item_counter). Horizon Total Cost of Substation Power: \$ $(round(data[:PSubsCost_allT_dollar], digits=2))")
        item_counter += 1

        # Example metrics using the iterator
        println(f, "$(item_counter). Horizon Total Line Loss: $(round(data[:PLoss_allT_kW], digits=2)) kW")
        item_counter += 1

        println(f, "$(item_counter). Horizon Total Substation Power: $(round(data[:PSubs_allT_kW], digits=2)) kW + $(round(data[:QSubs_allT_kVAr], digits=2)) kVAr")
        item_counter += 1

        println(f, "$(item_counter). Horizon Total Load: $(round(data[:load_real_power_allT_kW], digits=2)) kW + $(round(data[:load_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1

        println(f, "$(item_counter). Horizon Total Generation: $(round(data[:total_gen_real_power_allT_kW], digits=2)) kW + $(round(data[:total_gen_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1

        # Placeholder for additional metrics
        # Example for static capacitor power generation
        println(f, "$(item_counter). Horizon Total Static Capacitor Reactive Power Generation: $(round(data[:static_cap_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1

        # PV-related metrics
        println(f, "$(item_counter). Horizon Total PV Generation: $(round(data[:pv_real_power_allT_kW], digits=2)) kW + $(round(data[:pv_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1

        # Battery-related metrics
        println(f, "$(item_counter). Horizon Total Battery Generation: $(round(data[:battery_real_power_allT_kW], digits=2)) kW + $(round(data[:battery_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1
        
        println(f, "$(item_counter). Horizon Total Battery Transaction Magnitude: $(round(data[:battery_real_power_transaction_magnitude_allT_kW], digits=2)) kW + $(round(data[:battery_reactive_power_transaction_magnitude_allT_kVAr], digits=2)) kVAr")
        item_counter += 1

        # Example for end horizon SCD observed and energy deviation
        println(f, "$(item_counter). Horizon Total SCD Observed: $(round(data[:scd_allT_kW], digits=2)) kW")
        item_counter += 1

        println(f, "$(item_counter). Horizon-end Battery Energy Deviation from Reference: $(round(data[:terminal_soc_violation_kWh], digits=2)) kWh")
        item_counter += 1

        # Example for substation power peak
        println(f, "$(item_counter). Horizon-Total All time Substation Power Peak: $(round(data[:substation_real_power_peak_allT_kW], digits=2)) kW")
        item_counter += 1

        # Additional Simulation Metadata
        println(f, "---------------------------------------------")
        println(f, "$(item_counter). Number of Macro-Iterations: $(macroItrsCompleted+1)")
        item_counter += 1
        if !data[:temporal_decmp]
            @unpack num_decvars, num_lincons, num_nonlincons = data;
            println(f, "$(item_counter). Number of Decision Variables: $(num_decvars)")
            item_counter += 1
            println(f, "$(item_counter). Number of Linear Constraints: $(num_lincons)")
            item_counter += 1
            println(f, "$(item_counter). Number of Nonlinear Constraints: $(num_nonlincons)")
            item_counter += 1
        end
        println(f, "$(item_counter). Simulation Time: $(round(solution_time, digits=2)) s")
        item_counter += 1
        println(f, "$(item_counter). Time to solve with sequential (non-parallel) computation: $(round(solution_time, digits=2)) s")
        item_counter += 1  # Placeholder for non-parallel time
        println(f, "$(item_counter). Time to solve if OPF computation parallelized: $(round(solution_time, digits=2)) s")  # Placeholder for parallelized time
        println(f, "---------------------------------------------")
    end

    HF.myprintln(verbose, "Simulation key results exported to $filename")
end
#endregion

#region export_validation_decision_variables
"""
    export_validation_decision_variables(modelDict; verbose::Bool=false)

Export the decision variables from the OpenDSS powerflow validation of the optimization's suggested control actions.

This function extracts the decision variables from `modelDict` and writes them to CSV files specifically for validation purposes.
The 'validation' involves using the 'control actions' stored in `modelDict` in an OpenDSS powerflow simulation, 
and checking if the resulting powerflow state variables match those solved by the optimizer.
"""
function export_validation_decision_variables(modelDict; verbose::Bool=false)

    @unpack valdVals, data = modelDict
    # Define the path and filename based on the specified structure
    @unpack T, systemName, numAreas, gedAppendix, machine_ID, objfunConciseDescription, processedDataFolderPath, simNatureAppendix, solver, linearizedModelAppendix = data
    base_dir = joinpath(processedDataFolderPath, systemName, gedAppendix, "Horizon_$(T)", "numAreas_$(numAreas)", "decvars")

    # Create the directory if it doesn't exist
    if !isdir(base_dir)
        if verbose
            println("Creating directory: $base_dir")
        end
        mkpath(base_dir)
    end

    # Define the filename with the appropriate structure
    @unpack alphaAppendix, gammaAppendix, objfunConciseDescription = data
    filename = joinpath(base_dir, "validationDecisionVariables_for_$(objfunConciseDescription)_via_$(simNatureAppendix)_alpha_$(alphaAppendix)_gamma_$(gammaAppendix)_with_$(linearizedModelAppendix).txt")

    # Write the valdVals dictionary to a CSV file
    CSV.write(filename, valdVals)

    # Print confirmation if verbose is enabled
    if verbose
        println("Validation decision variables written to $filename")
    end
end
#endregion

#region export_validation_key_results
"""
    export_validation_key_results(modelDict; verbose::Bool=false)

Export key results obtained from validating the Optimization results against OpenDSS powerflow results.
"""
function export_validation_key_results(modelDict; filename::String="validation_results.txt", verbose::Bool=false,
    printEveryTimeStepPowerflow::Bool=true)
    @unpack valdVals, data = modelDict
    # Define the path and filename based on the specified structure
    @unpack T, systemName, numAreas, gedAppendix, machine_ID, simNatureAppendix, linearizedModelAppendix = data
    base_dir = joinpath("processedData", systemName, gedAppendix, "Horizon_$(T)", "numAreas_$(numAreas)")

    if !isdir(base_dir)
        if verbose
            println("Creating directory: $base_dir")
        end
        mkpath(base_dir)
    end

    @unpack solver, alphaAppendix, gammaAppendix, objfunConciseDescription = data;
    # Adjust filename based on printEveryTimeStepPowerflow flag
    if printEveryTimeStepPowerflow
        filename = joinpath(base_dir, "valdResults_for_$(objfunConciseDescription)_via_$(simNatureAppendix)_with_$(linearizedModelAppendix)_alpha_$(alphaAppendix)_gamma_$(gammaAppendix)_full.txt")
    else
        filename = joinpath(base_dir, "valdResults_for_$(objfunConciseDescription)_via_$(simNatureAppendix)_with_$(linearizedModelAppendix)_alpha_$(alphaAppendix)_gamma_$(gammaAppendix).txt")
    end

    # Open the file and write each section
    open(filename, "w") do f
        # Initialize output item counter
        item_counter = 1

        # Header Section
        println(f, "---------------------------------------------")
        println(f, "$(item_counter). Machine ID: $(data[:machine_ID])")
        item_counter += 1
        println(f, "$(item_counter). Solver Used: $(data[:solver])")
        item_counter += 1
        println(f, "$(item_counter). System Name: $(data[:systemName])")
        item_counter += 1
        println(f, "$(item_counter). Horizon Duration: $(data[:T])")
        item_counter += 1
        println(f, "$(item_counter). Nature of Validation Simulation: $(data[:simNatureString])")
        item_counter += 1
        println(f, "$(item_counter). Optimization Model: $(data[:linearizedModelString])")
        item_counter += 1  # Placeholder
        println(f, "---------------------------------------------")

        # Horizon Results
        println(f, "Full $(data[:T]) Hour Horizon Validation Results")
        println(f, "$(item_counter). Horizon Total Substation Power Cost: \$$(round(valdVals[:vald_PSubsCost_allT_dollar], digits=2))")
        item_counter += 1
        println(f, "$(item_counter). Horizon Total Line Loss: $(round(valdVals[:vald_PLoss_allT_kW], digits=2)) kW")
        item_counter += 1
        println(f, "$(item_counter). Horizon Total Substation Power: $(round(valdVals[:vald_PSubs_allT_kW], digits=2)) kW + $(round(valdVals[:vald_QSubs_allT_kVAr], digits=2)) kVAr")
        item_counter += 1
        println(f, "$(item_counter). Horizon Total Load: $(round(valdVals[:vald_load_real_power_allT_kW], digits=2)) kW + $(round(valdVals[:vald_load_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1
        println(f, "$(item_counter). Horizon Total Generation: $(round(valdVals[:vald_total_gen_real_power_allT_kW], digits=2)) kW + $(round(valdVals[:vald_total_gen_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1
        println(f, "$(item_counter). Horizon Total Static Capacitor Reactive Power Generation: $(round(valdVals[:vald_static_cap_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1
        println(f, "$(item_counter). Horizon Total PV Generation: $(round(valdVals[:vald_pv_real_power_allT_kW], digits=2)) kW + $(round(valdVals[:vald_pv_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1
        println(f, "$(item_counter). Horizon Total Battery Generation: $(round(valdVals[:vald_battery_real_power_allT_kW], digits=2)) kW + $(round(valdVals[:vald_battery_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1
        println(f, "$(item_counter). Horizon Total Battery Transaction Magnitude: $(round(valdVals[:vald_battery_real_power_transaction_magnitude_allT_kW], digits=2)) kW + $(round(valdVals[:vald_battery_reactive_power_transaction_magnitude_allT_kVAr], digits=2)) kVAr")
        item_counter += 1
        println(f, "$(item_counter). Horizon Total SCD Observed: N/A")
        item_counter += 1
        println(f, "$(item_counter). Horizon-end Battery Energy Deviation from Reference: $(round(valdVals[:vald_terminal_soc_violation_kWh], digits=2)) kWh")
        item_counter += 1
        println(f, "$(item_counter). Horizon-Total All Time Substation Power Peak: $(round(valdVals[:vald_substation_real_power_peak_allT_kW], digits=2)) kW")
        item_counter += 1

        # Discrepancies
        println(f, "---------------------------------------------")
        println(f, "Discrepancies (Maximum All Time):")
        println(f, "$(item_counter). Maximum All Time Voltage Discrepancy: $(round(valdVals[:disc_voltage_all_time_pu], digits=6)) pu")
        item_counter += 1
        println(f, "$(item_counter). Maximum All Time Line Loss Discrepancy: $(round(valdVals[:disc_line_loss_all_time_kW], digits=6)) kW")
        item_counter += 1
        println(f, "$(item_counter). Maximum All Time Substation Borrowed Real Power Discrepancy: $(round(valdVals[:disc_PSubs_all_time_kW], digits=6)) kW")
        item_counter += 1
        println(f, "$(item_counter). Maximum All Time Substation Borrowed Reactive Power Discrepancy: $(round(valdVals[:disc_QSubs_all_time_kVAr], digits=6)) kVAr")
        item_counter += 1

        # Additional Metadata
        println(f, "---------------------------------------------")
        println(f, "$(item_counter). Solution Time: Small")
        item_counter += 1

        # Optional: Print per-timestep power flow results if enabled
        if printEveryTimeStepPowerflow
            println(f, "\nPer-Timestep Power Flow Results:")
            for t in data[:Tset]
                println(f, "\n" * "*"^30)
                println(f, "   Time Step: $t")
                println(f, "*"^30)
                println(f, "   Power Loss              : $(valdVals[:vald_PLoss_vs_t_1toT_kW][t]) kW")
                println(f, "   Substation Power        : $(valdVals[:vald_PSubs_vs_t_1toT_kW][t]) kW")
                println(f, "   Reactive Power          : $(valdVals[:vald_QSubs_vs_t_1toT_kVAr][t]) kVAr")
                println(f, "   Total Load Power        : $(valdVals[:vald_load_real_power_vs_t_1toT_kW][t]) kW")
                println(f, "   Total Load Reactive Power: $(valdVals[:vald_load_reactive_power_vs_t_1toT_kVAr][t]) kVAr")
                println(f, "   Total PV Power          : $(valdVals[:vald_pv_real_power_vs_t_1toT_kW][t]) kW")
                println(f, "   Total PV Reactive Power : $(valdVals[:vald_pv_reactive_power_vs_t_1toT_kVAr][t]) kVAr")
                println(f, "   Total Battery Power     : $(valdVals[:vald_battery_real_power_vs_t_1toT_kW][t]) kW")
                println(f, "   Total Battery Reactive Power: $(valdVals[:vald_battery_reactive_power_vs_t_1toT_kVAr][t]) kVAr")
                println(f, "*"^30 * "\n")
            end
        end
    end

    if verbose
        println("Validation key results exported to $filename")
    end
end
#endregion

end  # module Exporter
