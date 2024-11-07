module Exporter

using XLSX
import JuMP: value
include("./helperFunctions.jl")
import .helperFunctions: myprintln  # Import myprintln from the helperFunctions module

using Parameters: @unpack

export export_decision_variables, export_optimization_model, export_simulation_key_results_txt

using DelimitedFiles  # To write CSV files

function export_optimization_model(model, data;
    verbose::Bool=false)
    # Define the path and filename based on the specified structure
    @unpack T, systemName, numAreas, gedAppendix, machine_ID, objfunAppendix, simNatureAppendix = data
    base_dir = joinpath("processedData", systemName, gedAppendix, "Horizon_$(T)", "numAreas_$(numAreas)")

    # Create the directory if it doesn't exist
    if !isdir(base_dir)
        println("Creating directory: $base_dir")
        mkpath(base_dir)
    end

    # Define the filename with the appropriate structure
    filename = joinpath(base_dir, "Horizon_$(T)_$(machine_ID)_optimizationModel_$(gedAppendix)_for_$(objfunAppendix)_via_$(simNatureAppendix).txt")

    # Check if the file already exists, and delete it if so
    if isfile(filename)
        rm(filename)
    end

    # Open a new file and write the model contents to it
    open(filename, "w") do f
        print(f, model)
    end

    myprintln(verbose, "Model successfully written to $filename")
end

function export_decision_variables(model, data;
    filename::String="decision_variables.csv",
    verbose::Bool=false)

    # Define the path and filename based on the specified structure
    @unpack T, systemName, numAreas, gedAppendix, machine_ID, objfunAppendix, simNatureAppendix = data
    base_dir = joinpath("processedData", systemName, gedAppendix, "Horizon_$(T)", "numAreas_$(numAreas)")

    if !isdir(base_dir)
        println("Creating directory: $base_dir")
        mkpath(base_dir)
    end

    ext = ".csv"
    filename = joinpath(base_dir, "Horizon_$(T)_$(machine_ID)_decisionVariables_$(gedAppendix)_for_$(objfunAppendix)_via_$(simNatureAppendix)"*ext)
    
    # Log current working directory
    myprintln(verbose, "Current working directory: $(pwd())")

    # Ensure the filename has a valid path or is saved in the current directory
    myprintln(verbose, "Saving to filename: $filename")

    # Unpack the necessary data from `data`
    Tset = sort(collect(data[:Tset]))  # Time steps
    LoadShapeLoad = data[:LoadShapeLoad]  # Load shape values
    LoadShapePV = data[:LoadShapePV]  # Irradiance values
    LoadShapeCost = data[:LoadShapeCost] * 100  # Convert from $/kWh to cents/kWh

    Lset = sort(collect(data[:Lset]))  # Line set
    Nset = sort(collect(data[:Nset]))  # Node set
    Bset = sort(collect(data[:Bset]))  # Battery set
    Dset = sort(collect(data[:Dset]))  # PV set

    P_Subs = model[:P_Subs]
    P = model[:P]
    Q = model[:Q]
    l = model[:l]
    v = model[:v]
    q_D = model[:q_D]
    q_B = model[:q_B]
    P_c = model[:P_c]
    P_d = model[:P_d]
    B = model[:B]

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
    push!(data_matrix, ["P_Subs"; [value(P_Subs[t]) for t in Tset]...])

    # Next rows: P_ij for all (i, j) pairs in Lset
    for (i, j) in Lset
        push!(data_matrix, ["P_ij_$(i)_$(j)"; [value(P[(i, j), t]) for t in Tset]...])
    end

    # Next rows: Q_ij for all (i, j) pairs in Lset
    for (i, j) in Lset
        push!(data_matrix, ["Q_ij_$(i)_$(j)"; [value(Q[(i, j), t]) for t in Tset]...])
    end

    # Next rows: l_ij for all (i, j) pairs in Lset
    for (i, j) in Lset
        push!(data_matrix, ["l_ij_$(i)_$(j)"; [value(l[(i, j), t]) for t in Tset]...])
    end

    # Next rows: v_j for all nodes in Nset
    for j in Nset
        push!(data_matrix, ["v_j_$(j)"; [value(v[j, t]) for t in Tset]...])
    end

    # Next rows: q_D_j for all PV buses in Dset
    for j in Dset
        push!(data_matrix, ["q_D_j_$(j)"; [value(q_D[j, t]) for t in Tset]...])
    end

    # Next rows: q_B_j for all battery buses in Bset
    for j in Bset
        push!(data_matrix, ["q_B_j_$(j)"; [value(q_B[j, t]) for t in Tset]...])
    end

    # Next rows: P_c_j for all battery buses in Bset
    for j in Bset
        push!(data_matrix, ["P_c_j_$(j)"; [value(P_c[j, t]) for t in Tset]...])
    end

    # Next rows: P_d_j for all battery buses in Bset
    for j in Bset
        push!(data_matrix, ["P_d_j_$(j)"; [value(P_d[j, t]) for t in Tset]...])
    end

    # Next rows: B_j for all battery buses in Bset
    for j in Bset
        push!(data_matrix, ["B_j_$(j)"; [value(B[j, t]) for t in Tset]...])
    end

    # Write the matrix to the CSV file
    try
        myprintln(verbose, "Opening file: $filename")
        writedlm(filename, data_matrix, ',')
        myprintln(verbose, "Data written successfully.")
    catch e
        myprintln(verbose, "Error occurred: $e")
    end

    # Check if file was created
    if isfile(filename)
        myprintln(verbose, "File exists: $filename")
    else
        myprintln(verbose, "File does not exist: $filename")
    end

    myprintln(verbose, "Decision variables exported to $filename")
end


function export_simulation_key_results_txt(model, data; filename::String="simulation_results.txt", verbose::Bool=false)

    # Define the path and filename based on the specified structure
    @unpack T, systemName, numAreas, gedAppendix, machine_ID, objfunAppendix, simNatureAppendix = data
    base_dir = joinpath("processedData", systemName, gedAppendix, "Horizon_$(T)", "numAreas_$(numAreas)")

    if !isdir(base_dir)
        myprintln(verbose, "Creating directory: $base_dir")
        mkpath(base_dir)
    end

    filename = joinpath(base_dir, "Horizon_$(T)_$(machine_ID)_results_$(gedAppendix)_for_$(objfunAppendix)_via_$(simNatureAppendix).txt")

    # Extract system information and parameters from `data`
    # system_name = data[:machine_ID]  # System name
    # horizon_duration = data[:T]  # Horizon duration
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
        println(f, "$(item_counter). Horizon Duration: $(data[:T])")
        item_counter += 1
        println(f, "$(item_counter). Nature of Simulation: $(data[:simNatureString])")
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

        # Example substation power cost
        println(f, "$(item_counter). Horizon Total Substation Power Cost: \$$(round(data[:PSubsCost_allT_dollar], digits=2))")
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
        println(f, "$(item_counter). Simulation Time: $(round(solution_time, digits=2)) s")
        item_counter += 1
        println(f, "$(item_counter). Time to solve with sequential (non-parallel) computation: $(round(solution_time, digits=2)) s")
        item_counter += 1  # Placeholder for non-parallel time
        println(f, "$(item_counter). Time to solve if OPF computation parallelized: $(round(solution_time, digits=2)) s")  # Placeholder for parallelized time
        println(f, "---------------------------------------------")
    end

    myprintln(verbose, "Simulation key results exported to $filename")
end

function export_key_validation_results(vald, data; filename::String="validation_results.txt", verbose::Bool=false)
    # Define the path and filename based on the specified structure
    @unpack T, systemName, numAreas, gedAppendix, machine_ID, simNatureAppendix = data
    base_dir = joinpath("processedData", systemName, gedAppendix, "Horizon_$(T)", "numAreas_$(numAreas)")

    if !isdir(base_dir)
        if verbose
            println("Creating directory: $base_dir")
        end
        mkpath(base_dir)
    end

    filename = joinpath(base_dir, "Horizon_$(T)_$(machine_ID)_valdResults_$(gedAppendix)_via_$(simNatureAppendix).txt")

    # Open the file and write each section
    open(filename, "w") do f
        # Initialize output item counter
        item_counter = 1

        # Header Section
        println(f, "---------------------------------------------")
        println(f, "$(item_counter). Machine ID: $(data[:machine_ID])")
        item_counter += 1
        println(f, "$(item_counter). Horizon Duration: $(data[:T])")
        item_counter += 1
        println(f, "$(item_counter). Nature of Simulation: Validation")
        item_counter += 1
        println(f, "---------------------------------------------")

        # Horizon Results
        println(f, "Full $(data[:T]) Hour Horizon Validation Results")

        # Line Loss
        println(f, "$(item_counter). Horizon Total Line Loss: $(round(vald[:vald_PLoss_allT_kW], digits=2)) kW")
        item_counter += 1

        # Substation Power
        println(f, "$(item_counter). Horizon Total Substation Power: $(round(vald[:vald_PSubs_allT_kW], digits=2)) kW + $(round(vald[:vald_QSubs_allT_kVAr], digits=2)) kVAr")
        item_counter += 1

        # Load and Generation
        println(f, "$(item_counter). Horizon Total Load: $(round(vald[:vald_load_real_power_allT_kW], digits=2)) kW + $(round(vald[:vald_load_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1
        println(f, "$(item_counter). Horizon Total Generation: $(round(vald[:vald_total_gen_real_power_allT_kW], digits=2)) kW + $(round(vald[:vald_total_gen_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1

        # Static Capacitor Power
        println(f, "$(item_counter). Horizon Total Static Capacitor Reactive Power Generation: $(round(vald[:vald_static_cap_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1

        # Substation Power Cost
        println(f, "$(item_counter). Horizon Total Substation Power Cost: \$$(round(vald[:vald_PSubsCost_allT_dollar], digits=2))")
        item_counter += 1

        # PV Generation
        println(f, "$(item_counter). Horizon Total PV Generation: $(round(vald[:vald_pv_real_power_allT_kW], digits=2)) kW + $(round(vald[:vald_pv_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1

        # Battery Generation and Transactions
        println(f, "$(item_counter). Horizon Total Battery Generation: $(round(vald[:vald_battery_real_power_allT_kW], digits=2)) kW + $(round(vald[:vald_battery_reactive_power_allT_kVAr], digits=2)) kVAr")
        item_counter += 1
        println(f, "$(item_counter). Horizon Total Battery Transaction Magnitude: $(round(vald[:vald_battery_real_power_transaction_magnitude_allT_kW], digits=2)) kW + $(round(vald[:vald_battery_reactive_power_transaction_magnitude_allT_kVAr], digits=2)) kVAr")
        item_counter += 1

        # SCD and Energy Deviation
        println(f, "$(item_counter). Horizon Total SCD Observed: N/A")
        item_counter += 1
        println(f, "$(item_counter). Horizon-end Battery Energy Deviation from Reference: $(round(vald[:vald_terminal_soc_violation_kWh], digits=2)) kWh")
        item_counter += 1

        # Peak Power
        println(f, "$(item_counter). Horizon-Total All Time Substation Power Peak: $(round(vald[:vald_substation_real_power_peak_allT_kW], digits=2)) kW")
        item_counter += 1

        # Discrepancies
        println(f, "---------------------------------------------")
        println(f, "Discrepancies (Maximum All Time):")
        println(f, "$(item_counter). Maximum All Time Voltage Discrepancy: $(round(vald[:disc_voltage_all_time], digits=6)) pu")
        item_counter += 1
        println(f, "$(item_counter). Maximum All Time Line Loss Discrepancy: $(round(vald[:disc_line_loss_all_time], digits=6)) kW")
        item_counter += 1
        println(f, "$(item_counter). Maximum All Time Substation Borrowed Real Power Discrepancy: $(round(vald[:disc_PSubs_all_time], digits=6)) kW")
        item_counter += 1
        println(f, "$(item_counter). Maximum All Time Substation Borrowed Reactive Power Discrepancy: $(round(vald[:disc_QSubs_all_time], digits=6)) kVAr")
        item_counter += 1

        # Additional Metadata
        println(f, "---------------------------------------------")
        println(f, "$(item_counter). Solution Time: $(round(vald[:solution_time], digits=2)) s")
        item_counter += 1

    end

    if verbose
        println("Validation key results exported to $filename")
    end
end

end  # module Exporter
