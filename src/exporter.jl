module Exporter

using XLSX
import JuMP: value
include("./helperFunctions.jl")
import .helperFunctions: myprintln  # Import myprintln from the helperFunctions module

using Parameters: @unpack

export export_decision_variables, export_simulation_key_results_txt

function export_decision_variables(model, data;
    filename::String="decision_variables.xlsx",
    verbose::Bool=false)

    # Log current working directory
    myprintln(verbose, "Current working directory: $(pwd())")

    # Ensure the filename has a valid path or is saved in the current directory
    myprintln(verbose, "Saving to filename: $filename")

    # Unpack the necessary data from `data`
    Tset = sort(collect(data[:Tset]))  # Time steps
    LoadShape = data[:LoadShape]  # Load shape values
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

    # Use the `do` block syntax to ensure that the file is closed automatically
    try
        myprintln(verbose, "Opening file: $filename")
        XLSX.openxlsx(filename, mode="w") do xf
            myprintln(verbose, "File opened successfully.")
            sheet = xf[1]

            # Add the headers before the decision variables

            # First header: "t"
            row_index = 1
            sheet[row_index, 1] = "t"
            for t in Tset
                sheet[row_index, t+1] = t
            end
            row_index += 1

            # Second header: LoadShape values (lambda)
            sheet[row_index, 1] = "lambda"
            for t in eachindex(LoadShape)
                sheet[row_index, t+1] = LoadShape[t]
            end
            row_index += 1

            # Third header: Irradiance (LoadShapePV)
            sheet[row_index, 1] = "Irrad"
            for t in eachindex(LoadShapePV)
                sheet[row_index, t+1] = LoadShapePV[t]
            end
            row_index += 1

            # Fourth header: cents/kWh (LoadShapeCost)
            sheet[row_index, 1] = "cents/kWh"
            for t in eachindex(LoadShapeCost)
                sheet[row_index, t+1] = LoadShapeCost[t]
            end
            row_index += 1

            # Now start the original decision variables

            # First row: P_Subs for all time intervals, sorted by t
            sheet[row_index, 1] = "P_Subs"
            for t in Tset
                sheet[row_index, t+1] = value(P_Subs[t])
            end
            row_index += 1

            # Next row: P_ij for all (i, j), with all time steps in one row, sorted by (i, j)
            for (i, j) in Lset
                sheet[row_index, 1] = "P_ij_$(i)_$(j)"
                for t in Tset
                    sheet[row_index, t+1] = value(P[(i, j), t])
                end
                row_index += 1
            end

            # Next row: Q_ij for all (i, j), with all time steps in one row, sorted by (i, j)
            for (i, j) in Lset
                sheet[row_index, 1] = "Q_ij_$(i)_$(j)"
                for t in Tset
                    sheet[row_index, t+1] = value(Q[(i, j), t])
                end
                row_index += 1
            end

            # Next row: l_ij for all (i, j), with all time steps in one row, sorted by (i, j)
            for (i, j) in Lset
                sheet[row_index, 1] = "l_ij_$(i)_$(j)"
                for t in Tset
                    sheet[row_index, t+1] = value(l[(i, j), t])
                end
                row_index += 1
            end

            # Next row: v_j for all buses, with all time steps in one row, sorted by j
            for j in Nset
                sheet[row_index, 1] = "v_j_$(j)"
                for t in Tset
                    sheet[row_index, t+1] = value(v[j, t])
                end
                row_index += 1
            end

            # Next row: q_D_j for all PV buses, with all time steps in one row, sorted by j
            for j in Dset
                sheet[row_index, 1] = "q_D_j_$(j)"
                for t in Tset
                    sheet[row_index, t+1] = value(q_D[j, t])
                end
                row_index += 1
            end

            # Next row: q_B_j for all battery buses, with all time steps in one row, sorted by j
            for j in Bset
                sheet[row_index, 1] = "q_B_j_$(j)"
                for t in Tset
                    sheet[row_index, t+1] = value(q_B[j, t])
                end
                row_index += 1
            end

            # Next row: P_c_j for all battery buses, with all time steps in one row, sorted by j
            for j in Bset
                sheet[row_index, 1] = "P_c_j_$(j)"
                for t in Tset
                    sheet[row_index, t+1] = value(P_c[j, t])
                end
                row_index += 1
            end

            # Next row: P_d_j for all battery buses, with all time steps in one row, sorted by j
            for j in Bset
                sheet[row_index, 1] = "P_d_j_$(j)"
                for t in Tset
                    sheet[row_index, t+1] = value(P_d[j, t])
                end
                row_index += 1
            end

            # Next row: B_j for all battery buses, with all time steps in one row, sorted by j
            for j in Bset
                sheet[row_index, 1] = "B_j_$(j)"
                for t in Tset
                    sheet[row_index, t+1] = value(B[j, t])
                end
                row_index += 1
            end

            myprintln(verbose, "Data written successfully.")
        end

        # Check if file was created
        if isfile(filename)
            myprintln(verbose, "File exists: $filename")
        else
            myprintln(verbose, "File does not exist: $filename")
        end

    catch e
        myprintln(verbose, "Error occurred: $e")
    end

    myprintln(verbose, "Decision variables exported to $filename")
end

function export_simulation_key_results_txt(model, data; filename::String="simulation_results.txt", verbose::Bool=false)

    # Define the path and filename based on the specified structure
    @unpack T, systemName, numAreas, gedAppendix, machine_ID, objfunAppendix, simNatureAppendix = data
    base_dir = joinpath("processedData", systemName, "numAreas_$(numAreas)", "Horizon_$(T)", gedAppendix)

    if !isdir(base_dir)
        println("Creating directory: $base_dir")
        mkpath(base_dir)
    end

    filename = joinpath(base_dir, "Horizon_$(T)_$(machine_ID)_results_$(gedAppendix)_for_$(objfunAppendix)_via_$(simNatureAppendix).txt")

    # Extract system information and parameters from `data`
    # system_name = data[:machine_ID]  # System name
    # horizon_duration = data[:T]  # Horizon duration
    macroItrsCompleted = get(data, :macroItrsCompleted, 1)  # Default to 1 if not set
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

        # Todo: Add PV outputs too
        
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
        println(f, "$(item_counter). Number of Macro-Iterations: $(macroItrsCompleted)")
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

end  # module Exporter
