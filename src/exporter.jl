module Exporter

using XLSX
import JuMP: value
include("./helperFunctions.jl")
import .helperFunctions: myprintln  # Import myprintln from the helperFunctions module

using Parameters: @unpack

export export_decision_variables

function export_decision_variables(model, data;
    filename::String="decision_variables.xlsx",
    verbose::Bool=false)

    # Log current working directory
    myprintln(verbose, "Current working directory: $(pwd())")

    # Ensure the filename has a valid path or is saved in the current directory
    myprintln(verbose, "Saving to filename: $filename")

    # Unpack the data
    # @unpack Tset, Lset, Nset, Bset, Dset = data
    Tset = sort(collect(data[:Tset]))  # Convert Set to array and sort
    Lset = sort(collect(data[:Lset]))  # Convert Set to array and sort
    Nset = sort(collect(data[:Nset]))  # Convert Set to array and sort
    Bset = sort(collect(data[:Bset]))  # Convert Set to array and sort
    Dset = sort(collect(data[:Dset]))  # Convert Set to array and sort

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

            # First row: P_Subs for all time intervals, sorted by t
            row_index = 1
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

end  # module Exporter
