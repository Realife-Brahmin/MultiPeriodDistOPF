module Exporter

using XLSX
import JuMP: value

export export_decision_variables

function export_decision_variables(model, data, filename::String="decision_variables.xlsx")

    Tset_sorted = sort(collect(data[:Tset]))  # Convert Set to array and sort
    Lset_sorted = sort(collect(data[:Lset]))  # Convert Set to array and sort
    Nset_sorted = sort(collect(data[:Nset]))  # Convert Set to array and sort
    Bset_sorted = sort(collect(data[:Bset]))  # Convert Set to array and sort
    Dset_sorted = sort(collect(data[:Dset]))  # Convert Set to array and sort

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

    xf = XLSX.openxlsx(filename, mode="w")  # Open the file explicitly

    try
        sheet = xf[1]

        # First row: P_Subs for all time intervals, sorted by t
        row_index = 1
        sheet[row_index, 1] = "P_Subs"
        for t in Tset_sorted
            sheet[row_index, t+1] = value(P_Subs[t])
        end
        row_index += 1

        # Next row: P_ij for all (i, j), with all time steps in one row, sorted by (i, j)
        for (i, j) in Lset_sorted
            sheet[row_index, 1] = "P_ij_$(i)_$(j)"
            for t in Tset_sorted
                sheet[row_index, t+1] = value(P[(i, j), t])
            end
            row_index += 1
        end

        # Next row: Q_ij for all (i, j), with all time steps in one row, sorted by (i, j)
        for (i, j) in Lset_sorted
            sheet[row_index, 1] = "Q_ij_$(i)_$(j)"
            for t in Tset_sorted
                sheet[row_index, t+1] = value(Q[(i, j), t])
            end
            row_index += 1
        end

        # Next row: l_ij for all (i, j), with all time steps in one row, sorted by (i, j)
        for (i, j) in Lset_sorted
            sheet[row_index, 1] = "l_ij_$(i)_$(j)"
            for t in Tset_sorted
                sheet[row_index, t+1] = value(l[(i, j), t])
            end
            row_index += 1
        end

        # Next row: v_j for all buses, with all time steps in one row, sorted by j
        for j in Nset_sorted
            sheet[row_index, 1] = "v_j_$(j)"
            for t in Tset_sorted
                sheet[row_index, t+1] = value(v[j, t])
            end
            row_index += 1
        end

        # Next row: q_D_j for all PV buses, with all time steps in one row, sorted by j
        for j in Dset_sorted
            sheet[row_index, 1] = "q_D_j_$(j)"
            for t in Tset_sorted
                sheet[row_index, t+1] = value(q_D[j, t])
            end
            row_index += 1
        end

        # Next row: q_B_j for all battery buses, with all time steps in one row, sorted by j
        for j in Bset_sorted
            sheet[row_index, 1] = "q_B_j_$(j)"
            for t in Tset_sorted
                sheet[row_index, t+1] = value(q_B[j, t])
            end
            row_index += 1
        end

        # Next row: P_c_j for all battery buses, with all time steps in one row, sorted by j
        for j in Bset_sorted
            sheet[row_index, 1] = "P_c_j_$(j)"
            for t in Tset_sorted
                sheet[row_index, t+1] = value(P_c[j, t])
            end
            row_index += 1
        end

        # Next row: P_d_j for all battery buses, with all time steps in one row, sorted by j
        for j in Bset_sorted
            sheet[row_index, 1] = "P_d_j_$(j)"
            for t in Tset_sorted
                sheet[row_index, t+1] = value(P_d[j, t])
            end
            row_index += 1
        end

        # Next row: B_j for all battery buses, with all time steps in one row, sorted by j
        for j in Bset_sorted
            sheet[row_index, 1] = "B_j_$(j)"
            for t in Tset_sorted
                sheet[row_index, t+1] = value(B[j, t])
            end
            row_index += 1
        end
    finally
        XLSX.close(xf)  # Ensure the file is closed
    end

    println("Decision variables exported to $filename")
end

end  # module Exporter
