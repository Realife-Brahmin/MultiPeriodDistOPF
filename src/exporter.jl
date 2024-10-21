module Exporter

using XLSX
import JuMP: value

export export_decision_variables

function export_decision_variables(model, data, filename::String="decision_variables.xlsx")

    Tset = data[:Tset]  # Set of time intervals
    Lset = data[:Lset]  # Set of branch pairs (i, j)
    Nset = data[:Nset]  # Set of buses
    Bset = data[:Bset]  # Set of battery buses
    Dset = data[:Dset]  # Set of PV buses

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

    # Open an Excel file to write
    XLSX.openxlsx(filename, mode="w") do xf
        sheet = xf[1]

        # First row: P_Subs for all time intervals
        row_index = 1
        sheet[row_index, 1] = "P_Subs"
        for t in Tset
            sheet[row_index, t+1] = value(P_Subs[t])
        end
        row_index += 1

        # Next row: P_ij for all (i, j), with all time steps in one row
        for (i, j) in Lset
            sheet[row_index, 1] = "P_ij_$(i)_$(j)"
            for t in Tset
                sheet[row_index, t+1] = value(P[(i, j), t])
            end
            row_index += 1
        end

        # Next row: Q_ij for all (i, j), with all time steps in one row
        for (i, j) in Lset
            sheet[row_index, 1] = "Q_ij_$(i)_$(j)"
            for t in Tset
                sheet[row_index, t+1] = value(Q[(i, j), t])
            end
            row_index += 1
        end

        # Next row: l_ij for all (i, j), with all time steps in one row
        for (i, j) in Lset
            sheet[row_index, 1] = "l_ij_$(i)_$(j)"
            for t in Tset
                sheet[row_index, t+1] = value(l[(i, j), t])
            end
            row_index += 1
        end

        # Next row: v_j for all buses, with all time steps in one row
        for j in Nset
            sheet[row_index, 1] = "v_j_$(j)"
            for t in Tset
                sheet[row_index, t+1] = value(v[j, t])
            end
            row_index += 1
        end

        # Next row: q_D_j for all PV buses, with all time steps in one row
        for j in Dset
            sheet[row_index, 1] = "q_D_j_$(j)"
            for t in Tset
                sheet[row_index, t+1] = value(q_D[j, t])
            end
            row_index += 1
        end

        # Next row: q_B_j for all battery buses, with all time steps in one row
        for j in Bset
            sheet[row_index, 1] = "q_B_j_$(j)"
            for t in Tset
                sheet[row_index, t+1] = value(q_B[j, t])
            end
            row_index += 1
        end

        # Next row: P_c_j for all battery buses, with all time steps in one row
        for j in Bset
            sheet[row_index, 1] = "P_c_j_$(j)"
            for t in Tset
                sheet[row_index, t+1] = value(P_c[j, t])
            end
            row_index += 1
        end

        # Next row: P_d_j for all battery buses, with all time steps in one row
        for j in Bset
            sheet[row_index, 1] = "P_d_j_$(j)"
            for t in Tset
                sheet[row_index, t+1] = value(P_d[j, t])
            end
            row_index += 1
        end

        # Next row: B_j for all battery buses, with all time steps in one row
        for j in Bset
            sheet[row_index, 1] = "B_j_$(j)"
            for t in Tset
                sheet[row_index, t+1] = value(B[j, t])
            end
            row_index += 1
        end

    end

    println("Decision variables exported to $filename")
end

end  # module Exporter
