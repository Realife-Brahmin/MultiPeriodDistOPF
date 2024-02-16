# https://discourse.julialang.org/t/how-to-append-a-single-row-to-an-excel-file-using-xlsx/56203/8
# function append_xl_row(workbook_path::String, sheet_name::String, row_data::Vector)
function append_xl_row(;workbook_path::String,
                        sheet_name::String   ,
                        row_data::Vector     ,
                        Append::Bool=true    )

# println("Writing to Excel: Append = ", string(Append),";\tData = ", string(collect(row_data[[1,6:end]])))
# println("Writing to Excel: Append = ", string(Append),";\tIter# = $(row_data[1]), Data = $(row_data[6:end])")
println("Writing to Excel Log: Append = ", string(Append),";\tIter#$(row_data[1]),\tData : ", join(row_data[6:end],", "))

XLSX.openxlsx(workbook_path, mode="rw") do xf
        sheet = xf[sheet_name]
        num_rows = XLSX.get_dimension(sheet).stop.row_number
        # println("Excel sheet has ",string(num_rows)," rows")

        # if num_rows == 1
        #     sheet[2,1] = row_data
        # else
        sheet[num_rows+Append,:] = row_data
        # end
    end
# println("Wrote to Excel")
end
