## -------------- This code reads the sceanrios file ---------------------------
##----------- and store the scenarios for both PV and Wind Power ---------------

function interface_scenarios_fun(filename)
sheetname  = "wind_scenarios";
fields     = ["Scenarios","t1 (MW)","t2 (MW)","t3 (MW)","t4 (MW)","t5 (MW)","t6 (MW)","t7 (MW)","t8 (MW)","t9 (MW)","t10 (MW)","t11 (MW)","t12 (MW)","t13 (MW)","t14 (MW)","t15 (MW)","t16 (MW)","t17 (MW)","t18 (MW)","t19 (MW)","t20 (MW)","t21 (MW)","t22 (MW)","t23 (MW)","t24 (MW)"]; # Fields that have to be read from the file
raw_data   = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array
header     = raw_data[1,:]
data       = raw_data[2:end,2:size(header,1)]
data       = convert(Array{Float64}, data)
# data_cont  = zeros(size(fields,1))               # data_cont = Data Container
nPprf_wind_sc  = Int64(size(data,1))

nSc = Int64(size(data,1))

rheader_pProfile_wind_sc = header    # Exporting raw header of buses sheet
rdata_pProfile_wind_sc   = data      # Exporting raw data of buses sheet
raw_data = nothing
header   = nothing
data     = nothing


sheetname  = "pv_scenarios";
fields     = ["Scenarios","t1 (MW)","t2 (MW)","t3 (MW)","t4 (MW)","t5 (MW)","t6 (MW)","t7 (MW)","t8 (MW)","t9 (MW)","t10 (MW)","t11 (MW)","t12 (MW)","t13 (MW)","t14 (MW)","t15 (MW)","t16 (MW)","t17 (MW)","t18 (MW)","t19 (MW)","t20 (MW)","t21 (MW)","t22 (MW)","t23 (MW)","t24 (MW)"]; # Fields that have to be read from the file
raw_data   = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array
header     = raw_data[1,:]
data       = raw_data[2:end,2:size(header,1)]
data       = convert(Array{Float64}, data)
# data_cont  = zeros(size(fields,1))               # data_cont = Data Container
nPprf_pv_sc  = Int64(size(data,1))

rheader_pProfile_pv_sc = header    # Exporting raw header of buses sheet
rdata_pProfile_pv_sc   = data      # Exporting raw data of buses sheet
raw_data = nothing
header   = nothing
data     = nothing


sheetname  = "probabilities";
fields     = ["Scenarios","prob"]; # Fields that have to be read from the file
raw_data   = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array
header     = raw_data[1,:]
data       = raw_data[2:end,2:size(header,1)]
data       = convert(Array{Float64}, data)
# data_cont  = zeros(size(fields,1))               # data_cont = Data Container
n_prob  = Int64(size(data,1))

rheader_prob_sc = header    # Exporting raw header of buses sheet
rdata_prob_sc   = data      # Exporting raw data of buses sheet
raw_data = nothing
header   = nothing
data     = nothing

return begin   nPprf_wind_sc, nSc, rheader_pProfile_wind_sc, rdata_pProfile_wind_sc,
               nPprf_pv_sc, rheader_pProfile_pv_sc, rdata_pProfile_pv_sc,
               n_prob, rheader_prob_sc, rdata_prob_sc end
end