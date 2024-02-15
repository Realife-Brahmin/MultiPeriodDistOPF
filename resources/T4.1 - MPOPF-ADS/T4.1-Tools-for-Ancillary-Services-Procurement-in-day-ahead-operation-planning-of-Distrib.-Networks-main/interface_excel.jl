# This code reads each individual sheet of an ODS file and then process the available
# data in order to save the inforamtion related to each bus, line, load, gen etc.

# For this program to work, the number as well type of enteries mentioned in the
# field variable and data_type struct variables MUST be the same.
#---------------------------Formatting the Buses Data---------------------------
##---------------------------- Data from T2.3 File -----------------------------

function interface_excel_fun(filename_mat,filename_addt)
sheetname_mat  = "Buses"
sheetname_addt = "Buses_Addt"
sheetname_oltc_addt = "Buses_OLTC_Addt"
fields_mat   = ["bus_i","type","area","Vm","Va","baseKV","zone","Vmax","Vmin"]; # Fields that have to be read from the file
fields_addt  = ["bus_i","hasGEN","isLOAD","SNOM_MVA","SX","SY","GX","GY"]; # Fields that have to be read from the file
fields       = permutedims(hcat(permutedims(fields_mat),permutedims(fields_addt[2:end])))
################################################################################
##--------------------- Data from Basic MATPOWER File --------------------------
################################################################################
raw_data_mat  = ods_readall(filename_mat;sheetsNames=[sheetname_mat],innerType="Matrix")
raw_data_mat  = raw_data_mat[sheetname_mat]   # Conversion from Dict to Array
header_mat    = raw_data_mat[1,:]
data_mat      = raw_data_mat[2:end,1:size(header_mat,1)]
data_mat      = convert(Array{Float64},data_mat)

################################################################################
##------------------ Data from Additional Felxibility File ---------------------
################################################################################
raw_data_addt = ods_readall(filename_addt;sheetsNames=[sheetname_addt],innerType="Matrix")
raw_data_addt = raw_data_addt[sheetname_addt]   # Conversion from Dict to Array
header_addt   = raw_data_addt[1,:]
data_addt     = raw_data_addt[2:end,1:size(header_addt,1)]
data_addt     = convert(Array{Float64},data_addt)
################################################################################
##------------------ Data from Additional Felxibility File ---------------------
##------------------------- Buses OLTC Data ------------------------------------
################################################################################
raw_data_oltc_addt = ods_readall(filename_addt;sheetsNames=[sheetname_oltc_addt],innerType="Matrix")
raw_data_oltc_addt = raw_data_oltc_addt[sheetname_oltc_addt]   # Conversion from Dict to Array
header_oltc_addt   = raw_data_oltc_addt[1,:]
data_oltc_addt     = raw_data_oltc_addt[2:end,1:size(header_oltc_addt,1)]
data_oltc_addt     = convert(Array{Float64},data_oltc_addt)

# header    = permutedims(hcat(permutedims(header_mat),permutedims(header_addt[2:end])))
# data      = hcat(data_mat,data_addt[:,2:end])

header    = header_mat
data      = data_mat

if !isempty(data_oltc_addt)
    display("No OLTC transfomer found in the original data T2.3")
    display("An OLTC transfomer is added in the network")
    data      = vcat(data_oltc_addt[:,1:9],data)
end

nBus          = Int64(size(data,1))
array_bus     = data_reader(nBus,fields_mat,header,data,Bus)
rheader_buses = header      # Exporting raw header of buses sheet
rdata_buses   = data        # Exporting raw data of buses sheet

raw_data_mat   = nothing
raw_data_addt  = nothing
raw_data_oltc_addt = nothing
header_mat     = nothing
header_addt    = nothing
header_oltc_addt = nothing
header         = nothing
data_mat       = nothing
data_addt      = nothing
data_oltc_addt = nothing
data           = nothing

#----------------------------Formatting the Lines Data--------------------------
sheetname_mat  = "Lines";
sheetname_addt = "Lines_Addt"
sheetname_oltc_addt = "Lines_OLTC_Addt"
fields_mat     = ["fbus","tbus","r","x","b","rateA","rateB","rateC","ratio","angle","status","angmin","angmax"] # Fields that have to be read from the file
fields_addt    = ["fbus","tbus","step_size","actTap","minTap","maxTap","normalTap","nominalRatio","r_ip","r_n","r0","x0","b0","length(meter)","NormSTAT","g","minTapratio","maxTapratio","status_oltc_flex"]
fields         = permutedims(hcat(permutedims(fields_mat),permutedims(fields_addt[3:end])))
################################################################################
##--------------------- Data from Basic MATPOWER File --------------------------
################################################################################
raw_data_mat  = ods_readall(filename_mat;sheetsNames=[sheetname_mat],innerType="Matrix")
raw_data_mat  = raw_data_mat[sheetname_mat]   # Conversion from Dict to Array
header_mat    = raw_data_mat[1,:]
data_mat      = raw_data_mat[2:end,1:size(header_mat,1)]
data_mat      = convert(Array{Float64},data_mat)
################################################################################
##------------------ Data from Additional Felxibility File ---------------------
################################################################################
raw_data_addt = ods_readall(filename_addt;sheetsNames=[sheetname_addt],innerType="Matrix")
raw_data_addt = raw_data_addt[sheetname_addt]   # Conversion from Dict to Array
header_addt   = raw_data_addt[1,:]
data_addt     = raw_data_addt[2:end,1:size(header_addt,1)]

data_addt     = convert(Array{Float64},data_addt)
################################################################################
##------------------ Data from Additional Felxibility File ---------------------
##------------------------- Buses OLTC Data ------------------------------------
################################################################################
raw_data_oltc_addt = ods_readall(filename_addt;sheetsNames=[sheetname_oltc_addt],innerType="Matrix")
raw_data_oltc_addt = raw_data_oltc_addt[sheetname_oltc_addt]   # Conversion from Dict to Array
header_oltc_addt   = raw_data_oltc_addt[1,:]
data_oltc_addt     = raw_data_oltc_addt[2:end,1:size(header_oltc_addt,1)]
data_oltc_addt     = convert(Array{Float64},data_oltc_addt)

header    = permutedims(hcat(permutedims(header_mat),permutedims(header_addt[3:end])))
data      = hcat(data_mat,data_addt[:,3:end])

if !isempty(data_oltc_addt)
    data      = vcat(data_oltc_addt,data)
end
data      = convert(Array{Float64}, data)
#
# if new_line_order == 1 # New_line order sheet is available in the given excel sheet
#     sheetname = "New_Line_Order";
#     fields    = ["From","To"]
#     raw_data_nlo  = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
#     raw_data_nlo  = raw_data_nlo[sheetname]   # Conversion from Dict to Array
#     header_nlo    = raw_data_nlo[1,:]
#     data_nlo      = raw_data_nlo[2:end,:]
#     data_nlo_tf   = hcat(data_nlo[:,2],data_nlo[:,1])
#
#     data_nlo_ft   = convert(Array{Float64},data_nlo)
#     data_nlo_tf   = convert(Array{Float64},data_nlo_tf)
#
#     nLines_nlo    = Int64(size(data_nlo_ft,1))
#     idx_zero = [0.0 0.0]
#     data_new = zeros(size(data))
#     data_new = Array{Any}(nothing,(size(data)))
#     data_new[:,1:2] = convert.(Int64,data_nlo_ft)
#     check_ft = []
#     check_tf = []
#     for i in 1:nLines_nlo
#         confg_nlo_ft = permutedims(data_nlo_ft[i,:]).-data[:,1:2]
#         confg_nlo_tf = permutedims(data_nlo_tf[i,:]).-data[:,1:2]
#
#         idx_line_ft  = findall(all(confg_nlo_ft.==idx_zero,dims=2))
#         idx_line_tf  = findall(all(confg_nlo_tf.==idx_zero,dims=2))
#
#         if !isempty(idx_line_ft) idx_line_ft = idx_line_ft[1][1] end
#         if !isempty(idx_line_tf) idx_line_tf = idx_line_tf[1][1] end
#
#         if !isempty(idx_line_ft) && isempty(idx_line_tf)
#             data_new[i,3:end] = permutedims(data[idx_line_ft,3:end])
#         elseif isempty(idx_line_ft) && !isempty(idx_line_tf)
#             data_new[i,3:end] = permutedims(data[idx_line_tf,3:end])
#         else
#             println("$i,ERROR in LINES DATA")
#         end
#     end
# end
# data = data_new
#------------ Code for handling the lines whose status is zero -----------------
#------- This code will eliminate those lines whose status is set to zero-------
idx_br_status = findall(x->x=="status",header)
open_branches = findall(x->x==0,data[:,idx_br_status])
idx_open_branches = []
for i in 1:size(open_branches,1)
    push!(idx_open_branches,open_branches[i][1])
end
idx_closed_branches = setdiff(collect(1:size(data,1)),idx_open_branches)
data = data[idx_closed_branches,:]
#------------------- Code for handling the parallel lines ----------------------
idx_plines = []                                                                  # Indices of parallel lines
for i in 1:size(data,1)
    line = transpose(data[i,1:2])
    plines = line.-data[:,1:2]
    i0_from = findall(x->x==0.0,plines[:,1])
    i0_to   = findall(x->x==0.0,plines[:,2])
    ilines  = intersect(i0_from,i0_to)
    # ilines  = findall(x->x==0.0,plines[:,1])
    if isempty(idx_plines)
        push!(idx_plines,ilines)
    elseif isempty(findall(x->x==ilines,idx_plines))
        push!(idx_plines,ilines)
    else
        # do nothing!
    end
end
new_line_data = []
for i in 1:size(idx_plines,1)
    ilines = idx_plines[i,1]
    z_eq = sum((data[ilines,3]+data[ilines,4]im).^-1,dims=1).^-1
    y_eq = sum(data[ilines,27]+data[ilines,5]im,dims=1)
    amp_eq_A = sum(data[ilines,6],dims=1)
    amp_eq_B = sum(data[ilines,7],dims=1)
    amp_eq_C = sum(data[ilines,8],dims=1)
    nline_data = [transpose(data[ilines[1,1],1:2]) real(z_eq) imag(z_eq) real(y_eq) imag(y_eq) amp_eq_A amp_eq_B amp_eq_C transpose(data[ilines[1,1],[9:13;28:30]])]
    push!(new_line_data,nline_data)
end
data = vcat(new_line_data...)
data[:,1] = convert.(Int64,data[:,1])
data[:,2] = convert.(Int64,data[:,2])

idx_ratio = findall(x->x!=0,data[:,10])

if !isempty(idx_ratio)
data[idx_ratio,10] = data[idx_ratio,10].^-1
end
header = header[[1:13;27:end],1]
fields = fields[[1:13;27:end],1]

header = permutedims(hcat(permutedims(header[1:4]),(header[14]),permutedims(header[5:13]),permutedims(header[15:17])))
fields = permutedims(hcat(permutedims(fields[1:4]),(fields[14]),permutedims(fields[5:13]),permutedims(fields[15:17])))
##----------------------------------------------------------------------------##
nLines      = Int64(size(data,1))
array_lines = data_reader(nLines,fields,header,data,Line)

rheader_lines = header    # Exporting raw header of lines sheet
rdata_lines   = data        # Exporting raw data of lines sheet

raw_data_mat   = nothing
raw_data_addt  = nothing
raw_data_oltc_addt = nothing
header_mat     = nothing
header_addt    = nothing
header_oltc_addt = nothing
header         = nothing
data_mat       = nothing
data_addt      = nothing
data_oltc_addt = nothing
data           = nothing
#---------------------- Formatting the Transformer Data ---------------------- #
sheetname_addt = "Lines_Addt"
sheetname_oltc_addt = "Lines_OLTC_Addt"
fields_addt    = ["fbus","tbus","step_size","minTap","maxTap","normalTap","nominalRatio","minTapratio","maxTapratio","status_oltc_flex"]
fields         = fields_addt
################################################################################
##------------------ Data from Additional Flexibility File ---------------------
################################################################################
raw_data_addt = ods_readall(filename_addt;sheetsNames=[sheetname_addt],innerType="Matrix")
raw_data_addt = raw_data_addt[sheetname_addt]   # Conversion from Dict to Array
header_addt   = raw_data_addt[1,:]
data_addt     = raw_data_addt[2:end,1:size(header_addt,1)]
data_addt     = convert(Array{Float64},data_addt)
################################################################################
##------------------ Data from Additional Flexibility File ---------------------
##------------------------- Lines OLTC Data ---------------------------------###
################################################################################
raw_data_oltc_addt = ods_readall(filename_addt;sheetsNames=[sheetname_oltc_addt],innerType="Matrix")
raw_data_oltc_addt = raw_data_oltc_addt[sheetname_oltc_addt]   # Conversion from Dict to Array
header_oltc_addt   = raw_data_oltc_addt[1,:]
data_oltc_addt     = raw_data_oltc_addt[2:end,1:size(header_oltc_addt,1)]
data_oltc_addt     = convert(Array{Float64},data_oltc_addt)

data_oltc_addt    =  hcat(data_oltc_addt[:,1:2],data_oltc_addt[:,14:end])

if !isempty(data_oltc_addt)
    data_addt = vcat(data_oltc_addt,data_addt)
end

idx_trsf  = findall(x->x!=0,data_addt[:,8])                                      # Column 8 is Nominal Ratio Column
data_addt = data_addt[idx_trsf,:]                                                # Updated Data contains information about only Transformers

idx_plines = []                                                                  # Indices of parallel lines
for i in 1:size(data_addt,1)
    line = transpose(data_addt[i,1:2])
    plines = line.-data_addt[:,1:2]
    i0_from = findall(x->x==0.0,plines[:,1])
    i0_to   = findall(x->x==0.0,plines[:,2])
    ilines  = intersect(i0_from,i0_to)
    # ilines  = findall(x->x==0.0,plines[:,1])
    if isempty(idx_plines)
        push!(idx_plines,ilines)
    elseif isempty(findall(x->x==ilines,idx_plines))
        push!(idx_plines,ilines)
    else
        # do nothing!
    end
end

new_trsf_data = []
for i in 1:size(idx_plines,1)
    ilines = idx_plines[i,1]
    ntrsf_data = transpose(data_addt[ilines[1,1],:])
    push!(new_trsf_data,ntrsf_data)
end

data_addt = vcat(new_trsf_data...)
nTrsf     = Int64(size(data_addt,1))          # Number of transformers in a network

idx_trsf_s1 = findall(x->x==1,data_addt[:,end])
nTrsf_s1  = size(idx_trsf_s1,1)                   # No of transformers whose status is set to 1 (which has OLTC option)

idx_ratio = findall(x->x!=0,data_addt[:,8])

if !isempty(idx_ratio)
    data_addt[idx_ratio,8] = data_addt[idx_ratio,8].^-1
end

header = header_addt
data = data_addt
array_trsf = data_reader(nTrsf,fields,header,data,Transformer)

rheader_trsf = header    # Exporting raw header of loads sheet
rdata_trsf   = data      # Exporting raw data of loads sheet

raw_data_addt = nothing
header_addt   = nothing
data_addt     = nothing
header        = nothing
data          = nothing
#------------------------ Formatting the Loads Data-----------------------------
sheetname_mat  = "Loads"
sheetname_addt = "Loads_Addt"
fields_mat  = ["bus_i","Pd","Qd","Gs","Bs"];                                         # Fields that have to be read from the file
fields_addt = ["bus_i","Status","cost_inc","cost_dec"]
fields      = permutedims(hcat(permutedims(fields_mat),permutedims(fields_addt[2:end])))
################################################################################
##--------------------- Data from Basic MATPOWER File --------------------------
################################################################################
raw_data_mat  = ods_readall(filename_mat;sheetsNames=[sheetname_mat],innerType="Matrix")
raw_data_mat  = raw_data_mat[sheetname_mat]   # Conversion from Dict to Array
header_mat    = raw_data_mat[1,:]
data_mat      = raw_data_mat[2:end,:]
data_mat      = convert(Array{Float64},data_mat)

################################################################################
##------------------ Data from Additional Felxibility File ---------------------
################################################################################
raw_data_addt = ods_readall(filename_addt;sheetsNames=[sheetname_addt],innerType="Matrix")
raw_data_addt = raw_data_addt[sheetname_addt]   # Conversion from Dict to Array
header_addt   = raw_data_addt[1,:]
data_addt     = raw_data_addt[2:end,1:size(header_addt,1)]
data_addt     = convert(Array{Float64},data_addt)
header    = permutedims(hcat(permutedims(header_mat),permutedims(header_addt[2:end])))
# @infiltrate
data      = hcat(data_mat,data_addt[:,2:end])

nLoads    = Int64(size(data,1))
array_loads = data_reader(nLoads,fields,header,data,Loads)

rheader_loads = header    # Exporting raw header of loads sheet
rdata_loads   = data        # Exporting raw data of loads sheet

raw_data_mat  = nothing
raw_data_addt = nothing
header_mat    = nothing
data_mat      = nothing
header_addt   = nothing
data_addt     = nothing
header        = nothing
data          = nothing
#------------------------ Formatting the Gens Data------------------------------
sheetname_mat        = "Gens"
sheetname_res_addt   = "RES_Addt"
sheetname_curt_addt  = "Gens_Addt"
fields_mat_res_addt  = ["bus_i","Pg","Qg","Qmax","Qmin","Vg","mBase","status","Pmax","Pmin","Pc1","Pc2",
                        "Qc1min","Qc1max","Qc2min","Qc2max","ramp_agc","ramp_10","ramp_30","ramp_q","apf"] # Fields that have to be read from the file
fields_addt = ["bus_i","CurtStatus","Dgtype"]
fields      = permutedims(hcat(permutedims(fields_mat_res_addt),permutedims(fields_addt[2:end])))
################################################################################
##--------------------- Data from Basic MATPOWER File --------------------------
################################################################################
raw_data_mat  = ods_readall(filename_mat;sheetsNames=[sheetname_mat],innerType="Matrix")
raw_data_mat  = raw_data_mat[sheetname_mat]   # Conversion from Dict to Array
header_mat    = raw_data_mat[1,:]
data_mat      = raw_data_mat[2:end,:]
data_mat      = convert(Array{Float64},data_mat)
################################################################################
##------------------ RES Data from Additional Felxibility File -----------------
################################################################################
raw_data_res_addt  = ods_readall(filename_addt;sheetsNames=[sheetname_res_addt],innerType="Matrix")
raw_data_res_addt  = raw_data_res_addt[sheetname_res_addt]   # Conversion from Dict to Array
header_res_addt    = raw_data_res_addt[1,:]
data_res_addt      = raw_data_res_addt[2:end,:]
data_res_addt      = convert(Array{Float64},data_res_addt)
################################################################################
##------- Curtailment Status Data from Additional Felxibility File -------------
################################################################################
raw_data_curt_addt = ods_readall(filename_addt;sheetsNames=[sheetname_curt_addt],innerType="Matrix")
raw_data_curt_addt = raw_data_curt_addt[sheetname_curt_addt]   # Conversion from Dict to Array
header_curt_addt   = raw_data_curt_addt[1,:]
data_curt_addt     = raw_data_curt_addt[2:end,1:size(header_curt_addt,1)]
data_curt_addt     = convert(Array{Float64},data_curt_addt)

header    = permutedims(hcat(permutedims(header_mat),permutedims(header_curt_addt[2:end])))
data      = vcat(data_mat,data_res_addt)

data      = hcat(data,data_curt_addt[:,2:end])
nGens     = Int64(size(data,1))
array_gens = data_reader(nGens,fields,header,data,Gens)

rheader_gens = header    # Exporting raw header of gens sheet
rdata_gens = data        # Exporting raw data of gens sheet
raw_data_mat       = nothing
raw_data_res_addt  = nothing
raw_data_curt_addt = nothing
header_mat         = nothing
header_res_addt    = nothing
header_curt_addt   = nothing
data_mat           = nothing
data_res_addt      = nothing
data_curt_addt     = nothing
header             = nothing
data               = nothing
#---------------------- Formatting the gen cost data ---------------------------
sheetname  = "Gens_cost_Addt";
fields  = ["bus_i","model","startUpCost","startDownCost","n","c2","c1","c0"]                                    # Fields that have to be read from the file
raw_data = ods_readall(filename_addt;sheetsNames=[sheetname],innerType="Matrix")
raw_data = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,[1;6:8]]
data      = raw_data[2:end,[1;6:8]]
data      = convert(Array{Float64}, data)

nGens         = Int64(size(data,1))
array_gcost   = data_reader(nGens,fields[[1;6:8],1],header,data,gen_cost)

rheader_gcost = header      # Exporting raw header of gens sheet
rdata_gcost   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

#------------------- Formatting the Active Profiles data -----------------------
nTP = 24
sheetname  = "pLoad_Profiles_Addt";
fields     = ["bus_i","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18","t19","t20","t21","t22","t23","t24"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename_addt;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array
header      = raw_data[1,:]
data        = raw_data[2:end,:]
# @infiltrate # This works. just use a fresh REPL
data        = convert(Array{Float64}, data)
nPprf_load  = Int64(size(data,1))

time_res_data = Int64(size(data[:,2:end],2))
data_new      = zeros(nPprf_load,nTP+1)
if time_res_data != nTP
    time_res_step = Int64(time_res_data/nTP)
    for i in 1:nTP
        data_new[:,i+1] = sum(data[:,((i-1)*time_res_step)+2:(i*time_res_step)+1],dims=2)/time_res_step
    end
    data_new[:,1] = data[:,1]
    data  = data_new
end
rheader_pProfile_load = header      # Exporting raw header of gens sheet
rdata_pProfile_load   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

#------------------- Formatting the Reactive Profiles data ---------------------
sheetname  = "qLoad_Profiles_Addt";
fields     = ["bus_i","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18","t19","t20","t21","t22","t23","t24"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename_addt;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
nQprf_load     = Int64(size(data,1))

time_res_data = Int64(size(data[:,2:end],2))
data_new      = zeros(nQprf_load,nTP+1)
if time_res_data != nTP
    time_res_step = Int64(time_res_data/nTP)
    for i in 1:nTP
        data_new[:,i+1] = sum(data[:,((i-1)*time_res_step)+2:(i*time_res_step)+1],dims=2)/time_res_step
    end
    data_new[:,1] = data[:,1]
    data  = data_new
end

rheader_qProfile_load = header      # Exporting raw header of gens sheet
rdata_qProfile_load   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing
#------------- Formatting the Active Generation Profiles data ------------------
sheetname  = "pGen_Profiles_Min_Addt";
fields     = ["bus_i","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18","t19","t20","t21","t22","t23","t24"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename_addt;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header     = raw_data[1,:]
data       = raw_data[2:end,:]
data       = convert(Array{Float64}, data)
nPprf_gen_min = Int64(size(data,1))

rheader_pProfile_gen_min = header      # Exporting raw header of gens sheet
rdata_pProfile_gen_min   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing
#------------- Formatting the Active Generation Profiles data ------------------
sheetname  = "pGen_Profiles_Max_Addt";
fields     = ["bus_i","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18","t19","t20","t21","t22","t23","t24"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename_addt;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header     = raw_data[1,:]
data       = raw_data[2:end,:]
data       = convert(Array{Float64}, data)
nPprf_gen_max = Int64(size(data,1))
rheader_pProfile_gen_max = header      # Exporting raw header of gens sheet
rdata_pProfile_gen_max   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

#------------- Formatting the Active Generation Profiles data ------------------
sheetname  = "qGen_Profiles_Min_Addt";
fields     = ["bus_i","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18","t19","t20","t21","t22","t23","t24"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename_addt;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header     = raw_data[1,:]
data       = raw_data[2:end,:]
data       = convert(Array{Float64}, data)
nQprf_gen_min = Int64(size(data,1))
rheader_qProfile_gen_min = header      # Exporting raw header of gens sheet
rdata_qProfile_gen_min   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing
#------------- Formatting the Active Generation Profiles data ------------------
sheetname  = "qGen_Profiles_Max_Addt";
fields     = ["bus_i","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18","t19","t20","t21","t22","t23","t24"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename_addt;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header     = raw_data[1,:]
data       = raw_data[2:end,:]
data       = convert(Array{Float64}, data)
nQprf_gen_max = Int64(size(data,1))
rheader_qProfile_gen_max = header      # Exporting raw header of gens sheet
rdata_qProfile_gen_max   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

#------------------- Formatting the Reactive Profiles data ---------------------
sheetname  = "Base_MVA";
fields     = ["base_MVA"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename_mat;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
nSbase    = Int64(size(data,1))
array_sbase  = data_reader(nSbase,fields,header,data,base_mva)
rheader_sbase = header      # Exporting raw header of gens sheet
rdata_sbase   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing
#----------------------- Formatting the Storage data ---------------------------
sheetname  = "Storage_Addt"
fields     = ["bus_i","Ps","Qs","energy","eRat","chRat","disRat","chEff","disEff","thermalRat","qmin","qmax","r","x","ploss","qloss","status","soc_0","socMin","socMax"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename_addt;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
nStr     = Int64(size(data,1))
array_storage   = data_reader(nStr,fields,header,data,energy_storage)
rheader_storage = header      # Exporting raw header of gens sheet
# i_soc_0 = findfirst(i-> isequal(i, "soc_0"), rheader_storage) # by baraa. debugging
rdata_storage   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

#---------------------- Formatting the gen cost data ---------------------------
sheetname  = "Storage_cost_Addt";
fields  = ["bus_i","c2","c1","c0"];                                    # Fields that have to be read from the file
raw_data = ods_readall(filename_addt;sheetsNames=[sheetname],innerType="Matrix")
raw_data = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
nStr      = Int64(size(data,1))
array_Strcost   = data_reader(nStr,fields,header,data,str_cost)
rheader_Strcost = header      # Exporting raw header of gens sheet
rdata_Strcost   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

###------------------------------ Slack node Profiles ------------------------##
sheetname_addt = "Profile_Slack"
raw_data_addt = ods_readall(filename_addt;sheetsNames=[sheetname_addt],innerType="Matrix")
raw_data_addt = raw_data_addt[sheetname_addt]   # Conversion from Dict to Array
header_addt   = raw_data_addt[1,:]
data_addt     = raw_data_addt[2:end,1:size(header_addt,1)]
data_addt     = convert(Array{Float64},data_addt)

nVprf  = Int64(size(data_addt,1))

time_res_data = Int64(size(data_addt[:,2:end],2))
data_new      = zeros(nVprf,nTP+1)
if time_res_data != nTP
    time_res_step = Int64(time_res_data/nTP)
    for i in 1:nTP
        data_new[:,i+1] = sum(data_addt[:,((i-1)*time_res_step)+2:(i*time_res_step)+1],dims=2)/time_res_step
    end
    data_new[:,1] = data_addt[:,1]
    data_addt  = data_new
end
rheader_vProfile = header_addt      # Exporting raw header of gens sheet
rdata_vProfile   = data_addt   # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

###--------------------------- OLTC Ratio Profiles --------------------------###
sheetname_addt = "Profile_oltc_ratio"
raw_data_addt = ods_readall(filename_addt;sheetsNames=[sheetname_addt],innerType="Matrix")
raw_data_addt = raw_data_addt[sheetname_addt]   # Conversion from Dict to Array
header_addt   = raw_data_addt[1,:]
data_addt     = raw_data_addt[2:end,1:size(header_addt,1)]
data_addt     = convert(Array{Float64},data_addt)

nOltc_ratio  = Int64(size(data_addt,1))

time_res_data = Int64(size(data_addt[:,2:end],2))
data_new      = zeros(nOltc_ratio,nTP+1)
if time_res_data != nTP
    time_res_step = Int64(time_res_data/nTP)
    for i in 1:nTP
        data_new[:,i+1] = sum(data_addt[:,((i-1)*time_res_step)+2:(i*time_res_step)+1],dims=2)/time_res_step
    end
    data_new[:,1] = data_addt[:,1]
    data_addt  = data_new
end
rheader_oltcProfile = header_addt      # Exporting raw header of gens sheet
rdata_oltcProfile   = data_addt   # Exporting raw data of gens sheet
rdata_oltcProfile[:,2:end] = rdata_oltcProfile[:,2:end].^-1
raw_data = nothing
header   = nothing
data     = nothing

return begin array_bus,
    rheader_buses,
    rdata_buses,
    array_lines,
    rheader_lines,
    rdata_lines,
    nTrsf_s1,
    array_trsf,
    rheader_trsf,
    rdata_trsf,
    array_loads,
    rheader_loads,
    rdata_loads,
    array_gens,
    rheader_gens,
    rdata_gens,
    array_gcost,
    rheader_gcost,
    rdata_gcost,
    rheader_pProfile_load,
    rdata_pProfile_load,
    rheader_qProfile_load,
    rdata_qProfile_load,
    rheader_pProfile_gen_min,
    rdata_pProfile_gen_min,
    nPprf_gen_max,
    rheader_pProfile_gen_max,
    rdata_pProfile_gen_max,
    nQprf_gen_min,
    rheader_qProfile_gen_min,
    rdata_qProfile_gen_min,
    nQprf_gen_max,
    rheader_qProfile_gen_max,
    rdata_qProfile_gen_max,
    array_sbase,
    rheader_sbase,
    rdata_sbase,
    array_storage,
    rheader_storage,
    rdata_storage,
    array_Strcost,
    rheader_Strcost,
    rdata_Strcost,
    rheader_vProfile,
    rdata_vProfile,
    rheader_oltcProfile,
    rdata_oltcProfile,
    nGens, nLines, nBus, nTrsf, nTP, nLoads,
    nOltc_ratio, nVprf, nStr, nSbase
    end
end
#--------------------------------- End -----------------------------------------
