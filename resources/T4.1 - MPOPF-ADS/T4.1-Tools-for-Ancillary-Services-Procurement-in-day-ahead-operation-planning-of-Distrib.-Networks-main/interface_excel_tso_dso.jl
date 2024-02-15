# This code reads each individual sheet of an ODS file and then process the available
# data in order to save the inforamtion related to each bus, line, load, gen etc.

# For this program to work, the number as well type of enteries mentioned in the
# field variable and data_type struct variables MUST be the same.
#---------------------------Formatting the Buses Data---------------------------
sheetname  = "Buses";
fields  = ["Bus_N","Bus_Type","Area","V_Mag (pu)","V_Ang (deg)","V_Nom (kV)","V_max (pu)","V_min (pu)"]; # Fields that have to be read from the file
raw_data = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data = raw_data[sheetname]   # Conversion from Dict to Array
header    = raw_data[1,:]
data      = raw_data[2:end,1:size(header,1)]
data      = convert(Array{Float64}, data)
# data_cont = zeros(size(fields,1))               # data_cont = Data Container
nBus      = Int64(size(data,1))
# global array_bus = Array{Bus}(undef,nBus,1)

# array_bus = @trace data_reader(array_bus,nBus,fields,header,data,data_cont,Bus)
array_bus =  data_reader(nBus,fields,header,data,Bus)
rheader_buses = header    # Exporting raw header of buses sheet
rdata_buses = data        # Exporting raw data of buses sheet
raw_data = nothing
header = nothing
data = nothing

#----------------------------Formatting the Lines Data--------------------------
sheetname  = "Lines";
fields  = ["From","To","r (pu)","x (pu)","g_sh (pu)","b_sh (pu)","RATE_A","RATE_B","RATE_C","tap_ratio","angle","br_status","ang_min","ang_max","tap_ratio_min","tap_ratio_max"]; # Fields that have to be read from the file
raw_data = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
#------------ Code for handling the lines whose status is zero -----------------
#------- This code will eliminate those lines whose status is set to zero-------
idx_br_status = findall(x->x=="br_status",header)
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
    y_eq = sum(data[ilines,5]+data[ilines,6]im,dims=1)
    amp_eq_A = sum(data[ilines,7],dims=1)
    amp_eq_B = sum(data[ilines,8],dims=1)
    amp_eq_C = sum(data[ilines,9],dims=1)
    nline_data = [transpose(data[ilines[1,1],1:2]) real(z_eq) imag(z_eq) real(y_eq) imag(y_eq) amp_eq_A amp_eq_B amp_eq_C transpose(data[ilines[1,1],10:end])]
    push!(new_line_data,nline_data)
end
data = vcat(new_line_data...)
##----------------------------------------------------------------------------##
# data_cont = zeros(size(fields,1))               # data_cont = Data Container
nLines    = Int64(size(data,1))
# global array_lines = Array{Line}(undef,nLines,1)

# array_lines = data_reader(array_lines,nLines,fields,header,data,data_cont,Line)
array_lines = data_reader(nLines,fields,header,data,Line)

rheader_lines = header    # Exporting raw header of lines sheet
rdata_lines = data        # Exporting raw data of lines sheet

raw_data = nothing
header = nothing
data = nothing
#------------------------ Formatting the Loads Data-----------------------------
sheetname  = "Loads";
fields  = ["Bus_N","Pd (MW)","Qd (MW)","Gs (MW)","Bs (MVAr)","Status","cost_inc","cost_dec"]; # Fields that have to be read from the file
raw_data = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
# data_cont = zeros(size(fields,1))               # data_cont = Data Container
nLoads    = Int64(size(data,1))
# global array_loads = Array{Loads}(undef,nLoads,1)

# array_loads = data_reader(array_loads,nLoads,fields,header,data,data_cont,Loads)
array_loads = data_reader(nLoads,fields,header,data,Loads)

rheader_loads = header    # Exporting raw header of loads sheet
rdata_loads   = data        # Exporting raw data of loads sheet
raw_data = nothing
header   = nothing
data     = nothing
#------------------------ Formatting the Gens Data------------------------------
sheetname  = "Gens";
fields  = ["Bus","Pg (MW)","Qg (MVAr)","Q_max (MVAr)","Q_min (MVAr)","V_set (pu)",
           "mBase (MVA)","P_max (MW)","P_min (MW)","Pc1 (MW)","Pc2 (MW)",
           "Qc1_min (MVAr)","Qc1_max (MVAr)","Qc2_min (MVAr)","Qc2_max (MVAr)",
           "ramp_agc (MW/min)","ramp_10 (MW)","ramp_30 (MW)","ramp_q (MVAr/min)","apf","status_curt","DG type"]; # Fields that have to be read from the file
raw_data = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
# data_cont = zeros(size(fields,1))               # data_cont = Data Container
nGens    = Int64(size(data,1))
# global array_gens = Array{Gens}(undef,nGens,1)

# array_gens = data_reader(array_gens,nGens,fields,header,data,data_cont,Gens)
array_gens = data_reader(nGens,fields,header,data,Gens)

rheader_gens = header    # Exporting raw header of gens sheet
rdata_gens = data        # Exporting raw data of gens sheet
raw_data = nothing
header = nothing
data = nothing

#---------------------- Formatting the gen cost data ---------------------------
sheetname  = "Gen_cost";
fields  = ["Bus","cost_a","cost_b","cost_c"];                                    # Fields that have to be read from the file
raw_data = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
# data_cont = zeros(size(fields,1))               # data_cont = Data Container
nGens     = Int64(size(data,1))
# global array_gcost = Array{gen_cost}(undef,nGens,1)

# array_gcost   = data_reader(array_gcost,nGens,fields,header,data,data_cont,gen_cost)
array_gcost   = data_reader(nGens,fields,header,data,gen_cost)

rheader_gcost = header      # Exporting raw header of gens sheet
rdata_gcost   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

#------------------- Formatting the Active Profiles data -----------------------
sheetname  = "P_Profiles_Load";
fields     = ["Bus","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18","t19","t20","t21","t22","t23","t24"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
# data_cont = zeros(size(fields,1))               # data_cont = Data Container
nPprf_load     = Int64(size(data,1))
# global array_pProfiles = Array{profile_P}(undef,nPprf,1)
#
# array_pProfiles  = data_reader(array_pProfiles,nPprf,fields,header,data,data_cont,profile_P)
rheader_pProfile_load = header      # Exporting raw header of gens sheet
rdata_pProfile_load   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

#------------------- Formatting the Reactive Profiles data ---------------------
sheetname  = "Q_Profiles_Load";
fields     = ["Bus","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18","t19","t20","t21","t22","t23","t24"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
# data_cont = zeros(size(fields,1))               # data_cont = Data Container
nQprf_load     = Int64(size(data,1))
# global array_qProfiles = Array{profile_Q}(undef,nQprf,1)
#
# array_pProfiles  = data_reader(array_qProfiles,nQprf,fields,header,data,data_cont,profile_Q)
rheader_qProfile_load = header      # Exporting raw header of gens sheet
rdata_qProfile_load   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

#------------------- Formatting the Active Profiles data -----------------------
sheetname  = "P_Profiles_Load_1041";
fields     = ["Bus","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18","t19","t20","t21","t22","t23","t24"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
# data_cont = zeros(size(fields,1))               # data_cont = Data Container
nPprf_load     = Int64(size(data,1))
# global array_pProfiles = Array{profile_P}(undef,nPprf,1)
#
# array_pProfiles  = data_reader(array_pProfiles,nPprf,fields,header,data,data_cont,profile_P)
rheader_pProfile_load_1041 = header      # Exporting raw header of gens sheet
rdata_pProfile_load_1041   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

#------------------- Formatting the Reactive Profiles data ---------------------
sheetname  = "Q_Profiles_Load_1041";
fields     = ["Bus","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18","t19","t20","t21","t22","t23","t24"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
# data_cont = zeros(size(fields,1))               # data_cont = Data Container
nQprf_load     = Int64(size(data,1))
# global array_qProfiles = Array{profile_Q}(undef,nQprf,1)
#
# array_pProfiles  = data_reader(array_qProfiles,nQprf,fields,header,data,data_cont,profile_Q)
rheader_qProfile_load_1041 = header      # Exporting raw header of gens sheet
rdata_qProfile_load_1041   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

#------------- Formatting the Active Generation Profiles data ------------------
sheetname  = "P_Profiles_Gen_Min";
fields     = ["Bus","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18","t19","t20","t21","t22","t23","t24"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header     = raw_data[1,:]
data       = raw_data[2:end,:]
data       = convert(Array{Float64}, data)
# data_cont  = zeros(size(fields,1))               # data_cont = Data Container
nPprf_gen_min = Int64(size(data,1))
# global array_qProfiles = Array{profile_Q}(undef,nQprf,1)
#
# array_pProfiles  = data_reader(array_qProfiles,nQprf,fields,header,data,data_cont,profile_Q)
rheader_pProfile_gen_min = header      # Exporting raw header of gens sheet
rdata_pProfile_gen_min   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing
#------------- Formatting the Active Generation Profiles data ------------------
sheetname  = "P_Profiles_Gen_Max";
fields     = ["Bus","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18","t19","t20","t21","t22","t23","t24"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header     = raw_data[1,:]
data       = raw_data[2:end,:]
data       = convert(Array{Float64}, data)
# data_cont  = zeros(size(fields,1))               # data_cont = Data Container
nPprf_gen_max = Int64(size(data,1))
# global array_qProfiles = Array{profile_Q}(undef,nQprf,1)
#
# array_pProfiles  = data_reader(array_qProfiles,nQprf,fields,header,data,data_cont,profile_Q)
rheader_pProfile_gen_max = header      # Exporting raw header of gens sheet
rdata_pProfile_gen_max   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

#------------- Formatting the Active Generation Profiles data ------------------
sheetname  = "Q_Profiles_Gen_Min";
fields     = ["Bus","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18","t19","t20","t21","t22","t23","t24"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header     = raw_data[1,:]
data       = raw_data[2:end,:]
data       = convert(Array{Float64}, data)
# data_cont  = zeros(size(fields,1))               # data_cont = Data Container
nQprf_gen_min = Int64(size(data,1))
# global array_qProfiles = Array{profile_Q}(undef,nQprf,1)
#
# array_pProfiles  = data_reader(array_qProfiles,nQprf,fields,header,data,data_cont,profile_Q)
rheader_qProfile_gen_min = header      # Exporting raw header of gens sheet
rdata_qProfile_gen_min   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing
#------------- Formatting the Active Generation Profiles data ------------------
sheetname  = "Q_Profiles_Gen_Max";
fields     = ["Bus","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14","t15","t16","t17","t18","t19","t20","t21","t22","t23","t24"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header     = raw_data[1,:]
data       = raw_data[2:end,:]
data       = convert(Array{Float64}, data)
# data_cont  = zeros(size(fields,1))               # data_cont = Data Container
nQprf_gen_max = Int64(size(data,1))
# global array_qProfiles = Array{profile_Q}(undef,nQprf,1)
#
# array_pProfiles  = data_reader(array_qProfiles,nQprf,fields,header,data,data_cont,profile_Q)
rheader_qProfile_gen_max = header      # Exporting raw header of gens sheet
rdata_qProfile_gen_max   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

#------------------- Formatting the Reactive Profiles data ---------------------
sheetname  = "Base_MVA";
fields     = ["base_MVA"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
# data_cont = zeros(size(fields,1))               # data_cont = Data Container
nSbase    = Int64(size(data,1))
# global array_sbase = Array{base_mva}(undef,nSbase,1)

# array_sbase  = data_reader(array_sbase,nSbase,fields,header,data,data_cont,base_mva)
array_sbase  = data_reader(nSbase,fields,header,data,base_mva)
rheader_sbase = header      # Exporting raw header of gens sheet
rdata_sbase   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing
#----------------------- Formatting the Storage data ---------------------------
sheetname  = "Storage";
fields     = ["Bus","Ps","Qs","Energy (MWh)","E_rating (MWh)","Charge_rating (MW)","Discharge_rating (MW)","Charge_efficiency","Discharge_efficiency","Thermal_rating (MVA)","Qmin (MVAr)","Qmax (MVAr)","R","X","P_loss","Q_loss","Status","soc_initial","soc_min","soc_max"];                                    # Fields that have to be read from the file
raw_data   = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data   = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
# data_cont = zeros(size(fields,1))               # data_cont = Data Container
nStr     = Int64(size(data,1))
# global array_storage = Array{energy_storage}(undef,nStr,1)
#
# array_storage   = data_reader(array_storage,nStr,fields,header,data,data_cont,energy_storage)
array_storage   = data_reader(nStr,fields,header,data,energy_storage)
rheader_storage = header      # Exporting raw header of gens sheet
rdata_storage   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

#---------------------- Formatting the gen cost data ---------------------------
sheetname  = "Storage_cost";
fields  = ["Bus","cost_a","cost_b","cost_c"];                                    # Fields that have to be read from the file
raw_data = ods_readall(filename;sheetsNames=[sheetname],innerType="Matrix")
raw_data = raw_data[sheetname]   # Conversion from Dict to Array

header    = raw_data[1,:]
data      = raw_data[2:end,:]
data      = convert(Array{Float64}, data)
# data_cont = zeros(size(fields,1))               # data_cont = Data Container
nStr      = Int64(size(data,1))
# global array_Strcost = Array{str_cost}(undef,nStr,1)

# array_Strcost   = data_reader(array_Strcost,nStr,fields,header,data,data_cont,str_cost)
array_Strcost   = data_reader(nStr,fields,header,data,str_cost)
rheader_Strcost = header      # Exporting raw header of gens sheet
rdata_Strcost   = data        # Exporting raw data of gens sheet
raw_data = nothing
header   = nothing
data     = nothing

#--------------------------------- End -----------------------------------------
