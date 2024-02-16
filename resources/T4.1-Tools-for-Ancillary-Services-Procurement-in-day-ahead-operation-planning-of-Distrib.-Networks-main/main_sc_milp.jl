using Pkg
Pkg.activate("smpopf_usman")

################################################################################
# This following code corresponds to the Ancillary Services Procurement Tool in
# the Day-Ahead Operation Planning of Active Distribution System. Please cite the
# following work if you use this algorithm in your work. Thanks!

# Usman, M., & Capitanescu, F. (2022). A Novel Tractable Methodology to Stochastic Multi-Period AC OPF in Active Distribution Systems Using Sequential Linearization Algorithm. IEEE Transactions on Power Systems.

# @article{usman2022novel,
#       title={A Novel Tractable Methodology to Stochastic Multi-Period AC OPF in Active Distribution Systems Using Sequential Linearization Algorithm},
#       author={Usman, Muhammad and Capitanescu, Florin},
#       journal={IEEE Transactions on Power Systems},
#       year={2022},
#       publisher={IEEE}
# }
################################################################################
using CPLEX
using Cbc
using Dates
using Ipopt
using JuMP
using Juniper
using LinearAlgebra
using MathOptInterface
using OdsIO
using Plots
using XLSX
using DataFrames

using Statistics
using Printf
using Alert
using WAV # for notification sound
if false
    using Infiltrator
    Infiltrator.clear_disabled!()
    Infiltrator.toggle_async_check(false)
    # To use infiltrator, just place "@infiltrate" alone in a separate line, where you want to break
    # and run in a fresh REPL. do NOT re-run in an existing REPL.
    # https://docs.juliahub.com/Infiltrator/ge3PS/0.3.0/
    # Infiltrator shortcuts
    # https://github.com/JuliaDebug/Infiltrator.jl
    # ?: Print this help text.
    # @trace: Print the current stack trace.
    # @locals: Print local variables. @locals x y only prints x and y.
    # @exfiltrate: Save all local variables into the store. @exfiltrate x y saves x and y; this variant can also exfiltrate variables defined in the infil> REPL.
    # @toggle: Toggle infiltrating at this @infiltrate spot (clear all with Infiltrator.clear_disabled!()).
    # @continue: Continue to the next infiltration point or exit (shortcut: Ctrl-D).
    # @doc symbol: Get help for symbol (same as in the normal Julia REPL).
    # @exit: Stop infiltrating for the remainder of this session and exit (on Julia versions prior to 1.5 this needs to be manually cleared with Infiltrator.end_session!()).
end

# To remove the warnings from the log, use REGEXP in notepad++ to remove the following:
# ^┌.+\r\n(.+\r\n)?└ @.+$\r\n

#-------------------Accessing current folder directory--------------------------
println(pwd());
cd(dirname(@__FILE__))
println(pwd());
################################################################################
#### -----------  Compilation of User-Defined Functions -------------------#####
################################################################################
include("ac_power_flow_func_check.jl")
include("ac_power_flow_model_check.jl")

include("ac_power_flow_model_opf_int.jl")
include("ac_power_flow_opf_functions_int.jl")

include("ac_power_flow_model_opf.jl")
include("ac_power_flow_opf_functions.jl")

include("data_nw_functions.jl")                                                  # Transformer data file is independant of multi-period and scenarios
include("power_flow_funcitons_milp.jl")
include("loop_finding.jl")

include("opf_functions.jl")
include("sql_loop.jl")
include("sql_opf_model.jl")
include("opf_model_sub_functions.jl")
include("model_sol_recovery.jl")

include("data_nw_functions_pf.jl")                                               # Transformer data file is indepdendant of multi-period and scenarios
include("non_linear_power_equations.jl")

include("data_types.jl")                                                         # Reading the structure of network related fields
include("data_reader.jl")                                                        # Function setting the data corresponding to each network quantity
include("Append2Excel.jl")                                                       # Added by Baraa

################################################################################
######################  Select the Input file ##################################
################################################################################
CasesList =    Dict([
    ("UK_DX_01_BaseCase"                                         , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/BaseCase.ods"                                                 ),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_S_"            , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2030_Distribution_Network_Urban_UK_S.ods"            ),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_S_Without_1_7_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2030_Distribution_Network_Urban_UK_S_Without_1_7.ods"),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_S_Without_2_7_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2030_Distribution_Network_Urban_UK_S_Without_2_7.ods"),

    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_W_"            , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2030_Distribution_Network_Urban_UK_W.ods"            ),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_W_Without_1_7_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2030_Distribution_Network_Urban_UK_W_Without_1_7.ods"),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_W_Without_2_7_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2030_Distribution_Network_Urban_UK_W_Without_2_7.ods"),

    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_S_"            , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2040_Distribution_Network_Urban_UK_S.ods"            ),
    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_S_Without_1_7_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2040_Distribution_Network_Urban_UK_S_Without_1_7.ods"),
    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_S_Without_2_7_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2040_Distribution_Network_Urban_UK_S_Without_2_7.ods"),

    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_W_"            , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2040_Distribution_Network_Urban_UK_W.ods"            ),
    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_W_Without_1_7_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2040_Distribution_Network_Urban_UK_W_Without_1_7.ods"),
    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_W_Without_2_7_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2040_Distribution_Network_Urban_UK_W_Without_2_7.ods"),

    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_S_"            , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2050_Distribution_Network_Urban_UK_S.ods"            ),
    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_S_Without_1_7_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2050_Distribution_Network_Urban_UK_S_Without_1_7.ods"),
    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_S_Without_2_7_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2050_Distribution_Network_Urban_UK_S_Without_2_7.ods"),

    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_W_"            , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2050_Distribution_Network_Urban_UK_W.ods"            ),
    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_W_Without_1_7_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2050_Distribution_Network_Urban_UK_W_Without_1_7.ods"), # CPLEX Error  3019: Failure to solve MIP subproblem
    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_W_Without_2_7_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/Obsolete_but_uploaded/UK_Dx_01_2050_Distribution_Network_Urban_UK_W_Without_2_7.ods"),

    ("HR_DX_01/BROWN_BaseCase"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/BROWN/BaseCase.ods"                       ),
    ("HR_Dx_01_2030_Su_Bd_brown_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/BROWN/HR_Dx_01_2030_Su_Bd_brown_WithFlex_BROWN.ods"),
    ("HR_Dx_01_2030_W_Bd_brown_" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/BROWN/HR_Dx_01_2030_W_Bd_brown_WithFlex_BROWN.ods" ), # sql finishes successfully. only FL actuated
    ("HR_Dx_01_2040_Su_Bd_brown_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/BROWN/HR_Dx_01_2040_Su_Bd_brown_WithFlex_BROWN.ods"), # sql_loop fails on vm_ac not recognized
    ("HR_Dx_01_2040_W_Bd_brown_" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/BROWN/HR_Dx_01_2040_W_Bd_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2050_Su_Bd_brown_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/BROWN/HR_Dx_01_2050_Su_Bd_brown_WithFlex_BROWN.ods"),
    ("HR_Dx_01_2050_W_Bd_brown_" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/BROWN/HR_Dx_01_2050_W_Bd_brown_WithFlex_BROWN.ods" ),

    ("HR_DX_01/RED_BaseCase"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/RED/BaseCase.ods"                         ),
    ("HR_Dx_01_2030_Su_Bd_red_"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/RED/HR_Dx_01_2030_Su_Bd_red_WithFlex_RED.ods"      ),
    ("HR_Dx_01_2030_W_Bd_red_"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/RED/HR_Dx_01_2030_W_Bd_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2040_Su_Bd_red_"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/RED/HR_Dx_01_2040_Su_Bd_red_WithFlex_RED.ods"      ),
    ("HR_Dx_01_2040_W_Bd_red_"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/RED/HR_Dx_01_2040_W_Bd_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2050_Su_Bd_red_"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/RED/HR_Dx_01_2050_Su_Bd_red_WithFlex_RED.ods"      ),
    ("HR_Dx_01_2050_W_Bd_red_"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/RED/HR_Dx_01_2050_W_Bd_red_WithFlex_RED.ods"       ),  # UndefVarError: vm_ac not defined

    ("HR_DX_01/GREEN_BaseCase"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/GREEN/BaseCase.ods"                       ),
    ("HR_Dx_01_2030_Su_Bd_green_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/GREEN/HR_Dx_01_2030_Su_Bd_green_WithFlex_GREEN.ods"),
    ("HR_Dx_01_2030_W_Bd_green_" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/GREEN/HR_Dx_01_2030_W_Bd_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2040_Su_Bd_green_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/GREEN/HR_Dx_01_2040_Su_Bd_green_WithFlex_GREEN.ods"),
    ("HR_Dx_01_2040_W_Bd_green_" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/GREEN/HR_Dx_01_2040_W_Bd_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2050_Su_Bd_green_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/GREEN/HR_Dx_01_2050_Su_Bd_green_WithFlex_GREEN.ods"),
    ("HR_Dx_01_2050_W_Bd_green_" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/Obsolete but uploaded/GREEN/HR_Dx_01_2050_W_Bd_green_WithFlex_GREEN.ods" ),

    ("ES_DX_03_BaseCase"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/ES DX 03/Newer and uploaded/BaseCase.ods"          ),
    ("ES_Dx_03_2030_S_Bd_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/ES DX 03/Newer and uploaded/ES_Dx_03_2030_S_Bd.ods"),
    ("ES_Dx_03_2030_W_Bd_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/ES DX 03/Newer and uploaded/ES_Dx_03_2030_W_Bd.ods"),
    ("ES_Dx_03_2040_S_Bd_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/ES DX 03/Newer and uploaded/ES_Dx_03_2040_S_Bd.ods"),
    ("ES_Dx_03_2040_W_Bd_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/ES DX 03/Newer and uploaded/ES_Dx_03_2040_W_Bd.ods"),
    ("ES_Dx_03_2050_S_Bd_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/ES DX 03/Newer and uploaded/ES_Dx_03_2050_S_Bd.ods"),
    ("ES_Dx_03_2050_W_Bd_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/ES DX 03/Newer and uploaded/ES_Dx_03_2050_W_Bd.ods"),

    ("PT_DX_01_BaseCase"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/OlderButUploadedToCloud/ReinforcedLines/BaseCase.ods"           ), # UndefVarError: vm_ac not defined
    ("PT_Dx_01_2030_Su_Bd_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/OlderButUploadedToCloud/ReinforcedLines/PT_Dx_01_2030_Su_Bd_WithFlex.ods"),
    ("PT_Dx_01_2030_W_Bd_" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/OlderButUploadedToCloud/ReinforcedLines/PT_Dx_01_2030_W_Bd_WithFlex.ods" ),
    ("PT_Dx_01_2040_Su_Bd_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/OlderButUploadedToCloud/ReinforcedLines/PT_Dx_01_2040_Su_Bd_WithFlex.ods"),
    ("PT_Dx_01_2040_W_Bd_" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/OlderButUploadedToCloud/ReinforcedLines/PT_Dx_01_2040_W_Bd_WithFlex.ods" ),
    ("PT_Dx_01_2050_Su_Bd_", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/OlderButUploadedToCloud/ReinforcedLines/PT_Dx_01_2050_Su_Bd_WithFlex.ods"),
    ("PT_Dx_01_2050_W_Bd_" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/OlderButUploadedToCloud/ReinforcedLines/PT_Dx_01_2050_W_Bd_WithFlex.ods" ),

    ("HR_Demo_TopoIsland1" , "input_data/DemonstrationData_July2023/kc35_ODS/DEMO2_kc35_TopoIsland_1.ods" ),
    ("HR_Demo_TopoIsland2" , "input_data/DemonstrationData_July2023/kc35_ODS/DEMO2_kc35_TopoIsland_2.ods" ),
    ("HR_Demo_TopoIsland3" , "input_data/DemonstrationData_July2023/kc35_ODS/DEMO2_kc35_TopoIsland_3.ods" ),
    # ("HR_Demo_None"        , "input_data/DemonstrationData_July2023/kc35_ODS/DEMO2_kc35.ods" ),

    #=
    ### multiline break
    remaining issues mentioned on August 17th:
    fix the dictionary of cases (CasesList) such that the files in ATTEST cloud drive would work right away
    Some cases fail with an error message: can't find tap settings or something like that. investigate this

    mount the files they sent and make sure that a simulation can start
    =#
    ])

IncludeStats    = false # Based upon the request of ATTEST partners, to remove last 4 lines of stats (MAX, MIN, AVG, TOT)
CaseNames       = collect(keys(CasesList))
CaseNames       = sort(CaseNames)[end:-1:1]
N_Cases         = length(CasesList) # sizeof(CasesList) here returns something else that i dont understand
# ==============================================================================
Flex_Case = ""
=
#ACTIVATE THIS WITH UK CASES
if true # Disabled temporarily. should be enabled.
    Flex_Case = "WithFlex" # Edit the case name. comment out if not applicable
else
    Flex_Case = "WOFlex" # Edit the case name. comment out if not applicable
end
Flex_Case = string("_",Flex_Case)
# =#

CustomCaseName = "" # Edit the case name. comment out if not applicable
# ==============================================================================
if  false # specify which case to run, by specifying the case's name
    TheKey      = "HR_Dx_01_2040_Su_Bd_green_" # PT_Dx_01_2040_W_Bd_Reinforced # PT_Dx_01_2040_Su_Bd_Reinforced
    i           = findfirst(isequal(TheKey),CaseNames)
    #### i = findall(x -> occursin(TheKey, x),CaseNames) # also works

else    # specify which case to run, by specifying the case's number
    # i         = 58 # there are 57 cases
    i = 1
    TheKey    = CaseNames[i]
end

if sizeof(CustomCaseName)>0
    CustomCaseName = string("_",CustomCaseName)
end
# ==============================================================================
OutLog      = "output_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/OutLog.xlsx"
LogThis(;CaseLog, Append=false) = append_xl_row(workbook_path=OutLog,
                                                sheet_name="Log",
                                                row_data=collect(CaseLog),
                                                Append=Append)

LogLabels   = ["#","CaseName","InputFile","OutputFile","Start_Time", "Init_PF","N_Violations","Opt_Started","Opt_Found"]
if isfile(OutLog)
    XLS_MODE = "rw"
else
    XLS_MODE = "w" # file does NOT exist. use write mode.
end
XLSX.openxlsx(OutLog, mode=XLS_MODE) do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "Log")
    sheet["A1"] = LogLabels
end
# ==============================================================================

filename_mat    = CasesList[TheKey]
filename_mat    = replace(filename_mat,'/'=>'\\')
filename_mat    = replace(filename_mat  , ".ods"      =>string(Flex_Case,".ods")) # added temporarily. should be disabled.
filename_addt   = replace(filename_mat  , ".ods"      =>"_flex.ods")
file_op         = replace(filename_mat  , "input_data"=>"output_data")
file_op         = replace(file_op       , ".ods"      =>string(Flex_Case,CustomCaseName,"_output.xlsx"))

@printf("Case #%3i/%i: %s\n", i, N_Cases, TheKey)
StartTime = now()

CaseLog = [i, string(TheKey,Flex_Case,CustomCaseName), filename_addt, file_op, string(StartTime), "Did not start", false, false, false]
# println("Writing to Excel: Append = False;\tData = ", string(hcat(CaseLog[[1,6:end]])))

# open(replace(file_op, ".ods" => ".txt"), "w") do out
#     redirect_stdout(out) do

println("Running the case: ",filename_mat, " @", StartTime)
# println("This is CaseLog: ", join(CaseLog,", "))
LogThis(CaseLog=CaseLog, Append=true)

################################################################################
### -------------------- Reading Input Data Files -------------------------- ###
################################################################################
include("interface_excel.jl")  # Reading Input Network data and additional flexibility data corresponding to the network
if isfile(filename_mat)==false
throw(error(string("Could not find network data: \"", filename_mat, "\" not found")))
elseif isfile(filename_addt)==false
    throw(error(string("Could not find flex data \"", filename_addt, "\" not found")))
end

(array_bus,
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
nOltc_ratio, nVprf, nStr, nSbase) = interface_excel_fun(filename_mat,filename_addt)

filename = "input_data/scenario_gen.ods"                                         # Reading Scenario files generated by Scenario Generation Tool Box
include("interface_scenarios.jl") # this file builds the time profile for loads.

(nPprf_wind_sc,
nSc,
rheader_pProfile_wind_sc,
rdata_pProfile_wind_sc,
nPprf_pv_sc,
rheader_pProfile_pv_sc,
rdata_pProfile_pv_sc,
n_prob,
rheader_prob_sc,
rdata_prob_sc) = interface_scenarios_fun(filename)

################################################################################
### ---------------------- Constants Declaration -------------------------- ####
################################################################################
# "INCLUDE" actually runs/executes the file. not just makes it available to be called.
include("constants.jl") # this file contains all the parameters selected by the programmer/user
include("constants_pf_opf.jl")

if  occursin("WOFlex", CaseLog[2])
    println("Case withOUT Flexibility. Skipping")
    # continue
    flex_fl   = 0
    fl_bin    = 0
    flex_str  = 0
    str_bin   = 0

    flex_oltc = 0    # On-load tap changing transfomer
    oltc_bin  = 0    # Binary variables for the modeling of OLTC transformer

    # There are others which i'm not changing. such as APC.
    println("case has no flexibility. Storage and Flexible_Load variables were set to zero")
else
    flex_fl   = 1
    fl_bin    = 1
    flex_str  = 1
    str_bin   = 1
end

## ------------------------------- Network Data ---------------------------#####
(idx_from_line,
 idx_to_line,
 yij_line,
 yij_line_sh,
 tap_ratio,
 tap_ratio_min,
 tap_ratio_max,
 nTrsf,
 dLines) = data_lines(array_lines,
                      nw_lines,
                      nLines)

idx_tap = data_trsf(dLines)
dLines = [dLines idx_tap]

node_data = data_nodes(nBus,
                       nw_buses,
                       idx_from_line,
                       idx_to_line,
                       yij_line,
                       yij_line_sh,
                       tap_ratio,
                       tap_ratio_min,
                       tap_ratio_max,
                       idx_tap,
                       node)

(tdata, itap_tr, itap_0) = data_trsf_tap(node_data)
itap_r = data_trsf_tap_index(tdata, itap_tr, itap_0)

float_itap = union(itap_tr, setdiff(itap_0, itap_r))
itap = [Int(i) for i in float_itap]
tdata = tdata[itap, :]

(trsf_fict_nodes,
 tap_ratio_range,
 tratio_init,
 n_tap) = fictitious_trsf_node_new(tap_ratio,
                                    rdata_buses,
                                    tdata,
                                    nw_trsf,
                                    rdata_trsf)

(cfl_inc,
 cfl_dec,
 cost_load_inc,
 cost_load_dec,
 bus_data_lsheet,
 idx_Gs_lsheet,
 idx_Bs_lsheet,
 idx_St_lsheet,
 nFl,
 iFl,
 nd_fl) = data_load(nw_loads,
                    rheader_loads,
                    rdata_loads,
                    nLoads,
                    nTP,
                    nSc,
                    sbase)

(bus_data_gsheet,
 i_ncurt_gens,
 i_curt_gens,
 nNcurt_gen,
 nCurt_gen,
 nd_ncurt_gen,
 nd_curt_gen,
 Pg_max,
 Pg_min,
 Qg_max,
 Qg_min,
 cA_gen,
 cB_gen,
 cC_gen) = data_gen(nw_gens,
                    nw_gcost,
                    rheader_gens,
                    rdata_gens,
                    nGens,
                    nTP)

(bus_data_Ssheet,
 bus_data_Strcost,
 idx_St_Strsheet,
 iStr_active,
 nStr_active,
 nd_Str_active,
 cA_str,
 cB_str,
 cC_str,
 cost_a_str,
 cost_b_str,
 cost_c_str,
 soc_0) = data_storage(rheader_storage,
                            rdata_storage,
                            nw_storage,
                            nw_Strcost,
                            nTP,
                            nSc,
                            sbase)

(p_load,
 q_load,
 pg_max,
 pg_min,
 qg_max,
 qg_min,
 cost_a_gen,
 cost_b_gen,
 cost_c_gen) = data_load_gen(nLoads,
                            nGens,
                            nTP,
                            nSc,
                            bus_data_lsheet,
                            bus_data_gsheet,
                            Pg_max,
                            Pg_min,
                            Qg_max,
                            Qg_min,
                            nw_pPrf_data_load,
                            nw_qPrf_data_load,
                            nw_pPrf_data_gen_max,
                            nw_pPrf_data_gen_min,
                            nw_qPrf_data_gen_max,
                            nw_qPrf_data_gen_min,
                            nw_loads,
                            sbase,
                            scenario,
                            scenario_data_p_max,
                            scenario_data_p_min,
                            scenario_data_q_max,
                            scenario_data_q_min,
                            cA_gen,
                            cB_gen,
                            cC_gen)

graph = rdata_lines[:, 1:2]
graph = convert.(Int64, graph)
cycleList = cycles_finding(graph, numCycles, cycles)                               # Loop finding Code
oltc_ratio = zeros(Float64, (nSc, nTP, nTrsf))

## -------------------------- AC Power Flow ---------------------------------###
# vol_nodes_mag_pf   = zeros(nSc_pf,nTP_pf,size(nw_buses,1))
# vol_nodes_theta_pf = zeros(nSc_pf,nTP_pf,size(nw_buses,1))
##-----------------Scheme to run different types of AC-Power Flow ----------####
# 1. Deterministic SP: scenario = 0, nSc_pf = 1, nTP_pf = 1
# 2. Deterministic MP: scenario = 0, nSc_pf = 1, nTP_pf = 24
# 3. Stochastic MP (Avg): scenario = 1, nSc_pf = 1, nTP_pf = 24
# 4. Stochastic MP (Full): scenario = 1, nSc_pf = 10, nTP_pf = 24
# 5. Stochastic SP is not programmed
##------------------------------------------------------------------------------
idx_curt_dg = zeros(size(i_curt_gens, 1))
for i in 1:size(i_curt_gens, 1)
    idx_curt_dg[i, 1] = i_curt_gens[i][1]
end
idx_curt_dg = convert.(Int64, idx_curt_dg)

p_curt_pf   = zeros(nSc_pf, nTP_pf, nCurt_gen) ./ sbase
p_dis_pf    = zeros(nSc_pf, nTP_pf, nStr_active) ./ sbase
p_ch_pf     = zeros(nSc_pf, nTP_pf, nStr_active) ./ sbase
p_od_pf     = zeros(nSc_pf, nTP_pf, nFl) ./ sbase
p_ud_pf     = zeros(nSc_pf, nTP_pf, nFl) ./ sbase
q_dg_pf     = tan(acos(dg_ac_pf)) .* (pg_max[:, :, idx_curt_dg] - p_curt_pf)
oltc_tap_init = tratio_init                         # Unless the profile of tap ratio is given, it is fixed to the same value in all time-periods and all scenarios
p_strg_pf = p_dis_pf - p_ch_pf

pg_max_pf   = pg_max[pf_scenario, :, :]
pg_min_pf   = pg_min[pf_scenario, :, :]
qg_max_pf   = qg_max[pf_scenario, :, :]
qg_min_pf   = pg_min[pf_scenario, :, :]

pg_max_pf   = reshape(pg_max_pf, (nSc_pf_new, nTP_pf_new, nGens))
pg_min_pf   = reshape(pg_min_pf, (nSc_pf_new, nTP_pf_new, nGens))
qg_max_pf   = reshape(qg_max_pf, (nSc_pf_new, nTP_pf_new, nGens))
qg_min_pf   = reshape(qg_min_pf, (nSc_pf_new, nTP_pf_new, nGens))
prob_scs_pf = prob_scs[pf_scenario, 1]

q_dg_pf     = tan(acos(dg_ac_pf)) .* (pg_max[:, :, idx_curt_dg] - p_curt_pf)
oltc_tap_init = tratio_init                         # Unless the profile of tap ratio is given, it is fixed to the same value in all time-periods and all scenarios

idx_status_oltc = findall(x -> x == 1, rdata_trsf[:, end])
idx_slack_bus_pf = findall(x -> x == slack_bus_type, rdata_buses[:, 2])
if !isempty(idx_status_oltc)
    for s in 1:nSc_pf
        for i in 1:size(idx_status_oltc, 1)
            oltc_tap_init[s, :, idx_status_oltc[i, 1]] = rdata_oltcProfile[idx_status_oltc[i, 1], 2:nTP+1]
        end
    end
end
vm_ac = zeros(nBus, nTP, nSc)
va_ac = zeros(nBus, nTP, nSc)

##################################################
display("Initial AC Optimum Power Flow before MILP Model")
CaseLog[6] = "Started only" # 6:"Init_PF"
LogThis(CaseLog=CaseLog, Append=false)
(vol_nodes_mag_pf,
vol_nodes_theta_pf,
vol_rect_pf,
ObjValue) = ac_power_flow_opf_int(nBus,
                                    nLines,
                                    nNcurt_gen,
                                    rdata_buses,
                                    nSc_pf,
                                    nTP_pf,
                                    prob_scs_pf,
                                    time_step,
                                    stoch_model,
                                    p_curt_pf,
                                    p_od_pf,
                                    p_ud_pf,
                                    p_strg_pf,
                                    q_dg_pf,
                                    nw_buses,
                                    nw_lines,
                                    rdata_loads,
                                    bus_data_lsheet,
                                    bus_data_Ssheet,
                                    node_data,
                                    nd_fl,
                                    nd_curt_gen,
                                    nd_ncurt_gen,
                                    p_load,
                                    q_load,
                                    idx_Gs_lsheet,
                                    idx_Bs_lsheet,
                                    yii_sh,
                                    i_curt_gens,
                                    sbase,
                                    flex_oltc,
                                    vbase,
                                    Ibase,
                                    pg_min_pf,
                                    pg_max_pf,
                                    qg_min_pf,
                                    qg_max_pf,
                                    yij_line,
                                    dLines,
                                    i_ncurt_gens,
                                    error_msg,
                                    pgen_tol,
                                    qgen_tol,
                                    oltc_tap_init,
                                    vol_cstr_tol)
CaseLog[6] = true # 6:"Init_PF"
LogThis(CaseLog=CaseLog, Append=false)

##################################################

(br_crnt_pu_pf_int,
 br_crnt_si_pf_int,
 br_pwr_pu_pf_int,
 br_pwr_si_pf_int) = recovery_branch_current_pf(nSc,
                                                nTP,
                                                nLines,
                                                nw_lines_pf,
                                                vol_rect_pf,
                                                yij_line,
                                                Ibase,
                                                sbase,
                                                ybase,
                                                vbase,
                                                oltc_tap_init,
                                                dLines)

(vol_viol_int,
 crnt_viol_int,
 max_vol_viol_int,
 max_crnt_viol_int,
 viol_nodes_int,
 nVol_viol_int,
 avg_vol_viol_int,
 above_avg_vol_viol_int,
 below_avg_vol_viol_int,
 nCrnt_viol_int,
 avg_crnt_viol_int,
 above_avg_crnt_viol_int,
 below_avg_crnt_viol_int,
 vol_viol_cn_lin_int) = constraint_violation_pf(nBus,
                                                nw_buses,
                                                nLines,
                                                nw_lines,
                                                nTP,
                                                nSc,
                                                vol_nodes_mag_pf,
                                                br_crnt_si_pf_int,
                                                I_rat,
                                                1)

N_violations = [nVol_viol_int,
                nCrnt_viol_int]

##################################################
################ Added by baraa ##################
function Append2Excel(xf, Data,SheetName)
        if !XLSX.hassheet(xf, SheetName)
            XLSX.addsheet!(xf, SheetName)
        end
        if size(Data,2)>16384
            error("Object has $(size(Data,2)) columns. Maximum allowed number of columns is 16384")
        end

        if size(Data,1)>0
            XLSX.writetable!(xf[SheetName],
                            Data,
                            ;anchor_cell=XLSX.CellRef("A1"))
        else
            println("Array is empty. Skipping")
        end
end

# --------------------------------------------
function Flatten3DArray(TheMatrix::Array{Float64, 3})
    PermutedMatrix = PermutedDimsArray(TheMatrix,(3,2,1))
    EmptyArray  = Array{Float64}(undef, 0, size(TheMatrix,2)+2)
    for i in 1:size(TheMatrix,1)
        Temp_Mat        = hcat(PermutedMatrix[:,1:2,i],PermutedMatrix[:,:,i])
        Temp_Mat[:,1]  .= i
        Temp_Mat[:,2]  .= 1:size(Temp_Mat,1)
        EmptyArray = vcat(EmptyArray, Temp_Mat)
    end
    # println(size(EmptyArray))
    return EmptyArray
end

# --------------------------------------------
PF_headers = vcat("Scen","ID",[string("T",i) for i in 1:24])

XLSX.writetable(file_op,
                VOLT      = (collect(eachcol(Flatten3DArray(vol_nodes_mag_pf))), PF_headers),
                Crnt_PU   = (collect(eachcol(Flatten3DArray(br_crnt_pu_pf_int))), PF_headers),
                Crnt_SI   = (collect(eachcol(Flatten3DArray(br_crnt_si_pf_int))), PF_headers),
                P_load    = (collect(eachcol(Flatten3DArray(p_load))), PF_headers),
                Q_load    = (collect(eachcol(Flatten3DArray(q_load))), PF_headers),
                Pg_max    = (collect(eachcol(Flatten3DArray(pg_max))), PF_headers),
                Qg_max    = (collect(eachcol(Flatten3DArray(qg_max))), PF_headers),
                ; overwrite=true)

XLSX.openxlsx(file_op, mode="rw") do xf
    if  sizeof(vol_viol_int) > 0
        WriteThis  = [map(string,vol_viol_int[i]) for i in 1:size(vol_viol_int,1)]
        WriteThis2 = hcat(collect.(WriteThis)...)
        # WriteThis3 = DataFrame(WriteThis2, :auto)
        WriteThis3 = DataFrame( collect.(eachrow(WriteThis2)), ["Itr","Num","Scen","T","Bus","V","Viol(V)"])
        Append2Excel(xf, WriteThis3, "Vlt_Viol")
    end

    if  sizeof(crnt_viol_int)>0
        WriteThis  = [map(string,crnt_viol_int[i]) for i in 1:size(crnt_viol_int,1)]
        WriteThis2 = hcat(collect.(WriteThis)...)
        # WriteThis3 = DataFrame(WriteThis2, :auto)
        WriteThis3 = DataFrame( collect.(eachrow(WriteThis2)), ["Itr","Num","Scen","T","F_Bus","T_Bus","S","S_limit","Viol(S)"])
        Append2Excel(xf, WriteThis3, "Crnt_Viol")
    end
end

println("Checkout initial OPF results @ ",file_op)
CaseLog[6] = true                   # 6:"Init_PF", 7:"Viols_Exist", 8:"Opt_Commenced", 9:"Opt_Found"
CaseLog[7] = sum(N_violations)      # 6:"Init_PF", 7:"Viols_Exist", 8:"Opt_Commenced", 9:"Opt_Found"

if sum(N_violations)==0
    alert("JULIA: NO violations. This case is a waste of time")
    CaseLog[8] = "Unnecessary"      # 6:"Init_PF", 7:"Viols_Exist", 8:"Opt_Commenced", 9:"Opt_Found"
    CaseLog[9] = "Unnecessary"      # 6:"Init_PF", 7:"Viols_Exist", 8:"Opt_Commenced", 9:"Opt_Found"
    LogThis(CaseLog=CaseLog, Append=false)
    # continue
    throw(error("No violations. This case is a waste of time"))

else
    println(string(nVol_viol_int),  " Voltage limit violations.\n",
    string(nCrnt_viol_int), " Current limit violations.\n",
    string(sum(N_violations)), " Total Violations were found.\n",
    )
    CaseLog[8] = true              # 6:"Init_PF", 7:"Viols_Exist", 8:"Opt_Commenced", 9:"Opt_Found"
    LogThis(CaseLog=CaseLog, Append=false)
    alert("JULIA: $(sum(N_violations)) Violations FOUND. Optimization will commence")
    println("Flexibility optimization will commence now @",Dates.format(now(), "HH:MM"),"\n",repeat("=",60),repeat("\n",30))
end

# break
############### Finished by baraa ################
##################################################

if !isempty(max_vol_viol_int)
    vol_viol_max_int = maximum(max_vol_viol_int)
    vol_viol_less_tol = size(findall(x -> x <= 1, max_vol_viol_int), 1)
else
    vol_viol_max_int = 0.001
    vol_viol_less_tol = 0.0
end

if !isempty(max_crnt_viol_int)
    crnt_viol_max_int = maximum(max_crnt_viol_int)
    crnt_viol_less_tol = size(findall(x -> x <= 1, max_crnt_viol_int), 1)
else
    crnt_viol_max_int = 0.001
    crnt_viol_less_tol = 0.0
end












































# Number Max Avg Above_Avg Below_Avg Less than 1%
vol_viol_nw_int =  [nVol_viol_int
                    vol_viol_max_int
                    avg_vol_viol_int
                    above_avg_vol_viol_int
                    below_avg_vol_viol_int
                    vol_viol_less_tol]

                    crnt_viol_nw_int = [nCrnt_viol_int
                    crnt_viol_max_int
                    avg_crnt_viol_int
                    above_avg_crnt_viol_int
                    below_avg_crnt_viol_int
                    crnt_viol_less_tol]

                    slack_nd = convert.(Int64, rdata_buses[idx_slack_bus_pf, 1])
slack_nd = slack_nd[1, 1]

idx_slack_from = findall(x -> x == slack_nd, rdata_lines[:, 1])
idx_slack_to = findall(x -> x == slack_nd, rdata_lines[:, 2])

num_slack_cnctd_nodes = size(idx_slack_from, 1) + size(idx_slack_to, 1)              # Number of nodes connecred to the slack = Number of times slack appear in the From and To columns
to_nodes = convert.(Int64, rdata_lines[idx_slack_from, 2])
from_nodes = convert.(Int64, rdata_lines[idx_slack_to, 1])
slack_cnctd_nodes = vcat(to_nodes, from_nodes)
slack_cnctd_nodes = vcat(slack_nd, slack_cnctd_nodes)
br_crnt_si_pf_cnv = zeros(nSc, nTP, nLines)
br_pwr_si_pf_cnv = zeros(nSc, nTP, nLines)

v_mag_pf_int = deepcopy(vol_nodes_mag_pf)
v_ang_pf_int = deepcopy(vol_nodes_theta_pf)


# (vol_viol_nw, # voltage violations of the network
# vol_viol_nw_max, #
# crnt_viol_nw, # current (amperes) violations
# crnt_viol_nw_max,
# vol_nodes_pf_cnv, # voltage magnitdues at the node after the PF problem converges (post optimization check/test)
# vol_nodes_ang_pf_cnv, # voltages at nodes, angles.
# br_crnt_si_pf_cnv, # branch current / flow
# sql_obj, # objective function value / cost funciton
# sql_time, # total time of solution
# rel_max_error, # relative max error
# rel_avg_error, # relative average error. explained in the paper [inner loop \delta S]
# norm_s_err_max_power, # normal apparent power error, max
# norm_s_err_avg, # normal apparent power error, average
# p_curt_gen, # curtailment of active power, per generator. this is an array. (gen/time/scenario)
# q_curt_gen, # curtailment of reactive power, indirect as a result of P curtailment.
# total_p_curt, # this is total curtailment. p_curt_gen is per generator. this is a scalar
# str_dis, # storage discharge, per (ESS/time/scenario), quantity is reported in KW
# str_ch,  # storage    charge, per (ESS/time/scenario), quantity is reported in KW.
# viol_nodes_new, # i dont know
# p_curt_lin, # p_curtail linear
# min_vol_limit, # minimum violation
# max_vol_limit,
# oltc_ratio, # on load tap changer ratio, per (transformer/time period/scenario)
# p_err_max_power, # not sure
# q_err_max_power,
# vm_ac, # voltage magnitude in AC. (AC means polar form, |V|/_angle)
# va_ac, # voltage angle in AC
# vol_viol_cn_milp,
# model_infs,
# br_crnt_opf_app,
# br_crnt_opf_ex,
# p_res_curt_op,
# q_res_op,
# fl_inc_op, # flexible load increment per (load unit / time period / scenario)
# fl_dec_op, # flexible load decrement
# p_ch_op,   # P charging optimum. check this against str_ch
# p_dis_op,  # P discharge optimum. check this against str_dis
# p_strg_op, # P_net = p_charge - p_disch
# p_slack_pf_ip, # slack power coming from PF solution based on IPOPT solver
# q_slack_pf_ip, # slack reactive power
# br_pwr_si_pf_cnv, # branch power, in SI units (MVA) (not in PU)
# v_rect_euler, # doesn't make sense. ignore
# v_rect_cmplx) = sql_loop(vol_viol_max,
#                         vol_viol_tol,
#                         crnt_viol_max,
#                         crnt_viol_tol,
#                         lin_itr,
#                         lin_itr_max,
#                         vol_nodes_mag_pf,
#                         vol_nodes_theta_pf,
#                         nSc,
#                         nTP,
#                         nBus,
#                         nNcurt_gen,
#                         nCurt_gen,
#                         nTrsf,
#                         nStr_active,
#                         nFl,
#                         nd_fl,
#                         flex_apc,
#                         flex_oltc,
#                         flex_adpf,
#                         flex_str,
#                         flex_fl,
#                         str_bin,
#                         fl_bin,
#                         rdata_buses,
#                         nLines,
#                         trsf_fict_nodes,
#                         tap_ratio_range,
#                         oltc_bin,
#                         prob_scs,
#                         time_step,
#                         tdata,
#                         bus_data_Ssheet,
#                         bus_data_lsheet,
#                         cost_a_str,
#                         cost_b_str,
#                         cost_c_str,
#                         cost_load_inc,
#                         cost_load_dec,
#                         nw_buses,
#                         rdata_loads,
#                         node_data,
#                         nd_curt_gen,
#                         nd_ncurt_gen,
#                         p_load,
#                         q_load,
#                         idx_Gs_lsheet,
#                         idx_Bs_lsheet,
#                         pg_max,
#                         yii_sh,
#                         i_curt_gens,
#                         sbase,
#                         dg_pf,
#                         iStr_active,
#                         yij_line,
#                         dLines,
#                         i_ncurt_gens,
#                         nw_lines,
#                         bus_data_gsheet,
#                         pg_min,
#                         qg_min,
#                         qg_max,
#                         pgen_tol,
#                         qgen_tol,
#                         p_tol,
#                         q_tol,
#                         rcvr_branch,
#                         rcvr_inj,
#                         rcvr_tso_dso,
#                         Ibase,
#                         idx_slack_bus_pf,
#                         solver_ipopt,
#                         solver_cbc,
#                         solver_bonmin,
#                         sql_itr,
#                         sql_itr_max,
#                         s_inj_error_rel,
#                         s_error_tol,
#                         angle_lb,
#                         angle_ub,
#                         cycleList,
#                         idx_slack_from,
#                         idx_slack_to,
#                         num_slack_cnctd_nodes,
#                         slack_cnctd_nodes,
#                         slack_nd,
#                         stoch_model,
#                         tratio_init,
#                         rdata_storage,
#                         load_theta,
#                         nd_Str_active,
#                         slack_bus_type,
#                         vol_ctrl_bus_type,
#                         load_bus_type,
#                         rdata_gens,
#                         nw_pPrf_data_load,
#                         nw_qPrf_data_load,
#                         nw_pPrf_data_gen_max,
#                         nw_qPrf_data_gen_max,
#                         scenario_data_p_min,
#                         scenario_data_p_max,
#                         scenario_data_q_min,
#                         scenario_data_q_max,
#                         nw_buses_pf,
#                         nw_lines_pf,
#                         nw_loads_pf,
#                         nw_gens_pf,
#                         nw_gcost_pf,
#                         nw_sbase_pf,
#                         v_initial_pf,
#                         v_mag_pf,
#                         v_angle_pf,
#                         max_mismatch_pf,
#                         epsilon,
#                         iteration,
#                         itr_max,
#                         ordata_buses_pf,
#                         ybase,
#                         vbase,
#                         I_rat,
#                         vol_viol_nw,
#                         vol_viol_nw_max,
#                         crnt_viol_nw,
#                         crnt_viol_nw_max,
#                         min_vol_limit,
#                         max_vol_limit,
#                         max_crnt_limit,
#                         term_status,
#                         vol_cstr_tol,
#                         br_crnt_si_pf_cnv,
#                         solver_cplex,
#                         int_tol,
#                         opt_tol,
#                         num_thread,
#                         rdata_trsf,
#                         nTrsf_s1,
#                         nw_trsf,
#                         n_tap,
#                         br_pwr_si_pf_cnv)

# throw(error("For no reason"))
InputDict = Dict("vol_viol_max"           =>   vol_viol_max             ,
                "vol_viol_tol"           =>   vol_viol_tol             ,
                "crnt_viol_max"          =>   crnt_viol_max            ,
                "crnt_viol_tol"          =>   crnt_viol_tol            ,
                "lin_itr"                =>   lin_itr                  ,
                "lin_itr_max"            =>   lin_itr_max              ,
                "vol_nodes_mag_pf"       =>   vol_nodes_mag_pf         ,
                "vol_nodes_theta_pf"     =>   vol_nodes_theta_pf       ,
                "nSc"                    =>   nSc                      ,
                "nTP"                    =>   nTP                      ,
                "nBus"                   =>   nBus                     ,
                "nNcurt_gen"             =>   nNcurt_gen               ,
                "nCurt_gen"              =>   nCurt_gen                ,
                "nTrsf"                  =>   nTrsf                    ,
                "nStr_active"            =>   nStr_active              ,
                "nFl"                    =>   nFl                      ,
                "nd_fl"                  =>   nd_fl                    ,
                "flex_apc"               =>   flex_apc                 ,
                "flex_oltc"              =>   flex_oltc                ,
                "flex_adpf"              =>   flex_adpf                ,
                "flex_str"               =>   flex_str                 ,
                "flex_fl"                =>   flex_fl                  ,
                "str_bin"                =>   str_bin                  ,
                "fl_bin"                 =>   fl_bin                   ,
                "rdata_buses"            =>   rdata_buses              ,
                "nLines"                 =>   nLines                   ,
                "trsf_fict_nodes"        =>   trsf_fict_nodes          ,
                "tap_ratio_range"        =>   tap_ratio_range          ,
                "oltc_bin"               =>   oltc_bin                 ,
                "prob_scs"               =>   prob_scs                 ,
                "time_step"              =>   time_step                ,
                "tdata"                  =>   tdata                    ,
                "bus_data_Ssheet"        =>   bus_data_Ssheet          ,
                "bus_data_lsheet"        =>   bus_data_lsheet          ,
                "cost_a_str"             =>   cost_a_str               ,
                "cost_b_str"             =>   cost_b_str               ,
                "cost_c_str"             =>   cost_c_str               ,
                "cost_load_inc"          =>   cost_load_inc            ,
                "cost_load_dec"          =>   cost_load_dec            ,
                "nw_buses"               =>   nw_buses                 ,
                "rdata_loads"            =>   rdata_loads              ,
                "node_data"              =>   node_data                ,
                "nd_curt_gen"            =>   nd_curt_gen              ,
                "nd_ncurt_gen"           =>   nd_ncurt_gen             ,
                "p_load"                 =>   p_load                   ,
                "q_load"                 =>   q_load                   ,
                "idx_Gs_lsheet"          =>   idx_Gs_lsheet            ,
                "idx_Bs_lsheet"          =>   idx_Bs_lsheet            ,
                "pg_max"                 =>   pg_max                   ,
                "yii_sh"                 =>   yii_sh                   ,
                "i_curt_gens"            =>   i_curt_gens              ,
                "sbase"                  =>   sbase                    ,
                "dg_pf"                  =>   dg_pf                    ,
                "iStr_active"            =>   iStr_active              ,
                "yij_line"               =>   yij_line                 ,
                "dLines"                 =>   dLines                   ,
                "i_ncurt_gens"           =>   i_ncurt_gens             ,
                "nw_lines"               =>   nw_lines                 ,
                "bus_data_gsheet"        =>   bus_data_gsheet          ,
                "pg_min"                 =>   pg_min                   ,
                "qg_min"                 =>   qg_min                   ,
                "qg_max"                 =>   qg_max                   ,
                "pgen_tol"               =>   pgen_tol                 ,
                "qgen_tol"               =>   qgen_tol                 ,
                "p_tol"                  =>   p_tol                    ,
                "q_tol"                  =>   q_tol                    ,
                "rcvr_branch"            =>   rcvr_branch              ,
                "rcvr_inj"               =>   rcvr_inj                 ,
                "rcvr_tso_dso"           =>   rcvr_tso_dso             ,
                "Ibase"                  =>   Ibase                    ,
                "idx_slack_bus_pf"       =>   idx_slack_bus_pf         ,
                "solver_ipopt"           =>   solver_ipopt             ,
                "solver_cbc"             =>   solver_cbc               ,
                "solver_bonmin"          =>   solver_bonmin            ,
                "sql_itr"                =>   sql_itr                  ,
                "sql_itr_max"            =>   sql_itr_max              ,
                "s_inj_error_rel"        =>   s_inj_error_rel          ,
                "s_error_tol"            =>   s_error_tol              ,
                "angle_lb"               =>   angle_lb                 ,
                "angle_ub"               =>   angle_ub                 ,
                "cycleList"              =>   cycleList                ,
                "idx_slack_from"         =>   idx_slack_from           ,
                "idx_slack_to"           =>   idx_slack_to             ,
                "num_slack_cnctd_nodes"  =>   num_slack_cnctd_nodes    ,
                "slack_cnctd_nodes"      =>   slack_cnctd_nodes        ,
                "slack_nd"               =>   slack_nd                 ,
                "stoch_model"            =>   stoch_model              ,
                "tratio_init"            =>   tratio_init              ,
                "rdata_storage"          =>   rdata_storage            ,
                "load_theta"             =>   load_theta               ,
                "nd_Str_active"          =>   nd_Str_active            ,
                "slack_bus_type"         =>   slack_bus_type           ,
                "vol_ctrl_bus_type"      =>   vol_ctrl_bus_type        ,
                "load_bus_type"          =>   load_bus_type            ,
                "rdata_gens"             =>   rdata_gens               ,
                "nw_pPrf_data_load"      =>   nw_pPrf_data_load        ,
                "nw_qPrf_data_load"      =>   nw_qPrf_data_load        ,
                "nw_pPrf_data_gen_max"   =>   nw_pPrf_data_gen_max     ,
                "nw_qPrf_data_gen_max"   =>   nw_qPrf_data_gen_max     ,
                "scenario_data_p_min"    =>   scenario_data_p_min      ,
                "scenario_data_p_max"    =>   scenario_data_p_max      ,
                "scenario_data_q_min"    =>   scenario_data_q_min      ,
                "scenario_data_q_max"    =>   scenario_data_q_max      ,
                "nw_buses_pf"            =>   nw_buses_pf              ,
                "nw_lines_pf"            =>   nw_lines_pf              ,
                "nw_loads_pf"            =>   nw_loads_pf              ,
                "nw_gens_pf"             =>   nw_gens_pf               ,
                "nw_gcost_pf"            =>   nw_gcost_pf              ,
                "nw_sbase_pf"            =>   nw_sbase_pf              ,
                "v_initial_pf"           =>   v_initial_pf             ,
                "v_mag_pf"               =>   v_mag_pf                 ,
                "v_angle_pf"             =>   v_angle_pf               ,
                "max_mismatch_pf"        =>   max_mismatch_pf          ,
                "epsilon"                =>   epsilon                  ,
                "iteration"              =>   iteration                ,
                "itr_max"                =>   itr_max                  ,
                "ordata_buses_pf"        =>   ordata_buses_pf          ,
                "ybase"                  =>   ybase                    ,
                "vbase"                  =>   vbase                    ,
                "I_rat"                  =>   I_rat                    ,
                "vol_viol_nw"            =>   vol_viol_nw              ,
                "vol_viol_nw_max"        =>   vol_viol_nw_max          ,
                "crnt_viol_nw"           =>   crnt_viol_nw             ,
                "crnt_viol_nw_max"       =>   crnt_viol_nw_max         ,
                "min_vol_limit"          =>   min_vol_limit            ,
                "max_vol_limit"          =>   max_vol_limit            ,
                "max_crnt_limit"         =>   max_crnt_limit           ,
                "term_status"            =>   term_status              ,
                "vol_cstr_tol"           =>   vol_cstr_tol             ,
                "br_crnt_si_pf_cnv"      =>   br_crnt_si_pf_cnv        ,
                "solver_cplex"           =>   solver_cplex             ,
                "int_tol"                =>   int_tol                  ,
                "opt_tol"                =>   opt_tol                  ,
                "num_thread"             =>   num_thread               ,
                "rdata_trsf"             =>   rdata_trsf               ,
                "nTrsf_s1"               =>   nTrsf_s1                 ,
                "nw_trsf"                =>   nw_trsf                  ,
                "n_tap"                  =>   n_tap                    ,
                "br_pwr_si_pf_cnv"       =>   br_pwr_si_pf_cnv         ,
                "soc_0"                  =>   soc_0,
                )

# let OutDict = Dict()
# try
SQL_Output  = sql_loop(InputDict);
OutDict     = SQL_Output;
# global OutDict     = SQL_Output

println(string("Exited SQL_Loop. with exitflag (term_status) = ", OutDict["term_status"]))
# catch Err
#     if Err == UndefVarError(:vm_ac)
#         println("Case is infeasible. Will proceed to next iteration")
#         # continue
#     else
        # @infiltrate
        # println("Error details are: ", string(Err))
        # rethrow(Err)
        # CaseLog[9] = 0.5  ### 6:"Init_PF", 7:"Viols_Exist", 8:"Opt_Commenced", 9:"Opt_Found"
        # LogThis(CaseLog=CaseLog, Append=false)
        # @printf("---===<<<\tError in:#%3i/%i \"%s\"\t>>>===---\n", i, N_Cases, TheKey)
#     end
# end # END OF TRY
# end # END OF LET OUTDICT=....

# println("Outside Try-Break", string(OutDict))

# vol_viol_nw              =    OutDict["vol_viol_nw"]
# vol_viol_nw_max          =    OutDict["vol_viol_nw_max"]
# crnt_viol_nw             =    OutDict["crnt_viol_nw"]
# crnt_viol_nw_max         =    OutDict["crnt_viol_nw_max"]
# vol_nodes_pf_cnv         =    OutDict["vol_nodes_pf_cnv"]
# vol_nodes_ang_pf_cnv     =    OutDict["vol_nodes_ang_pf_cnv"]
# br_crnt_si_pf_cnv        =    OutDict["br_crnt_si_pf_cnv"]
# sql_obj                  =    OutDict["sql_obj"]
# sql_time                 =    OutDict["sql_time"]
# rel_max_error            =    OutDict["rel_max_error"]
# rel_avg_error            =    OutDict["rel_avg_error"]
# norm_s_err_max_power     =    OutDict["norm_s_err_max_power"]
# norm_s_err_avg           =    OutDict["norm_s_err_avg"]
# p_curt_gen               =    OutDict["p_curt_gen"]
# q_curt_gen               =    OutDict["q_curt_gen"]
# total_p_curt             =    OutDict["total_p_curt"]
# str_dis                  =    OutDict["str_dis"]
# str_ch                   =    OutDict["str_ch"]
# viol_nodes_new           =    OutDict["viol_nodes_new"]
# p_curt_lin               =    OutDict["p_curt_lin"]
# min_vol_limit            =    OutDict["min_vol_limit"]
# max_vol_limit            =    OutDict["max_vol_limit"]
# oltc_ratio               =    OutDict["oltc_ratio"]
# p_err_max_power          =    OutDict["p_err_max_power"]
# q_err_max_power          =    OutDict["q_err_max_power"]
# vm_ac                    =    OutDict["vm_ac"]
# va_ac                    =    OutDict["va_ac"]
# vol_viol_cn_milp         =    OutDict["vol_viol_cn_milp"]
# model_infs               =    OutDict["model_infs"]
# br_crnt_opf_app          =    OutDict["br_crnt_opf_app"]
# br_crnt_opf_ex           =    OutDict["br_crnt_opf_ex"]
# p_res_curt_op            =    OutDict["p_res_curt_op"]
# q_res_op                 =    OutDict["q_res_op"]
# fl_inc_op                =    OutDict["fl_inc_op"]
# fl_dec_op                =    OutDict["fl_dec_op"]
# p_ch_op                  =    OutDict["p_ch_op"]
# p_dis_op                 =    OutDict["p_dis_op"]
# p_strg_op                =    OutDict["p_strg_op"]
# p_slack_pf_ip            =    OutDict["p_slack_pf_ip"]
# q_slack_pf_ip            =    OutDict["q_slack_pf_ip"]
# br_pwr_si_pf_cnv         =    OutDict["br_pwr_si_pf_cnv"]
# v_rect_euler             =    OutDict["v_rect_euler"]
# v_rect_cmplx             =    OutDict["v_rect_cmplx"]

# break
# After this point, all code is just formatting the results, presentation
################################################################################
###################### Output Data to Comillas Tools ###########################
################################################################################

# severity = br_crnt_si_pf_cnv./max_crnt_limit_si
# severity = maximum(permutedims(reshape(severity,nTP,nLines)),dims=2)
#
# br_pwr = permutedims(reshape(br_pwr_si_pf_cnv,nTP,nLines))
# p_pwr = real(br_pwr)
# q_pwr = imag(br_pwr)

################################################################################
######################### Output Data to ICENT Tools ###########################
################################################################################

### ---- Curtailed Active Power and Output Reactive Power from RES Data -----###
nodes_curt_gen = Int64.(rdata_gens[idx_curt_dg, 1])
avg_p_gen_curt = reshape(sum(prob_scs .* OutDict["p_res_curt_op"], dims=1), (nTP, nCurt_gen))
avg_q_gen_curt = reshape(sum(prob_scs .* OutDict["q_res_op"], dims=1), (nTP, nCurt_gen))

avg_p_gen_curt = permutedims(avg_p_gen_curt)
avg_q_gen_curt = permutedims(avg_q_gen_curt)

# avg_p_gen_curt = hcat(nodes_curt_gen, avg_p_gen_curt)
# avg_q_gen_curt = hcat(nodes_curt_gen, avg_q_gen_curt)

# ################################################################################
# ##### ----------------- Flexible Load Data --------------------------------#####
# ################################################################################
nodes_fl = Int64.(rdata_loads[iFl, 1])
avg_fl_inc = prob_scs .* OutDict["fl_inc_op"]
avg_fl_dec = prob_scs .* OutDict["fl_dec_op"]

avg_fl_inc_op = permutedims(avg_fl_inc[sc_max_prob, :, :])
avg_fl_dec_op = permutedims(avg_fl_dec[sc_max_prob, :, :])

# avg_fl_inc_op = hcat(nodes_fl, avg_fl_inc_op)
# avg_fl_dec_op = hcat(nodes_fl, avg_fl_dec_op)

# ################################################################################
# ##### ----------------- Energy Storage Data --------------------------------####
# ################################################################################
nodes_str = Int64.(rdata_storage[iStr_active, 1])
avg_p_ch_nd = prob_scs .* OutDict["p_ch_op"]
avg_p_dch_nd = prob_scs .* OutDict["p_dis_op"]

avg_p_ch_op  = permutedims(avg_p_ch_nd[sc_max_prob, :, :])
avg_p_dch_op = permutedims(avg_p_dch_nd[sc_max_prob, :, :])

#added by baraa. The multiplication by prob_scs applies on each scenario separately. scenarios are NOT summed after that.
#But i dont know why he is taking th scenario with highest probability: sc_max_prob
# avg_str_soc_nd   = prob_scs .* OutDict["str_soc"]
avg_str_soc_nd   = OutDict["str_soc"]
avg_str_soc_op   = permutedims(avg_str_soc_nd[sc_max_prob, :, :])

# avg_p_ch_op = hcat(nodes_str, avg_p_ch_op)
# avg_p_dch_op = hcat(nodes_str, avg_p_dch_op)
# ###############################################################################
#
header_op = [string("t", i) for i in 1:24]
# header_op[1] = "nodes"
# header_op = ["nodes", "t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10", "t11",
#     "t12", "t13", "t14", "t15", "t16", "t17", "t18", "t19", "t20", "t21", "t22", "t23", "t24"]

nodes_curt_gen_df = DataFrames.DataFrame(nodes=nodes_curt_gen)

df_p_curt = DataFrame(avg_p_gen_curt, header_op)
df_p_curt = hcat(nodes_curt_gen_df,df_p_curt)

df_q_res = DataFrame(avg_q_gen_curt, header_op)
df_q_res = hcat(nodes_curt_gen_df,df_q_res)
#---
nodes_curt_fl_df = DataFrames.DataFrame(nodes=nodes_fl)

df_fl_inc = DataFrame(avg_fl_inc_op, header_op)
df_fl_inc = hcat(nodes_curt_fl_df,df_fl_inc)

df_fl_dec = DataFrame(avg_fl_dec_op, header_op)
df_fl_dec = hcat(nodes_curt_fl_df,df_fl_dec)
#---
nodes_curt_str_df = DataFrames.DataFrame(nodes=nodes_str)

df_str_ch = DataFrame(avg_p_ch_op, header_op)
df_str_ch = hcat(nodes_curt_str_df,df_str_ch)

df_str_dch = DataFrame(avg_p_dch_op, header_op)
df_str_dch = hcat(nodes_curt_str_df,df_str_dch)

df_str_soc = DataFrame(avg_str_soc_op, header_op)
df_str_soc = hcat(nodes_curt_str_df,df_str_soc)

# ################################################################################
# ################## Cost of Procurement of Flexibility ##########################
# ################################################################################
header_cost = ["cost_flex"]
data_obj = reshape(OutDict["sql_obj"], (sql_itr_max * lin_itr_max))
idx_0 = findall(x -> x != 0, data_obj)
if sizeof(idx_0)>0
    cost_flex_proc = [data_obj[idx_0[end]]]
else
    # cost_flex_proc = NaN # a value of NaN or INF cause problems to EXCEL which is expecting a finite number
    cost_flex_proc =data_obj[1]
end
df_cost_pr = DataFrame(cost=cost_flex_proc)
#
# ################################################################################
# ####------------- Writing Output Data to Excel File ---------------------- #####
# ################################################################################

function AppendStats(DF, IncludeStats=false)
    if (size(DF,1) > 1) & (IncludeStats==true)
        newDF = deepcopy(DF)
        push!(newDF,
            append!([-Inf], [minimum(DF[:,i]) for i in 2:size(DF,2)]), # -INF refers to min_norm
            append!([ Inf], [maximum(DF[:,i]) for i in 2:size(DF,2)]), # INF refers to max_norm
            append!([ Inf], [   mean(DF[:,i]) for i in 2:size(DF,2)]), # 0.5 refers to avg
            append!([ Inf], [    sum(DF[:,i]) for i in 2:size(DF,2)])) # 707 looks like TOT, or you could just use 1 to refer to 1st norm
        newDF.nodes = [string(i) for i in newDF.nodes]
        newDF.nodes[size(DF,1)+1] = "Min"
        newDF.nodes[size(DF,1)+2] = "Max"
        newDF.nodes[size(DF,1)+3] = "Avg"
        newDF.nodes[size(DF,1)+4] = "Sum"
        return newDF
    else
        return DF
    end
end

println("Finished running && Writing to Excel @", Dates.format(now(), "HH:MM"))
XLSX.openxlsx(file_op, mode="rw") do xf
    Append2Excel(xf, AppendStats(df_p_curt,IncludeStats), "APC_MW")
    Append2Excel(xf, AppendStats(df_str_ch,IncludeStats), "EES_CH_MW")
    Append2Excel(xf, AppendStats(df_str_dch,IncludeStats),"EES_DCH_MW")
    Append2Excel(xf, AppendStats(df_str_soc,IncludeStats),"EES_SoC_MWh")
    Append2Excel(xf, AppendStats(df_fl_inc,IncludeStats), "FL_OD_MW")
    Append2Excel(xf, AppendStats(df_fl_dec,IncludeStats), "FL_UD_MW")
    Append2Excel(xf, df_cost_pr, "COST")
end

if OutDict["term_status"] == 1
    MILP_term_status = "Success"
elseif OutDict["term_status"] == 3
    MILP_term_status = "INFEASIBLE"
else
    ErrMsg = string("OutDict[\"term_status\"] = ",OutDict["term_status"],", OutDict[\"term_status\"] can either be 1 or 3")
    throw(error(ErrMsg))
    MILP_term_status = "INFEASIBLE"
end

if any(i-> i!=0, avg_p_ch_op+avg_p_dch_op)
    alert("Storage action detected. look at the results")
    # using WAV
    y, fs = wavread(raw"C:\Windows\Media\Ring04.wav")
    # y, fs = wavread(raw"C:\Windows\Media\tada.wav")
    wavplay(y, fs)
end
if cost_flex_proc[1]>5e-2
    alert("Non-trivial solution found")
    # using WAV
    y, fs = wavread(raw"C:\Windows\Media\tada.wav")
    # y, fs = wavread(raw"C:\Windows\Media\tada.wav")
    wavplay(y, fs)
end
bool_trivial = (cost_flex_proc[1]>5e-2) ? "non-trivial"  : "trivial"
bool_trivial = (MILP_term_status=="INFEASIBLE") ? "" : bool_trivial

CaseLog[9] = string(MILP_term_status,(sizeof(bool_trivial)>0 ? " - " : ""),bool_trivial)  # 6:"Init_PF", 7:"Viols_Exist", 8:"Opt_Commenced", 9:"Opt_Found"
LogThis(CaseLog=CaseLog, Append=false)
alert(string("Finished optimization. Results are ",bool_trivial," and were exported to Excel: ", TheKey))
println("Finished writing final results to excel @ ",Dates.format(now(), "HH:MM"),".\nFind final results @ ",file_op)
println("Results are ", bool_trivial)
# run(`$(pwd())\\$file_op`) # doesnt work because file_op is a relative path, not the complete path with drive (e.g. C:\))
# zft=open(file_op)
# continue
# XLSX.writetable(file_op,
#                 APC_MW      = (collect(DataFrames.eachcol(df_p_curt)),  DataFrames.names(df_p_curt)),
#                 EES_CH_MW   = (collect(DataFrames.eachcol(df_str_ch)),  names(df_str_ch)),
#                 EES_DCH_MW  = (collect(DataFrames.eachcol(df_str_dch)), names(df_str_dch)),
#                 FL_OD_MW    = (collect(DataFrames.eachcol(df_fl_inc)),  names(df_fl_inc)),
#                 FL_UD_MW    = (collect(DataFrames.eachcol(df_fl_dec)),  names(df_fl_dec)),
#                 COST        = (collect(DataFrames.eachcol(df_cost_pr)), names(df_cost_pr)),
#                 ; overwrite=true)

# ################################################################################
# end
# end
