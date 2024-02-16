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
using DataFrames
using Dates
using Ipopt
using JuMP
using Juniper
using LinearAlgebra
using MathOptInterface
using OdsIO
using Plots
using XLSX
using Statistics

using Infiltrator
using Alert
using Printf
# using DefaultApplication

# To remove the warnings from the log, use REGEXP in notepad++ to remove the following: ^┌.+\r\n(.+\r\n)?└ @.+$\r\n

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
# CasesList = Dict([
#                   ("UK_DX_30", "input_data/CIRED/UK_Dx_01_2030.ods"),
#                   ("UK_DX_40", "input_data/CIRED/UK_Dx_01_2040.ods"),
#                   ("UK_DX_50", "input_data/CIRED/UK_Dx_01_2050.ods"),

#                   ("PT_DX_30", "input_data/CIRED/PT_Dx_01_2030.ods"),
#                   ("PT_DX_40", "input_data/CIRED/PT_Dx_01_2040.ods"),
#                   ("PT_DX_50", "input_data/CIRED/PT_Dx_01_2050.ods"),

#                   ("ES_DX_30", "input_data/CIRED/ES_Dx_03_2030.ods"),
#                   ("ES_DX_40", "input_data/CIRED/ES_Dx_03_2040.ods"),
#                   ("ES_DX_50", "input_data/CIRED/ES_Dx_03_2050.ods"),

#                   ("HR_DX_30:Br", "input_data/CIRED/HR_Dx_01_2030_BROWN.ods"),
#                   ("HR_DX_40:Br", "input_data/CIRED/HR_Dx_01_2040_BROWN.ods"),
#                   ("HR_DX_50:Br", "input_data/CIRED/HR_Dx_01_2050_BROWN.ods"),

#                   ("HR_DX_30:G", "input_data/CIRED/HR_Dx_01_2030_GREEN.ods"),
#                   ("HR_DX_40:G", "input_data/CIRED/HR_Dx_01_2040_GREEN.ods"),
#                   ("HR_DX_50:G", "input_data/CIRED/HR_Dx_01_2050_GREEN.ods"),

#                   ("HR_DX_30:R", "input_data/CIRED/HR_Dx_01_2030_RED.ods"),
#                   ("HR_DX_40:R", "input_data/CIRED/HR_Dx_01_2040_RED.ods"),
#                   ("HR_DX_50:R", "input_data/CIRED/HR_Dx_01_2050_RED.ods"),

#                   ("UK"      , "input_data/uk_dx_01_2020.ods"),
#                   ("PT"      , "input_data/pt_dx_01_2020.ods"),
#                   ("ES"      , "input_data/es_dx_01_2020.ods"),
#                   ("HR:Br"   , "input_data/hr_dx_01_2020_brown.ods"),
#                   ("HR:G"    , "input_data/hr_dx_01_2020_green.ods"),
#                   ("HR:R"    , "input_data/hr_dx_01_2020_red.ods"),

#                   ("WP5:HR:K", "input_data/hr_dx_01_2020_brown.ods"),
#                   ("WP5:HR:R", "input_data/WP5/hr_kpc_10_red.ods"),
#                   ("WP5:HR:B", "input_data/WP5/hr_kpc_10_blue.ods"),
#                   ("WP5:HR:G", "input_data/WP5/hr_kpc_10_green.ods"),
#                   ])

CasesList =    Dict([
    ("PT_Dx_01_2030_A_Bd_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_A_Bd_WithFlex.ods" ),
    ("PT_Dx_01_2030_A_Bd_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_A_Bd_WOFlex.ods"   ),
    ("PT_Dx_01_2030_A_Sa_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_A_Sa_WithFlex.ods" ),
    ("PT_Dx_01_2030_A_Sa_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_A_Sa_WOFlex.ods"   ),
    ("PT_Dx_01_2030_A_Su_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_A_Su_WithFlex.ods" ),
    ("PT_Dx_01_2030_A_Su_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_A_Su_WOFlex.ods"   ),

    ("PT_Dx_01_2030_Su_Bd_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_Su_Bd_WithFlex.ods"),
    ("PT_Dx_01_2030_Su_Bd_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_Su_Bd_WOFlex.ods"  ),
    ("PT_Dx_01_2030_Su_Sa_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_Su_Sa_WithFlex.ods"),
    ("PT_Dx_01_2030_Su_Sa_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_Su_Sa_WOFlex.ods"  ),
    ("PT_Dx_01_2030_Su_Su_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_Su_Su_WithFlex.ods"),
    ("PT_Dx_01_2030_Su_Su_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_Su_Su_WOFlex.ods"  ),

    ("PT_Dx_01_2030_S_Bd_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_S_Bd_WithFlex.ods" ),
    ("PT_Dx_01_2030_S_Bd_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_S_Bd_WOFlex.ods"   ),
    ("PT_Dx_01_2030_S_Sa_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_S_Sa_WithFlex.ods" ),
    ("PT_Dx_01_2030_S_Sa_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_S_Sa_WOFlex.ods"   ),
    ("PT_Dx_01_2030_S_Su_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_S_Su_WithFlex.ods" ),
    ("PT_Dx_01_2030_S_Su_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_S_Su_WOFlex.ods"   ),

    ("PT_Dx_01_2030_W_Bd_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_W_Bd_WithFlex.ods" ),
    ("PT_Dx_01_2030_W_Bd_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_W_Bd_WOFlex.ods"   ),
    ("PT_Dx_01_2030_W_Sa_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_W_Sa_WithFlex.ods" ),
    ("PT_Dx_01_2030_W_Sa_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_W_Sa_WOFlex.ods"   ),
    ("PT_Dx_01_2030_W_Su_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_W_Su_WithFlex.ods" ),
    ("PT_Dx_01_2030_W_Su_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2030_W_Su_WOFlex.ods"   ),

    ("PT_Dx_01_2040_A_Bd_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_A_Bd_WithFlex.ods" ),
    ("PT_Dx_01_2040_A_Bd_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_A_Bd_WOFlex.ods"   ),
    ("PT_Dx_01_2040_A_Sa_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_A_Sa_WithFlex.ods" ),
    ("PT_Dx_01_2040_A_Sa_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_A_Sa_WOFlex.ods"   ),
    ("PT_Dx_01_2040_A_Su_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_A_Su_WithFlex.ods" ),
    ("PT_Dx_01_2040_A_Su_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_A_Su_WOFlex.ods"   ),

    ("PT_Dx_01_2040_Su_Bd_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_Su_Bd_WithFlex.ods"),
    ("PT_Dx_01_2040_Su_Bd_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_Su_Bd_WOFlex.ods"  ),
    ("PT_Dx_01_2040_Su_Sa_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_Su_Sa_WithFlex.ods"),
    ("PT_Dx_01_2040_Su_Sa_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_Su_Sa_WOFlex.ods"  ),
    ("PT_Dx_01_2040_Su_Su_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_Su_Su_WithFlex.ods"),
    ("PT_Dx_01_2040_Su_Su_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_Su_Su_WOFlex.ods"  ),

    ("PT_Dx_01_2040_S_Bd_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_S_Bd_WithFlex.ods" ),
    ("PT_Dx_01_2040_S_Bd_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_S_Bd_WOFlex.ods"   ),
    ("PT_Dx_01_2040_S_Sa_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_S_Sa_WithFlex.ods" ),
    ("PT_Dx_01_2040_S_Sa_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_S_Sa_WOFlex.ods"   ),
    ("PT_Dx_01_2040_S_Su_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_S_Su_WithFlex.ods" ),
    ("PT_Dx_01_2040_S_Su_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_S_Su_WOFlex.ods"   ),

    ("PT_Dx_01_2040_W_Bd_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_W_Bd_WithFlex.ods" ),
    ("PT_Dx_01_2040_W_Bd_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_W_Bd_WOFlex.ods"   ),
    ("PT_Dx_01_2040_W_Sa_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_W_Sa_WithFlex.ods" ),
    ("PT_Dx_01_2040_W_Sa_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_W_Sa_WOFlex.ods"   ),
    ("PT_Dx_01_2040_W_Su_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_W_Su_WithFlex.ods" ),
    ("PT_Dx_01_2040_W_Su_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2040_W_Su_WOFlex.ods"   ),

    ("PT_Dx_01_2050_A_Bd_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_A_Bd_WithFlex.ods" ),
    ("PT_Dx_01_2050_A_Bd_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_A_Bd_WOFlex.ods"   ),
    ("PT_Dx_01_2050_A_Sa_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_A_Sa_WithFlex.ods" ),
    ("PT_Dx_01_2050_A_Sa_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_A_Sa_WOFlex.ods"   ),
    ("PT_Dx_01_2050_A_Su_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_A_Su_WithFlex.ods" ),
    ("PT_Dx_01_2050_A_Su_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_A_Su_WOFlex.ods"   ),

    ("PT_Dx_01_2050_Su_Bd_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_Su_Bd_WithFlex.ods"),
    ("PT_Dx_01_2050_Su_Bd_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_Su_Bd_WOFlex.ods"  ),
    ("PT_Dx_01_2050_Su_Sa_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_Su_Sa_WithFlex.ods"),
    ("PT_Dx_01_2050_Su_Sa_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_Su_Sa_WOFlex.ods"  ),
    ("PT_Dx_01_2050_Su_Su_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_Su_Su_WithFlex.ods"),
    ("PT_Dx_01_2050_Su_Su_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_Su_Su_WOFlex.ods"  ),

    ("PT_Dx_01_2050_S_Bd_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_S_Bd_WithFlex.ods" ),
    ("PT_Dx_01_2050_S_Bd_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_S_Bd_WOFlex.ods"   ),
    ("PT_Dx_01_2050_S_Sa_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_S_Sa_WithFlex.ods" ),
    ("PT_Dx_01_2050_S_Sa_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_S_Sa_WOFlex.ods"   ),
    ("PT_Dx_01_2050_S_Su_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_S_Su_WithFlex.ods" ),
    ("PT_Dx_01_2050_S_Su_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_S_Su_WOFlex.ods"   ),

    ("PT_Dx_01_2050_W_Bd_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_W_Bd_WithFlex.ods" ),
    ("PT_Dx_01_2050_W_Bd_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_W_Bd_WOFlex.ods"   ),
    ("PT_Dx_01_2050_W_Sa_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_W_Sa_WithFlex.ods" ),
    ("PT_Dx_01_2050_W_Sa_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_W_Sa_WOFlex.ods"   ),
    ("PT_Dx_01_2050_W_Su_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_W_Su_WithFlex.ods" ),
    ("PT_Dx_01_2050_W_Su_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/PT DX 01/PT_Dx_01_2050_W_Su_WOFlex.ods"   ),

    ("UK_DX_01_BaseCase"                                                 , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/BaseCase.ods"                                                          ),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_S_WithFlex"            , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2030_Distribution_Network_Urban_UK_S_WithFlex.ods"            ),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_S_WOFlex"              , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2030_Distribution_Network_Urban_UK_S_WOFlex.ods"              ),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_S_Without_1_7_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2030_Distribution_Network_Urban_UK_S_Without_1_7_WithFlex.ods"),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_S_Without_1_7_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2030_Distribution_Network_Urban_UK_S_Without_1_7_WOFlex.ods"  ),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_S_Without_2_7_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2030_Distribution_Network_Urban_UK_S_Without_2_7_WithFlex.ods"),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_S_Without_2_7_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2030_Distribution_Network_Urban_UK_S_Without_2_7_WOFlex.ods"  ),

    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_W_WithFlex"            , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2030_Distribution_Network_Urban_UK_W_WithFlex.ods"            ),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_W_WOFlex"              , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2030_Distribution_Network_Urban_UK_W_WOFlex.ods"              ),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_W_Without_1_7_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2030_Distribution_Network_Urban_UK_W_Without_1_7_WithFlex.ods"),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_W_Without_1_7_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2030_Distribution_Network_Urban_UK_W_Without_1_7_WOFlex.ods"  ),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_W_Without_2_7_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2030_Distribution_Network_Urban_UK_W_Without_2_7_WithFlex.ods"),
    ("UK_Dx_01_2030_Distribution_Network_Urban_UK_W_Without_2_7_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2030_Distribution_Network_Urban_UK_W_Without_2_7_WOFlex.ods"  ),

    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_S_WithFlex"            , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2040_Distribution_Network_Urban_UK_S_WithFlex.ods"            ),
    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_S_WOFlex"              , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2040_Distribution_Network_Urban_UK_S_WOFlex.ods"              ),
    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_S_Without_1_7_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2040_Distribution_Network_Urban_UK_S_Without_1_7_WithFlex.ods"),
    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_S_Without_1_7_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2040_Distribution_Network_Urban_UK_S_Without_1_7_WOFlex.ods"  ),
    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_S_Without_2_7_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2040_Distribution_Network_Urban_UK_S_Without_2_7_WithFlex.ods"),
    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_S_Without_2_7_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2040_Distribution_Network_Urban_UK_S_Without_2_7_WOFlex.ods"  ),

    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_W_WithFlex"            , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2040_Distribution_Network_Urban_UK_W_WithFlex.ods"            ),
    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_W_WOFlex"              , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2040_Distribution_Network_Urban_UK_W_WOFlex.ods"              ),
    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_W_Without_1_7_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2040_Distribution_Network_Urban_UK_W_Without_1_7_WithFlex.ods"),
    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_W_Without_1_7_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2040_Distribution_Network_Urban_UK_W_Without_1_7_WOFlex.ods"  ),
    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_W_Without_2_7_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2040_Distribution_Network_Urban_UK_W_Without_2_7_WithFlex.ods"),
    ("UK_Dx_01_2040_Distribution_Network_Urban_UK_W_Without_2_7_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2040_Distribution_Network_Urban_UK_W_Without_2_7_WOFlex.ods"  ),

    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_S_WithFlex"            , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2050_Distribution_Network_Urban_UK_S_WithFlex.ods"            ),
    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_S_WOFlex"              , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2050_Distribution_Network_Urban_UK_S_WOFlex.ods"              ),
    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_S_Without_1_7_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2050_Distribution_Network_Urban_UK_S_Without_1_7_WithFlex.ods"),
    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_S_Without_1_7_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2050_Distribution_Network_Urban_UK_S_Without_1_7_WOFlex.ods"  ),
    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_S_Without_2_7_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2050_Distribution_Network_Urban_UK_S_Without_2_7_WithFlex.ods"),
    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_S_Without_2_7_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2050_Distribution_Network_Urban_UK_S_Without_2_7_WOFlex.ods"  ),

    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_W_WithFlex"            , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2050_Distribution_Network_Urban_UK_W_WithFlex.ods"            ),
    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_W_WOFlex"              , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2050_Distribution_Network_Urban_UK_W_WOFlex.ods"              ),
    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_W_Without_1_7_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2050_Distribution_Network_Urban_UK_W_Without_1_7_WithFlex.ods"),
    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_W_Without_1_7_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2050_Distribution_Network_Urban_UK_W_Without_1_7_WOFlex.ods"  ),
    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_W_Without_2_7_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2050_Distribution_Network_Urban_UK_W_Without_2_7_WithFlex.ods"),
    ("UK_Dx_01_2050_Distribution_Network_Urban_UK_W_Without_2_7_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/UK DX 01/UK_Dx_01_2050_Distribution_Network_Urban_UK_W_Without_2_7_WOFlex.ods"  ),

    ("HR_DX_01_BROWN_BaseCase"           , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/BaseCase.ods"                                ),
    ("HR_DX_01_RED_BaseCase"             , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/BaseCase.ods"                                  ),
    ("HR_DX_01_GREEN_BaseCase"           , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/BaseCase.ods"                                ),

    ("HR_Dx_01_2030_A_Bd_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_A_Bd_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2030_A_Bd_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_A_Bd_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2030_A_Sa_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_A_Sa_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2030_A_Sa_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_A_Sa_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2030_A_Su_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_A_Su_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2030_A_Su_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_A_Su_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2030_Su_Bd_brown_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_Su_Bd_brown_WithFlex_BROWN.ods"),
    ("HR_Dx_01_2030_Su_Bd_brown_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_Su_Bd_brown_WOFlex_BROWN.ods"  ),
    ("HR_Dx_01_2030_Su_Sa_brown_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_Su_Sa_brown_WithFlex_BROWN.ods"),
    ("HR_Dx_01_2030_Su_Sa_brown_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_Su_Sa_brown_WOFlex_BROWN.ods"  ),
    ("HR_Dx_01_2030_Su_Su_brown_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_Su_Su_brown_WithFlex_BROWN.ods"),
    ("HR_Dx_01_2030_Su_Su_brown_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_Su_Su_brown_WOFlex_BROWN.ods"  ),
    ("HR_Dx_01_2030_S_Bd_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_S_Bd_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2030_S_Bd_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_S_Bd_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2030_S_Sa_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_S_Sa_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2030_S_Sa_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_S_Sa_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2030_S_Su_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_S_Su_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2030_S_Su_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_S_Su_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2030_W_Bd_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_W_Bd_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2030_W_Bd_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_W_Bd_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2030_W_Sa_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_W_Sa_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2030_W_Sa_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_W_Sa_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2030_W_Su_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_W_Su_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2030_W_Su_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2030_W_Su_brown_WOFlex_BROWN.ods"   ),

    ("HR_Dx_01_2030_A_Bd_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_A_Bd_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2030_A_Bd_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_A_Bd_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2030_A_Sa_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_A_Sa_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2030_A_Sa_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_A_Sa_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2030_A_Su_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_A_Su_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2030_A_Su_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_A_Su_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2030_Su_Bd_red_WithFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_Su_Bd_red_WithFlex_RED.ods"      ),
    ("HR_Dx_01_2030_Su_Bd_red_WOFlex"    , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_Su_Bd_red_WOFlex_RED.ods"        ),
    ("HR_Dx_01_2030_Su_Sa_red_WithFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_Su_Sa_red_WithFlex_RED.ods"      ),
    ("HR_Dx_01_2030_Su_Sa_red_WOFlex"    , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_Su_Sa_red_WOFlex_RED.ods"        ),
    ("HR_Dx_01_2030_Su_Su_red_WithFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_Su_Su_red_WithFlex_RED.ods"      ),
    ("HR_Dx_01_2030_Su_Su_red_WOFlex"    , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_Su_Su_red_WOFlex_RED.ods"        ),
    ("HR_Dx_01_2030_S_Bd_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_S_Bd_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2030_S_Bd_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_S_Bd_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2030_S_Sa_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_S_Sa_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2030_S_Sa_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_S_Sa_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2030_S_Su_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_S_Su_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2030_S_Su_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_S_Su_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2030_W_Bd_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_W_Bd_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2030_W_Bd_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_W_Bd_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2030_W_Sa_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_W_Sa_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2030_W_Sa_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_W_Sa_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2030_W_Su_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_W_Su_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2030_W_Su_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2030_W_Su_red_WOFlex_RED.ods"         ),

    ("HR_Dx_01_2030_A_Bd_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_A_Bd_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2030_A_Bd_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_A_Bd_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2030_A_Sa_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_A_Sa_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2030_A_Sa_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_A_Sa_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2030_A_Su_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_A_Su_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2030_A_Su_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_A_Su_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2030_Su_Bd_green_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_Su_Bd_green_WithFlex_GREEN.ods"),
    ("HR_Dx_01_2030_Su_Bd_green_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_Su_Bd_green_WOFlex_GREEN.ods"  ),
    ("HR_Dx_01_2030_Su_Sa_green_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_Su_Sa_green_WithFlex_GREEN.ods"),
    ("HR_Dx_01_2030_Su_Sa_green_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_Su_Sa_green_WOFlex_GREEN.ods"  ),
    ("HR_Dx_01_2030_Su_Su_green_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_Su_Su_green_WithFlex_GREEN.ods"),
    ("HR_Dx_01_2030_Su_Su_green_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_Su_Su_green_WOFlex_GREEN.ods"  ),
    ("HR_Dx_01_2030_S_Bd_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_S_Bd_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2030_S_Bd_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_S_Bd_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2030_S_Sa_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_S_Sa_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2030_S_Sa_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_S_Sa_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2030_S_Su_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_S_Su_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2030_S_Su_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_S_Su_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2030_W_Bd_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_W_Bd_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2030_W_Bd_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_W_Bd_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2030_W_Sa_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_W_Sa_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2030_W_Sa_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_W_Sa_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2030_W_Su_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_W_Su_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2030_W_Su_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2030_W_Su_green_WOFlex_GREEN.ods"   ),

    ("HR_Dx_01_2040_A_Bd_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_A_Bd_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2040_A_Bd_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_A_Bd_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2040_A_Sa_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_A_Sa_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2040_A_Sa_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_A_Sa_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2040_A_Su_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_A_Su_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2040_A_Su_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_A_Su_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2040_Su_Bd_brown_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_Su_Bd_brown_WithFlex_BROWN.ods"),
    ("HR_Dx_01_2040_Su_Bd_brown_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_Su_Bd_brown_WOFlex_BROWN.ods"  ),
    ("HR_Dx_01_2040_Su_Sa_brown_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_Su_Sa_brown_WithFlex_BROWN.ods"),
    ("HR_Dx_01_2040_Su_Sa_brown_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_Su_Sa_brown_WOFlex_BROWN.ods"  ),
    ("HR_Dx_01_2040_Su_Su_brown_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_Su_Su_brown_WithFlex_BROWN.ods"),
    ("HR_Dx_01_2040_Su_Su_brown_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_Su_Su_brown_WOFlex_BROWN.ods"  ),
    ("HR_Dx_01_2040_S_Bd_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_S_Bd_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2040_S_Bd_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_S_Bd_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2040_S_Sa_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_S_Sa_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2040_S_Sa_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_S_Sa_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2040_S_Su_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_S_Su_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2040_S_Su_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_S_Su_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2040_W_Bd_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_W_Bd_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2040_W_Bd_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_W_Bd_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2040_W_Sa_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_W_Sa_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2040_W_Sa_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_W_Sa_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2040_W_Su_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_W_Su_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2040_W_Su_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2040_W_Su_brown_WOFlex_BROWN.ods"   ),

    ("HR_Dx_01_2040_A_Bd_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_A_Bd_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2040_A_Bd_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_A_Bd_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2040_A_Sa_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_A_Sa_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2040_A_Sa_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_A_Sa_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2040_A_Su_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_A_Su_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2040_A_Su_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_A_Su_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2040_Su_Bd_red_WithFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_Su_Bd_red_WithFlex_RED.ods"      ),
    ("HR_Dx_01_2040_Su_Bd_red_WOFlex"    , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_Su_Bd_red_WOFlex_RED.ods"        ),
    ("HR_Dx_01_2040_Su_Sa_red_WithFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_Su_Sa_red_WithFlex_RED.ods"      ),
    ("HR_Dx_01_2040_Su_Sa_red_WOFlex"    , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_Su_Sa_red_WOFlex_RED.ods"        ),
    ("HR_Dx_01_2040_Su_Su_red_WithFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_Su_Su_red_WithFlex_RED.ods"      ),
    ("HR_Dx_01_2040_Su_Su_red_WOFlex"    , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_Su_Su_red_WOFlex_RED.ods"        ),
    ("HR_Dx_01_2040_S_Bd_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_S_Bd_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2040_S_Bd_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_S_Bd_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2040_S_Sa_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_S_Sa_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2040_S_Sa_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_S_Sa_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2040_S_Su_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_S_Su_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2040_S_Su_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_S_Su_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2040_W_Bd_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_W_Bd_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2040_W_Bd_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_W_Bd_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2040_W_Sa_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_W_Sa_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2040_W_Sa_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_W_Sa_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2040_W_Su_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_W_Su_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2040_W_Su_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2040_W_Su_red_WOFlex_RED.ods"         ),

    ("HR_Dx_01_2040_A_Bd_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_A_Bd_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2040_A_Bd_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_A_Bd_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2040_A_Sa_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_A_Sa_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2040_A_Sa_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_A_Sa_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2040_A_Su_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_A_Su_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2040_A_Su_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_A_Su_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2040_Su_Bd_green_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_Su_Bd_green_WithFlex_GREEN.ods"),
    ("HR_Dx_01_2040_Su_Bd_green_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_Su_Bd_green_WOFlex_GREEN.ods"  ),
    ("HR_Dx_01_2040_Su_Sa_green_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_Su_Sa_green_WithFlex_GREEN.ods"),
    ("HR_Dx_01_2040_Su_Sa_green_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_Su_Sa_green_WOFlex_GREEN.ods"  ),
    ("HR_Dx_01_2040_Su_Su_green_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_Su_Su_green_WithFlex_GREEN.ods"),
    ("HR_Dx_01_2040_Su_Su_green_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_Su_Su_green_WOFlex_GREEN.ods"  ),
    ("HR_Dx_01_2040_S_Bd_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_S_Bd_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2040_S_Bd_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_S_Bd_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2040_S_Sa_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_S_Sa_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2040_S_Sa_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_S_Sa_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2040_S_Su_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_S_Su_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2040_S_Su_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_S_Su_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2040_W_Bd_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_W_Bd_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2040_W_Bd_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_W_Bd_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2040_W_Sa_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_W_Sa_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2040_W_Sa_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_W_Sa_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2040_W_Su_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_W_Su_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2040_W_Su_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2040_W_Su_green_WOFlex_GREEN.ods"   ),

    ("HR_Dx_01_2050_A_Bd_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_A_Bd_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2050_A_Bd_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_A_Bd_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2050_A_Sa_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_A_Sa_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2050_A_Sa_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_A_Sa_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2050_A_Su_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_A_Su_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2050_A_Su_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_A_Su_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2050_Su_Bd_brown_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_Su_Bd_brown_WithFlex_BROWN.ods"),
    ("HR_Dx_01_2050_Su_Bd_brown_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_Su_Bd_brown_WOFlex_BROWN.ods"  ),
    ("HR_Dx_01_2050_Su_Sa_brown_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_Su_Sa_brown_WithFlex_BROWN.ods"),
    ("HR_Dx_01_2050_Su_Sa_brown_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_Su_Sa_brown_WOFlex_BROWN.ods"  ),
    ("HR_Dx_01_2050_Su_Su_brown_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_Su_Su_brown_WithFlex_BROWN.ods"),
    ("HR_Dx_01_2050_Su_Su_brown_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_Su_Su_brown_WOFlex_BROWN.ods"  ),
    ("HR_Dx_01_2050_S_Bd_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_S_Bd_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2050_S_Bd_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_S_Bd_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2050_S_Sa_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_S_Sa_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2050_S_Sa_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_S_Sa_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2050_S_Su_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_S_Su_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2050_S_Su_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_S_Su_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2050_W_Bd_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_W_Bd_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2050_W_Bd_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_W_Bd_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2050_W_Sa_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_W_Sa_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2050_W_Sa_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_W_Sa_brown_WOFlex_BROWN.ods"   ),
    ("HR_Dx_01_2050_W_Su_brown_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_W_Su_brown_WithFlex_BROWN.ods" ),
    ("HR_Dx_01_2050_W_Su_brown_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/BROWN/HR_Dx_01_2050_W_Su_brown_WOFlex_BROWN.ods"   ),

    ("HR_Dx_01_2050_A_Bd_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_A_Bd_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2050_A_Bd_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_A_Bd_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2050_A_Sa_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_A_Sa_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2050_A_Sa_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_A_Sa_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2050_A_Su_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_A_Su_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2050_A_Su_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_A_Su_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2050_Su_Bd_red_WithFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_Su_Bd_red_WithFlex_RED.ods"      ),
    ("HR_Dx_01_2050_Su_Bd_red_WOFlex"    , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_Su_Bd_red_WOFlex_RED.ods"        ),
    ("HR_Dx_01_2050_Su_Sa_red_WithFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_Su_Sa_red_WithFlex_RED.ods"      ),
    ("HR_Dx_01_2050_Su_Sa_red_WOFlex"    , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_Su_Sa_red_WOFlex_RED.ods"        ),
    ("HR_Dx_01_2050_Su_Su_red_WithFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_Su_Su_red_WithFlex_RED.ods"      ),
    ("HR_Dx_01_2050_Su_Su_red_WOFlex"    , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_Su_Su_red_WOFlex_RED.ods"        ),
    ("HR_Dx_01_2050_S_Bd_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_S_Bd_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2050_S_Bd_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_S_Bd_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2050_S_Sa_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_S_Sa_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2050_S_Sa_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_S_Sa_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2050_S_Su_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_S_Su_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2050_S_Su_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_S_Su_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2050_W_Bd_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_W_Bd_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2050_W_Bd_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_W_Bd_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2050_W_Sa_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_W_Sa_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2050_W_Sa_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_W_Sa_red_WOFlex_RED.ods"         ),
    ("HR_Dx_01_2050_W_Su_red_WithFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_W_Su_red_WithFlex_RED.ods"       ),
    ("HR_Dx_01_2050_W_Su_red_WOFlex"     , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/RED/HR_Dx_01_2050_W_Su_red_WOFlex_RED.ods"         ),

    ("HR_Dx_01_2050_A_Bd_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_A_Bd_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2050_A_Bd_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_A_Bd_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2050_A_Sa_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_A_Sa_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2050_A_Sa_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_A_Sa_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2050_A_Su_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_A_Su_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2050_A_Su_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_A_Su_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2050_Su_Bd_green_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_Su_Bd_green_WithFlex_GREEN.ods"),
    ("HR_Dx_01_2050_Su_Bd_green_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_Su_Bd_green_WOFlex_GREEN.ods"  ),
    ("HR_Dx_01_2050_Su_Sa_green_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_Su_Sa_green_WithFlex_GREEN.ods"),
    ("HR_Dx_01_2050_Su_Sa_green_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_Su_Sa_green_WOFlex_GREEN.ods"  ),
    ("HR_Dx_01_2050_Su_Su_green_WithFlex", "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_Su_Su_green_WithFlex_GREEN.ods"),
    ("HR_Dx_01_2050_Su_Su_green_WOFlex"  , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_Su_Su_green_WOFlex_GREEN.ods"  ),
    ("HR_Dx_01_2050_S_Bd_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_S_Bd_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2050_S_Bd_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_S_Bd_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2050_S_Sa_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_S_Sa_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2050_S_Sa_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_S_Sa_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2050_S_Su_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_S_Su_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2050_S_Su_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_S_Su_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2050_W_Bd_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_W_Bd_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2050_W_Bd_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_W_Bd_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2050_W_Sa_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_W_Sa_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2050_W_Sa_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_W_Sa_green_WOFlex_GREEN.ods"   ),
    ("HR_Dx_01_2050_W_Su_green_WithFlex" , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_W_Su_green_WithFlex_GREEN.ods" ),
    ("HR_Dx_01_2050_W_Su_green_WOFlex"   , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/HR DX 01/GREEN/HR_Dx_01_2050_W_Su_green_WOFlex_GREEN.ods"   ),

    #=
    ### multiline break
    ("ES_Dx_03_2030_S_Bd"       , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/ES DX 03/ES_Dx_03_2030_S_Bd.ods"),
    ("ES_Dx_03_2030_W_Bd"       , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/ES DX 03/ES_Dx_03_2030_W_Bd.ods"),
    ("ES_Dx_03_2040_S_Bd"       , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/ES DX 03/ES_Dx_03_2040_S_Bd.ods"),
    ("ES_Dx_03_2040_W_Bd"       , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/ES DX 03/ES_Dx_03_2040_W_Bd.ods"),
    ("ES_Dx_03_2050_S_Bd"       , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/ES DX 03/ES_Dx_03_2050_S_Bd.ods"),
    ("ES_Dx_03_2050_W_Bd"       , "input_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/ES DX 03/ES_Dx_03_2050_W_Bd.ods"),
    =#
    ])

CaseNames   = collect(keys(CasesList))
N_Cases     = length(CasesList) # sizeof(CasesList) here returns something else that i dont understand
OutLog      = "output_data/March2023_LoadDataFromAttestCloud/EditedInputFiles/OutLog.xlsx"
LogThis(;CaseLog, Append=false) = append_xl_row(workbook_path=OutLog,
                                                sheet_name="Log",
                                                row_data=collect(CaseLog),
                                                Append=Append)

LogLabels   = ["#","CaseName","InputFile","OutputFile","Start_Time", "Init_PF","Viols_Exist","Opt_Commenced","Opt_Found"]
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

TheKey    = "HR_Dx_01_2050_A_Bd_green_WithFlex"

# ==============================================================================
# for TheKey in CaseNames # Added by Baraa
#### i = findall(x -> occursin(TheKey, x),CaseNames) # also works
i               = findfirst(isequal(TheKey),CaseNames)

# for i in range(1, N_Cases) # Added by Baraa
TheKey          = CaseNames[i]

filename_mat    = CasesList[TheKey]
filename_addt   = replace(filename_mat  , ".ods"      =>"_flex.ods")
file_op         = replace(filename_mat  , "input_data"=>"output_data")
file_op         = replace(file_op       , ".ods"      =>"_output.xlsx")

@printf("Case #%3i/%i: %s\n", i, N_Cases, TheKey)
StartTime = now()
CaseLog = Array([i, TheKey, filename_addt, file_op, string(StartTime), false, false, false, false])
# open(replace(file_op, ".ods" => ".txt"), "w") do out
#     redirect_stdout(out) do

println("Running the case: ",filename_mat, " @", StartTime)
LogThis(CaseLog=CaseLog, Append=true)
# break

################################################################################
### -------------------- Reading Input Data Files -------------------------- ###
################################################################################
include("interface_excel.jl")  # Reading Input Network data and additional flexibility data corresponding to the network
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
(nPprf_wind_sc, nSc, rheader_pProfile_wind_sc, rdata_pProfile_wind_sc,
nPprf_pv_sc, rheader_pProfile_pv_sc, rdata_pProfile_pv_sc,
n_prob, rheader_prob_sc, rdata_prob_sc) = interface_scenarios_fun(filename)

################################################################################
### ---------------------- Constants Declaration -------------------------- ####
################################################################################
# "INCLUDE" actually runs/executes the file. not just makes it available to be called.
include("constants.jl") # this file contains all the parameters selected by the programmer/user
include("constants_pf_opf.jl")

if  occursin("WOFlex", filename_mat)
    println("Case withOUT Flexibility. Skipping")
    # continue
    flex_fl   = 0
    fl_bin    = 0
    flex_str  = 0
    str_bin   = 0
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
 cost_c_str) = data_storage(rheader_storage,
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

p_curt_pf = zeros(nSc_pf, nTP_pf, nCurt_gen) ./ sbase
p_dis_pf = zeros(nSc_pf, nTP_pf, nStr_active) ./ sbase
p_ch_pf = zeros(nSc_pf, nTP_pf, nStr_active) ./ sbase
p_od_pf = zeros(nSc_pf, nTP_pf, nFl) ./ sbase
p_ud_pf = zeros(nSc_pf, nTP_pf, nFl) ./ sbase
q_dg_pf = tan(acos(dg_ac_pf)) .* (pg_max[:, :, idx_curt_dg] - p_curt_pf)
oltc_tap_init = tratio_init                         # Unless the profile of tap ratio is given, it is fixed to the same value in all time-periods and all scenarios
p_strg_pf = p_dis_pf - p_ch_pf

pg_max_pf = pg_max[pf_scenario, :, :]
pg_min_pf = pg_min[pf_scenario, :, :]
qg_max_pf = qg_max[pf_scenario, :, :]
qg_min_pf = pg_min[pf_scenario, :, :]

pg_max_pf = reshape(pg_max_pf, (nSc_pf_new, nTP_pf_new, nGens))
pg_min_pf = reshape(pg_min_pf, (nSc_pf_new, nTP_pf_new, nGens))
qg_max_pf = reshape(qg_max_pf, (nSc_pf_new, nTP_pf_new, nGens))
qg_min_pf = reshape(qg_min_pf, (nSc_pf_new, nTP_pf_new, nGens))
prob_scs_pf = prob_scs[pf_scenario, 1]


q_dg_pf = tan(acos(dg_ac_pf)) .* (pg_max[:, :, idx_curt_dg] - p_curt_pf)
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

display("Initial AC-Power Flow before MILP Model")
(vol_nodes_mag_pf,
vol_nodes_theta_pf,
vol_rect_pf) = ac_power_flow_opf_int(nBus,
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

        XLSX.writetable!(xf[SheetName],
                         Data,
                         ;anchor_cell=XLSX.CellRef("A1"))
    # end
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
    if  sizeof(vol_viol_int)>0
        WriteThis  = [map(string,vol_viol_int[i]) for i in 1:size(vol_viol_int,1)]
        WriteThis2 = hcat(collect.(WriteThis)...)
        WriteThis3 = DataFrame(WriteThis2, :auto)
        Append2Excel(xf, WriteThis3, "Vlt_Viol")
    end

    if  sizeof(crnt_viol_int)>0
        WriteThis  = [map(string,crnt_viol_int[i]) for i in 1:size(crnt_viol_int,1)]
        WriteThis2 = hcat(collect.(WriteThis)...)
        WriteThis3 = DataFrame(WriteThis2, :auto)
        Append2Excel(xf, WriteThis3, "Crnt_Viol")
    end
end

println("Checkout initial PF results @ ",file_op)
CaseLog[6] = 1                  # 6:"Init_PF", 7:"Viols_Exist", 8:"Opt_Commenced", 9:"Opt_Found"
CaseLog[7] = sum(N_violations)  # 6:"Init_PF", 7:"Viols_Exist", 8:"Opt_Commenced", 9:"Opt_Found"

if sum(N_violations)==0
    alert("JULIA: No violations. This case is a waste of time")
    LogThis(CaseLog=CaseLog, Append=false)
    # continue
    throw(error("No violations. This case is a waste of time"))

else
    println(string(nVol_viol_int),  " Voltage limit violations.\n",
    string(nCrnt_viol_int), " Current limit violations.\n",
    string(sum(N_violations)), " Total Violations were found.\n",
    )
    CaseLog[8] = 1              # 6:"Init_PF", 7:"Viols_Exist", 8:"Opt_Commenced", 9:"Opt_Found"
    LogThis(CaseLog=CaseLog, Append=false)
    alert("JULIA: Violations FOUND. Optimization will commence")
    println("Flexibility optimization will commence\n",repeat("=",60),repeat("\n",30))
end

# throw(error("Breaking for no reason"))
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











































try
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


(vol_viol_nw, # voltage violations of the network
vol_viol_nw_max, #
crnt_viol_nw, # current (amperes) violations
crnt_viol_nw_max,
vol_nodes_pf_cnv, # voltage magnitdues at the node after the PF problem converges (post optimization check/test)
vol_nodes_ang_pf_cnv, # voltages at nodes, angles.
br_crnt_si_pf_cnv, # branch current / flow
sql_obj, # objective function value / cost funciton
sql_time, # total time of solution
rel_max_error, # relative max error
rel_avg_error, # relative average error. explained in the paper [inner loop \delta S]
norm_s_err_max_power, # normal apparent power error, max
norm_s_err_avg, # normal apparent power error, average
p_curt_gen, # curtailment of active power, per generator. this is an array. (gen/time/scenario)
q_curt_gen, # curtailment of reactive power, indirect as a result of P curtailment.
total_p_curt, # this is total curtailment. p_curt_gen is per generator. this is a scalar
str_dis, # storage discharge, per (ESS/time/scenario), quantity is reported in KW
str_ch,  # storage    charge, per (ESS/time/scenario), quantity is reported in KW.
viol_nodes_new, # i dont know
p_curt_lin, # p_curtail linear
min_vol_limit, # minimum violation
max_vol_limit,
oltc_ratio, # on load tap changer ratio, per (transformer/time period/scenario)
p_err_max_power, # not sure
q_err_max_power,
vm_ac, # voltage magnitude in AC. (AC means polar form, |V|/_angle)
va_ac, # voltage angle in AC
vol_viol_cn_milp,
model_infs,
br_crnt_opf_app,
br_crnt_opf_ex,
p_res_curt_op,
q_res_op,
fl_inc_op, # flexible load increment per (load unit / time period / scenario)
fl_dec_op, # flexible load decrement
p_ch_op,   # P charging optimum. check this against str_ch
p_dis_op,  # P discharge optimum. check this against str_dis
p_strg_op, # P_net = p_charge - p_disch
p_slack_pf_ip, # slack power coming from PF solution based on IPOPT solver
q_slack_pf_ip, # slack reactive power
br_pwr_si_pf_cnv, # branch power, in SI units (MVA) (not in PU)
v_rect_euler, # doesn't make sense. ignore
v_rect_cmplx) = sql_loop(vol_viol_max,
                        vol_viol_tol,
                        crnt_viol_max,
                        crnt_viol_tol,
                        lin_itr,
                        lin_itr_max,
                        vol_nodes_mag_pf,
                        vol_nodes_theta_pf,
                        nSc,
                        nTP,
                        nBus,
                        nNcurt_gen,
                        nCurt_gen,
                        nTrsf,
                        nStr_active,
                        nFl,
                        nd_fl,
                        flex_apc,
                        flex_oltc,
                        flex_adpf,
                        flex_str,
                        flex_fl,
                        str_bin,
                        fl_bin,
                        rdata_buses,
                        nLines,
                        trsf_fict_nodes,
                        tap_ratio_range,
                        oltc_bin,
                        prob_scs,
                        time_step,
                        tdata,
                        bus_data_Ssheet,
                        bus_data_lsheet,
                        cost_a_str,
                        cost_b_str,
                        cost_c_str,
                        cost_load_inc,
                        cost_load_dec,
                        nw_buses,
                        rdata_loads,
                        node_data,
                        nd_curt_gen,
                        nd_ncurt_gen,
                        p_load,
                        q_load,
                        idx_Gs_lsheet,
                        idx_Bs_lsheet,
                        pg_max,
                        yii_sh,
                        i_curt_gens,
                        sbase,
                        dg_pf,
                        iStr_active,
                        yij_line,
                        dLines,
                        i_ncurt_gens,
                        nw_lines,
                        bus_data_gsheet,
                        pg_min,
                        qg_min,
                        qg_max,
                        pgen_tol,
                        qgen_tol,
                        p_tol,
                        q_tol,
                        rcvr_branch,
                        rcvr_inj,
                        rcvr_tso_dso,
                        Ibase,
                        idx_slack_bus_pf,
                        solver_ipopt,
                        solver_cbc,
                        solver_bonmin,
                        sql_itr,
                        sql_itr_max,
                        s_inj_error_rel,
                        s_error_tol,
                        angle_lb,
                        angle_ub,
                        cycleList,
                        idx_slack_from,
                        idx_slack_to,
                        num_slack_cnctd_nodes,
                        slack_cnctd_nodes,
                        slack_nd,
                        stoch_model,
                        tratio_init,
                        rdata_storage,
                        load_theta,
                        nd_Str_active,
                        slack_bus_type,
                        vol_ctrl_bus_type,
                        load_bus_type,
                        rdata_gens,
                        nw_pPrf_data_load,
                        nw_qPrf_data_load,
                        nw_pPrf_data_gen_max,
                        nw_qPrf_data_gen_max,
                        scenario_data_p_min,
                        scenario_data_p_max,
                        scenario_data_q_min,
                        scenario_data_q_max,
                        nw_buses_pf,
                        nw_lines_pf,
                        nw_loads_pf,
                        nw_gens_pf,
                        nw_gcost_pf,
                        nw_sbase_pf,
                        v_initial_pf,
                        v_mag_pf,
                        v_angle_pf,
                        max_mismatch_pf,
                        epsilon,
                        iteration,
                        itr_max,
                        ordata_buses_pf,
                        ybase,
                        vbase,
                        I_rat,
                        vol_viol_nw,
                        vol_viol_nw_max,
                        crnt_viol_nw,
                        crnt_viol_nw_max,
                        min_vol_limit,
                        max_vol_limit,
                        max_crnt_limit,
                        term_status,
                        vol_cstr_tol,
                        br_crnt_si_pf_cnv,
                        solver_cplex,
                        int_tol,
                        opt_tol,
                        num_thread,
                        rdata_trsf,
                        nTrsf_s1,
                        nw_trsf,
                        n_tap,
                        br_pwr_si_pf_cnv)

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
avg_p_gen_curt = reshape(sum(prob_scs .* p_res_curt_op, dims=1), (nTP, nCurt_gen))
avg_q_gen_curt = reshape(sum(prob_scs .* q_res_op, dims=1), (nTP, nCurt_gen))

avg_p_gen_curt = permutedims(avg_p_gen_curt)
avg_q_gen_curt = permutedims(avg_q_gen_curt)

avg_p_gen_curt = hcat(nodes_curt_gen, avg_p_gen_curt)
avg_q_gen_curt = hcat(nodes_curt_gen, avg_q_gen_curt)

# ################################################################################
# ##### ----------------- Flexible Load Data --------------------------------#####
# ################################################################################
nodes_fl = Int64.(rdata_loads[iFl, 1])
avg_fl_inc = prob_scs .* fl_inc_op
avg_fl_dec = prob_scs .* fl_dec_op

avg_fl_inc_op = permutedims(avg_fl_inc[sc_max_prob, :, :])
avg_fl_dec_op = permutedims(avg_fl_dec[sc_max_prob, :, :])

avg_fl_inc_op = hcat(nodes_fl, avg_fl_inc_op)
avg_fl_dec_op = hcat(nodes_fl, avg_fl_dec_op)

# ################################################################################
# ##### ----------------- Energy Storage Data --------------------------------####
# ################################################################################
nodes_str = Int64.(rdata_storage[iStr_active, 1])
avg_p_ch_nd = prob_scs .* p_ch_op
avg_p_dch_nd = prob_scs .* p_dis_op

avg_p_ch_op = permutedims(avg_p_ch_nd[sc_max_prob, :, :])
avg_p_dch_op = permutedims(avg_p_dch_nd[sc_max_prob, :, :])

avg_p_ch_op = hcat(nodes_str, avg_p_ch_op)
avg_p_dch_op = hcat(nodes_str, avg_p_dch_op)
# ###############################################################################
#
header_op = [string("t", i) for i in 0:24]
header_op[1] = "nodes"
# header_op = ["nodes", "t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10", "t11",
#     "t12", "t13", "t14", "t15", "t16", "t17", "t18", "t19", "t20", "t21", "t22", "t23", "t24"]

df_p_curt = DataFrame(avg_p_gen_curt, header_op)
df_q_res = DataFrame(avg_q_gen_curt, header_op)

df_fl_inc = DataFrame(avg_fl_inc_op, header_op)
df_fl_dec = DataFrame(avg_fl_dec_op, header_op)

df_str_ch = DataFrame(avg_p_ch_op, header_op)
df_str_dch = DataFrame(avg_p_dch_op, header_op)
#
# ################################################################################
# ################## Cost of Procurement of Flexibility ##########################
# ################################################################################
header_cost = ["cost_flex"]
data_obj = reshape(sql_obj, (sql_itr_max * lin_itr_max))
idx_0 = findall(x -> x != 0, data_obj)
cost_flex_proc = [data_obj[idx_0[end]]]
df_cost_pr = DataFrame(cost=cost_flex_proc)
#
# ################################################################################
# ####------------- Writing Output Data to Excel File ---------------------- #####
# ################################################################################

function AppendStats(DF)
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
end

println("Finished running && Writing to Excel @", now())
XLSX.openxlsx(file_op, mode="rw") do xf
    Append2Excel(xf, AppendStats(df_p_curt), "APC_MW")
    Append2Excel(xf, AppendStats(df_str_ch), "EES_CH_MW")
    Append2Excel(xf, AppendStats(df_str_dch), "EES_DCH_MW")
    Append2Excel(xf, AppendStats(df_fl_inc), "FL_OD_MW")
    Append2Excel(xf, AppendStats(df_fl_dec), "FL_UD_MW")
    Append2Excel(xf, df_cost_pr, "COST")
end
println("Finished writing final results to excel.\nFind final results @ ",file_op)
CaseLog[9] = 1  # 6:"Init_PF", 7:"Viols_Exist", 8:"Opt_Commenced", 9:"Opt_Found"
LogThis(CaseLog=CaseLog, Append=false)
alert(string("Finished optimization. Results exported to Excel", TheKey))
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

catch Err
    println(Err)
    CaseLog[9] = -1  # 6:"Init_PF", 7:"Viols_Exist", 8:"Opt_Commenced", 9:"Opt_Found"
    @printf("---===<<<\tError in:#%3i/%i \"%s\"\t>>>===---\n", i, N_Cases, TheKey)

end # END OF TRY
# end # END OF: "for (key, value) in CasesList"
