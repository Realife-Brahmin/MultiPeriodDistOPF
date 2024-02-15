################################################################################
##---------- Constants to be used in the AC Power Flow Model -------------------
################################################################################
# v_initial_pf    = 1.00
max_mismatch_pf = 1
epsilon         = 1e-6
iteration       = 0
itr_max         = 10
nTP             = 24
scenario        = 1     # Set scenario = 1 for S-MP-OPF problem and scenario = 0 for D-SP or D-MP OPF models
################################################################################
### --------- Parameters for Proposed Sequential Linear Algorithm ---------- ###
################################################################################
sql_itr         = 1  # starting inner-loop iteration of SLA
lin_itr         = 1  # starting outer-loop iteration of SLA
sql_itr_max     = 1  # Maximum number of inner-loop iterations
lin_itr_max     = 1  # maximum number of outer-loop iterations

vol_viol_tol    = 1 # Feasibility violation threshold. The value is given in %
crnt_viol_tol   = 1 # feasibility violation threshold. The value is given in %

################################################################################
### ---------- Parameters for various flexibility options ------------------ ###
################################################################################
load_inc_prct   = 0.50   # % increase in the flexible load
load_dec_prct   = 0.50   # % decrease in the flexible load
pf              = 0.95   # power factor of a load
dg_pf           = 0.90   # Used to set the injected reactive power of RENEWABLE DG'S. Must be set to 1 if no reactive power is injected from these DGs
ratio_init      = 1.00   #
# n_tap         =  09    # (tap ratio to be divided into 10 step)
################################################################################
# Parameters for Setting the Convergence of Optimization based AC Power Flow #
################################################################################
vol_cstr_tol    = 1*1e-3   # cannot be less than 1e-3
p_cstr_fctr     = 1*1e-1
q_cstr_fctr     = 1*1e-1
################################################################################
### ----------------- Selection of Flexibility Options ----------------------###
################################################################################
flex_apc  = 1     # Active Power Curtailment. This option will always be set to 1.

flex_adpf = 1    # Reactive power provision from RES

flex_str  = 1    # Storages
str_bin   = 1    # binary variable modeling the electrical energy storage

flex_fl   = 1    # Flexible Load
fl_bin    = 1    # binary variable modeling the fleixble load

flex_oltc = 1    # On-load tap changing transfomer
oltc_bin  = 1    # Binary variables for the modeling of OLTC transformer

################################################################################
### ---------------------- Selection of Solver ----------------------------- ###
################################################################################
solver_ipopt  = 0
solver_cbc    = 0
solver_bonmin = 0
solver_cplex  = 1

## ------------------------- CPLEX Paramter Values ------------------------- ###
int_tol    = 1e-5          # Parameter for setting integrality tolerance (applicable to all integer problems)
opt_tol    = 1e-5          # Optimality Tolerance
num_thread = 1             # number of cores to be used in parallel

################################################################################
### --------------------- Scenarios Data -------------------------------- ######
################################################################################
if scenario != 0
    (value,idx) = findmax(rdata_prob_sc)
    sc_max_prob = idx[1][1]
elseif scenario == 0
    nSc = 1
    sc_max_prob = 1
end
