# test_dss_parser.jl
# Test script for new DSS-based parsing functions

include("./src/setupMultiPeriodDistOPF.jl")


systemName = "ads10A_1ph"
T = 24

println("="^80)
println("Testing new DSS-based parser for system: $systemName")
println("="^80)

# Test 1: Load system into DSS
println("\n[Test 1] Loading system into OpenDSS...")
dss_file = Parser.load_system_in_dss(systemName)
println("✓ Successfully loaded: $dss_file")

# Test 2: Get full system configuration
println("\n[Test 2] Extracting system configuration from DSS...")
data = Parser.get_system_config_from_dss(systemName, T,
    linearizedModel=false,
    temporal_decmp=false,
    solver="Ipopt",
    tSOC_hard=false,
    relax_terminal_soc_constraint=true)

println("✓ Successfully extracted data")

# Test 3: Verify key data structures
println("\n[Test 3] Verifying extracted data...")
@unpack Nset, Lset, NLset, Dset, Bset, T, kVA_B, kV_B = data

println("  System: $systemName")
println("  Time steps: $T")
println("  Base values: kVA_B = $kVA_B, kV_B = $kV_B")
println("  Number of buses: $(length(Nset))")
println("  Number of branches: $(length(Lset))")
println("  Number of loads: $(length(NLset))")
println("  Number of PVs: $(length(Dset))")
println("  Number of batteries: $(length(Bset))")

# Test 4: Verify data compatibility with existing model builder
println("\n[Test 4] Testing compatibility with model builder...")
try
    @unpack kVA_B_dict, kV_B_dict, MVA_B_dict, rdict, xdict, rdict_pu, xdict_pu = data
    println("  ✓ Base value dictionaries present")
    
    @unpack p_L_R, q_L_R, p_L, q_L = data
    println("  ✓ Load data present")
    
    @unpack p_D_R, S_D_R, p_D = data
    println("  ✓ PV data present")
    
    @unpack B0, B_R, P_B_R, eta_C, eta_D = data
    println("  ✓ Battery data present")
    
    @unpack Vminpu, Vmaxpu = data
    println("  ✓ Voltage limit data present")
    
    println("\n✓ All required data structures are present and compatible!")
    
catch e
    println("✗ Error verifying data structures: $e")
end

println("\n" * "="^80)
println("Testing complete!")
println("="^80)

# Optional: Try building a model with the extracted data
println("\n[Optional] Attempting to build MPOPF model...")
try
    if data[:linearizedModel]
        modelDict = MB.build_MPOPF_1ph_L_model_t_in_Tset(data)
        println("✓ Successfully built linearized MPOPF model!")
    else
        modelDict = MB.build_MPOPF_1ph_NL_model_t_in_Tset(data)
        println("✓ Successfully built nonlinear MPOPF model!")
    end
    
    @unpack model = modelDict
    println("  Model has $(num_variables(model)) variables")
    println("  Model has $(num_constraints(model)) constraints")
catch e
    println("✗ Could not build model: $e")
    println("  (This might be expected if there are missing data fields)")
end
