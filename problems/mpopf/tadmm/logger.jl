#=
Logger utilities for tADMM
Functions for writing model summaries, data validation, and diagnostics
=#

using JuMP

"""
    write_model_summary(model, data, output_file)

Write a comprehensive model summary to a text file including:
- Organized constraint listing
- Variable and constraint counts
- Data validation checks
- Power balance verification
- Battery parameters (if present)

# Arguments
- `model`: JuMP model
- `data`: Dictionary containing system data with keys:
    - :NLset, :Dset, :Bset, :Nset, :Nm1set, :L1set, :Lm1set (node/branch sets)
    - :p_L_pu, :p_D_pu (load and PV profiles in p.u.)
    - :B0_pu, :B_R_pu, :P_B_R_pu, :soc_min, :soc_max (battery parameters)
    - :P_BASE (base power in kW)
- `output_file`: Path to output file
"""
function write_model_summary(model, data, output_file)
    # Extract data fields directly
    NLset = data[:NLset]
    Dset = data[:Dset]
    Bset = data[:Bset]
    Nset = data[:Nset]
    Nm1set = data[:Nm1set]
    L1set = data[:L1set]
    Lm1set = data[:Lm1set]
    p_L_pu = data[:p_L_pu]
    p_D_pu = data[:p_D_pu]
    P_BASE = data[:kVA_B]  # Base power in kW
    B0_pu = data[:B0_pu]
    B_R_pu = data[:B_R_pu]
    P_B_R_pu = data[:P_B_R_pu]
    soc_min = data[:soc_min]
    soc_max = data[:soc_max]
    
    open(output_file, "w") do io
        println(io, "="^80)
        println(io, "MODEL SUMMARY")
        println(io, "="^80)
        
        # Manually print all constraints in organized order
        println(io, "\n" * "="^80)
        println(io, "ALL CONSTRAINTS (ORGANIZED BY TYPE AND TIME)")
        println(io, "="^80)
        
        # Define constraint name patterns in desired order
        constraint_groups = [
            ("SUBSTATION POWER BALANCE", r"RealPowerBalance_Substation_t\d+"),
            ("NODAL REAL POWER BALANCE", r"RealPowerBalance_Node\d+_t\d+"),
            ("NODAL REACTIVE POWER BALANCE (Substation)", r"ReactivePowerBalance_Substation_t\d+"),
            ("NODAL REACTIVE POWER BALANCE", r"ReactivePowerBalance_Node\d+_t\d+"),
            ("KVL CONSTRAINTS", r"KVL_Branch_\d+_\d+_t\d+"),
            ("FIXED SUBSTATION VOLTAGE", r"FixedSubstationVoltage_t\d+"),
            ("BATTERY SOC TRAJECTORY", r"BatterySOC.*_t\d+"),
            ("VOLTAGE LIMITS", r"VoltageLimits_Node\d+_t\d+"),
            ("PV REACTIVE POWER LIMITS", r"PVReactiveLimits_DER\d+_t\d+"),
            ("BATTERY SOC LIMITS", r"BatterySOCLimits_\d+_t\d+"),
            ("BATTERY POWER LIMITS", r"BatteryPowerLimits_\d+_t\d+"),
        ]
        
        # Collect all constraint references
        all_constraint_refs = []
        for (F, S) in list_of_constraint_types(model)
            for con in all_constraints(model, F, S)
                push!(all_constraint_refs, (name(con), con))
            end
        end
        
        # Print constraints by group
        for (group_name, pattern) in constraint_groups
            println(io, "\n" * "-"^80)
            println(io, group_name)
            println(io, "-"^80)
            
            # Filter and sort constraints matching this pattern
            matching_cons = filter(x -> occursin(pattern, x[1]), all_constraint_refs)
            sort!(matching_cons, by=x -> x[1])  # Alphabetical sort within group
            
            for (con_name, con) in matching_cons
                println(io, con_name, ": ", con)
            end
            
            println(io, "Total: $(length(matching_cons)) constraints")
        end
        
        # Print any remaining constraints not captured above
        printed_names = Set{String}()
        for (group_name, pattern) in constraint_groups
            for (con_name, con) in all_constraint_refs
                if occursin(pattern, con_name)
                    push!(printed_names, con_name)
                end
            end
        end
        
        remaining = filter(x -> !(x[1] in printed_names), all_constraint_refs)
        if !isempty(remaining)
            println(io, "\n" * "-"^80)
            println(io, "OTHER CONSTRAINTS")
            println(io, "-"^80)
            for (con_name, con) in remaining
                println(io, con_name, ": ", con)
            end
        end
        
        println(io, "\n--- VARIABLE BOUNDS CHECK ---")
        println(io, "Number of variables: ", num_variables(model))
        println(io, "Number of constraints: ", num_constraints(model; count_variable_in_set_constraints=true))
        
        # Check for potential issues
        println(io, "\n--- DATA VALIDATION ---")
        println(io, "Load nodes (NLset): ", NLset)
        println(io, "PV nodes (Dset): ", Dset)
        println(io, "Battery nodes (Bset): ", Bset)
        println(io, "All nodes (Nset): ", Nset)
        println(io, "Non-substation nodes (Nm1set): ", Nm1set)
        println(io, "Substation branches (L1set): ", L1set)
        println(io, "Other branches (Lm1set): ", Lm1set)
        
        # Check power balance feasibility
        println(io, "\n--- POWER BALANCE CHECK (t=1) ---")
        total_load_t1 = sum(p_L_pu[j, 1] for j in NLset)
        total_pv_t1 = sum(p_D_pu[j, 1] for j in Dset)
        println(io, "Total load at t=1: $(total_load_t1 * P_BASE) kW")
        println(io, "Total PV at t=1: $(total_pv_t1 * P_BASE) kW")
        println(io, "Net demand: $((total_load_t1 - total_pv_t1) * P_BASE) kW")
        
        # Check battery parameters
        if !isempty(Bset)
            println(io, "\n--- BATTERY PARAMETERS ---")
            for b in Bset
                println(io, "Battery $b:")
                println(io, "  Initial SOC: $(B0_pu[b] * P_BASE) kWh")
                println(io, "  Capacity: $(B_R_pu[b] * P_BASE) kWh")
                println(io, "  Power rating: $(P_B_R_pu[b] * P_BASE) kW")
                println(io, "  SOC min: $(soc_min[b] * 100)%")
                println(io, "  SOC max: $(soc_max[b] * 100)%")
                println(io, "  Min energy: $(soc_min[b] * B_R_pu[b] * P_BASE) kWh")
                println(io, "  Max energy: $(soc_max[b] * B_R_pu[b] * P_BASE) kWh")
            end
        end
    end
    
    println("Model summary saved to: $output_file")
end
