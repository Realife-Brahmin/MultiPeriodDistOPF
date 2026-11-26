"""
OpenDSS Validator for tADMM
Simple powerflow validation utilities
"""

using OpenDSSDirect
using Printf

"""
    print_powerflow_summary(data::Dict)

Print a summary of the OpenDSS powerflow results.
"""
function print_powerflow_summary(data::Dict)
    # Run powerflow
    OpenDSSDirect.Text.Command("Solve")
    
    # Check if converged
    converged = OpenDSSDirect.Solution.Converged()
    
    println("\n", "="^40)
    println("Initial Powerflow")
    println("="^40)
    
    if converged
        println("Powerflow completed successfully!")
    else
        println("⚠ Powerflow did not converge!")
    end
    
    # Get total powers
    total_power_kW = 0.0
    total_power_kVAR = 0.0
    
    # Sum up all load powers
    load_names = OpenDSSDirect.Loads.AllNames()
    for load_name in load_names
        OpenDSSDirect.Loads.Name(load_name)
        total_power_kW += OpenDSSDirect.CktElement.Powers()[1]  # Real power
        total_power_kVAR += OpenDSSDirect.CktElement.Powers()[2]  # Reactive power
    end
    
    # Get losses
    losses = OpenDSSDirect.Circuit.Losses()
    loss_kW = real(losses[1]) / 1000  # Convert W to kW, take real part
    
    # Get PV power
    pv_power_kW = 0.0
    gen_names = OpenDSSDirect.Generators.AllNames()
    if !isempty(gen_names) && gen_names[1] != "NONE"
        for gen_name in gen_names
            OpenDSSDirect.Generators.Name(gen_name)
            pv_power_kW += abs(OpenDSSDirect.CktElement.Powers()[1])
        end
    end
    
    # Get battery power
    batt_power_kW = 0.0
    storage_names = OpenDSSDirect.Storages.AllNames()
    if !isempty(storage_names) && storage_names[1] != "NONE"
        for storage_name in storage_names
            OpenDSSDirect.Storages.Name(storage_name)
            batt_power_kW += OpenDSSDirect.CktElement.Powers()[1]
        end
    end
    
    # Calculate substation power: load + losses - generation - battery
    subs_power_kW = abs(total_power_kW) + loss_kW - pv_power_kW - batt_power_kW
    
    voltages = OpenDSSDirect.Circuit.AllBusMagPu()

    # Print summary
    println()
    @printf "[OpenDSS] Total system load: %.2f kW\n" abs(total_power_kW)
    @printf "[OpenDSS] Substation power: %.2f kW\n" subs_power_kW
    @printf "[OpenDSS] Total PV dispatch: %.2f kW\n" pv_power_kW
    @printf "[OpenDSS] Total battery dispatch: %.2f kW\n" batt_power_kW
    @printf "[OpenDSS] Total system losses: %.2f kW\n" loss_kW
    # @printf "[OpenDSS] Voltage magnitudes: \n" voltages
    # for (bus, v) in zip(OpenDSSDirect.Circuit.AllBusNames(), voltages)
    #     @printf "  - %s: %.4f pu\n" bus v
    # end
    println()
end
