"""
Convert IEEE 2552 node system data from MATLAB format to OpenDSS format.
This script reads the linedata.txt and powerdata.txt files from the ieee2552_rahul/NOPV folder
and generates OpenDSS files similar to the ieee123A_1ph system.
"""

using DelimitedFiles

# Define paths - script is in rawData directory
script_dir = @__DIR__
input_dir = joinpath(script_dir, "ieee2552_rahul", "NOPV")
output_dir = joinpath(script_dir, "ieee2552_1ph")

# Create output directory if it doesn't exist
mkpath(output_dir)

println("Reading data from: $input_dir")
println("Writing OpenDSS files to: $output_dir")

# Read powerdata.txt
# Columns: bus_i, type, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV, zone, Vmax, Vmin
powerdata = readdlm(joinpath(input_dir, "powerdata.txt"), '\t', Float64)
println("Read $(size(powerdata, 1)) buses from powerdata.txt")

# Read linedata.txt  
# Columns: fbus, tbus, r, x
linedata = readdlm(joinpath(input_dir, "linedata.txt"), '\t', Float64)
println("Read $(size(linedata, 1)) lines from linedata.txt")

# Extract parameters from the data
num_buses = size(powerdata, 1)
num_lines = size(linedata, 1)

# Base voltage (kV line-to-line) from powerdata
base_kv_ll = powerdata[1, 10]  # Column 10 is baseKV
# For 1-phase, use line-to-neutral voltage
base_kv_ln = base_kv_ll / sqrt(3)

println("\nSystem parameters:")
println("  Number of buses: $num_buses")
println("  Number of lines: $num_lines")
println("  Base voltage (L-L): $base_kv_ll kV")
println("  Base voltage (L-N): $(round(base_kv_ln, digits=4)) kV")

# ===== Generate Master.dss =====
println("\nGenerating Master.dss...")
open(joinpath(output_dir, "Master.dss"), "w") do f
    write(f, """Clear

New Circuit.ieee2552_1ph  

Edit "Vsource.source" bus1=1 pu=1.05 R0=0 X0=0.00000001 phases=1 basekv=$(round(base_kv_ln, digits=4)) R1=0 X1=0.00000001

Redirect BranchData.dss

Redirect ../LoadShapeDefault.dss
Redirect ../LoadShapePVDefault.dss
Redirect Loads.dss     ! Balanced Loads
// Redirect  Capacitor.dss
Redirect PVSystem.dss
Redirect Storage.dss

Set VoltageBases = [$base_kv_ll]
Set mode = Daily 
Set stepsize = 1h
Set number = 1

! Let DSS estimate the voltage bases
CalcVoltageBases     ! This also establishes the bus list

""")
end

# ===== Generate BranchData.dss =====
println("Generating BranchData.dss...")
open(joinpath(output_dir, "BranchData.dss"), "w") do f
    for i in 1:num_lines
        fbus = Int(linedata[i, 1])
        tbus = Int(linedata[i, 2])
        r = linedata[i, 3]
        x = linedata[i, 4]
        
        # Ensure minimum impedance values to avoid numerical issues
        r = max(r, 1e-06)
        x = max(x, 1e-06)
        
        write(f, "New Line.L$i Phases=1 Bus1=$fbus Bus2=$tbus r1=$r x1=$x\n")
    end
end

# ===== Generate Loads.dss =====
println("Generating Loads.dss...")
open(joinpath(output_dir, "Loads.dss"), "w") do f
    write(f, "! Specifying Loads as Model 1 and 2 type loads. For no CVR, Model 2 loads will be blank\n")
    
    load_count = 0
    for i in 1:num_buses
        bus_num = Int(powerdata[i, 1])
        pd_mw = powerdata[i, 3]  # Active power in MW
        qd_mvar = powerdata[i, 4]  # Reactive power in MVAr
        
        # Only create load if there is non-zero power
        if pd_mw > 0 || qd_mvar > 0
            # Convert MW to kW
            pd_kw = pd_mw * 1000
            qd_kvar = qd_mvar * 1000
            
            write(f, "New load.S$(bus_num)P Bus=$bus_num.1 Phases=1 Model=1 ")
            write(f, "kV=$(round(base_kv_ln, digits=4)) ")
            write(f, "kW=$(round(pd_kw, digits=6)) ")
            write(f, "kVAr=$(round(qd_kvar, digits=6)) ")
            write(f, "Vminpu=0.95 Vmaxpu=1.05 Daily=LoadShapeLoadDefault\n")
            
            load_count += 1
        end
    end
    
    println("  Created $load_count loads")
end

# ===== Generate PVSystem.dss (empty for now) =====
println("Generating PVSystem.dss (empty)...")
open(joinpath(output_dir, "PVSystem.dss"), "w") do f
    write(f, "! PV systems can be added here if needed\n")
    write(f, "! Format: New PVsystem.PV<bus> irrad=1.0 Phases=1 Bus1=<bus>.1 kV=$(round(base_kv_ln, digits=4)) kVA=<kVA> Pmpp=<Pmpp> %cutin=0.001 %cutout=0.001\n")
end

# ===== Generate Storage.dss (empty for now) =====
println("Generating Storage.dss (empty)...")
open(joinpath(output_dir, "Storage.dss"), "w") do f
    write(f, "! Storage systems can be added here if needed\n")
    write(f, "! Format: New Storage.Battery<bus> phases=1 Bus1=<bus> kV=$(round(base_kv_ln, digits=4)) kVA=<kVA> kWrated=<kW> kWhrated=<kWh> %stored=62.5 %reserve=30 %EffCharge=95 %EffDischarge=95 %idlingkW=0 Vminpu=0.95 Vmaxpu=1.05 DispMode=External\n")
end

# ===== Generate Capacitor.dss (empty for now) =====
println("Generating Capacitor.dss (empty)...")
open(joinpath(output_dir, "Capacitor.dss"), "w") do f
    write(f, "! Capacitor banks can be added here if needed\n")
    write(f, "! Format: New generator.<bus>C Phases=1 Bus1=<bus> kV=$(round(base_kv_ln, digits=4)) Model=1 kW=1e-06 kVAr=<kVAr>\n")
end

println("\n✓ Conversion complete!")
println("OpenDSS files written to: $output_dir")
println("\nGenerated files:")
println("  - Master.dss")
println("  - BranchData.dss")
println("  - Loads.dss")
println("  - PVSystem.dss")
println("  - Storage.dss")
println("  - Capacitor.dss")
