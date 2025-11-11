# opendss_validation.jl
# Validates OPF results against OpenDSS simulation
# Extracted from test_multi_poi.jl for cleaner separation

using Crayons
import OpenDSSDirect as dss
using Printf

"""
    validate_opf_with_opendss(modelDict, data)

Validates the OPF solution by:
1. Running delta sweep in OpenDSS
2. Comparing with BFM-NL results
3. Running generator simulation for verification
"""
function validate_opf_with_opendss(modelDict, data)
    systemName = "small2poi_1ph"
    rawDataFolder = "rawData/"
    systemFile = rawDataFolder * systemName * "/Master.dss"
    
    kVA_B = data[:kVA_B]
    
    println("\n" * "="^80)
    println("OPENDSS VALIDATION")
    println("="^80)
    
    dss.Text.Command("Clear")
    dss.Text.Command("Redirect " * systemFile)
    
    # Define delta sweep
    deltas = vcat(
        [-8, -5],
        collect(-2:0.5:-0.5),
        collect(-0.4:0.1:0.4),
        collect(0.5:0.5:2),
        [5, 8],
        [0.12452868971865128],
    )
    
    PSubs1 = []
    QSubs1 = []
    PSubs2 = []
    QSubs2 = []
    
    # Run delta sweep
    for delta in deltas
        dss.Text.Command("Edit Vsource.grid2 angle=" * string(-delta))
        dss.Text.Command("Solve")
        
        dss.Circuit.SetActiveElement("Vsource.grid1")
        powers1 = dss.CktElement.Powers()
        push!(PSubs1, -real(powers1[1]))
        push!(QSubs1, -imag(powers1[1]))
        
        dss.Circuit.SetActiveElement("Vsource.grid2")
        powers2 = dss.CktElement.Powers()
        push!(PSubs2, -real(powers2[1]))
        push!(QSubs2, -imag(powers2[1]))
    end
    
    # Get system information
    Vsource_names = sort(setdiff(dss.Vsources.AllNames(), ["source"]))
    Vsource_pus = Dict{String, Float64}()
    Vsource_basekv = Dict{String, Float64}()
    for vname in Vsource_names
        dss.Vsources.Name(vname)
        Vsource_pus[vname] = dss.Vsources.PU()
        Vsource_basekv[vname] = dss.Vsources.BasekV()
    end
    
    # Get total loading
    load_names = sort(dss.Loads.AllNames())
    total_P = 0.0
    total_Q = 0.0
    for lname in load_names
        dss.Loads.Name(lname)
        total_P += dss.Loads.kW()
        total_Q += dss.Loads.kvar()
    end
    P_L_kW = total_P
    Q_L_kVAr = total_Q
    
    # Print results table
    header = @sprintf("%-6s | %-12s | %-12s | %-12s | %-12s", 
                     "δ [°]", "PSubs1 [kW]", "PSubs2 [kW]", "QSubs1 [kVAr]", "QSubs2 [kVAr]")
    separator = "-"^length(header)
    
    println("\nDelta Sweep Results:")
    println("Total loading: $(round(P_L_kW, digits=2)) kW, $(round(Q_L_kVAr, digits=2)) kVAr")
    println("Substation voltages (pu): " *
        join(["$(name)=$(@sprintf("%.2f", Vsource_pus[name]))" for name in Vsource_names], ", "))
    println("Substation BasekV: " *
        join(["$(name)=$(@sprintf("%.4f", Vsource_basekv[name]))" for name in Vsource_names], ", "))
    println(separator)
    println(header)
    println(separator)
    
    # Print delta sweep rows
    for i in eachindex(deltas)
        p1 = PSubs1[i]
        p2 = PSubs2[i]
        p1_str = p1 > 0 ? Crayon(foreground = :green)(@sprintf("%-12.2f", p1)) : 
                          Crayon(foreground = :red)(@sprintf("%-12.2f", p1))
        p2_str = p2 > 0 ? Crayon(foreground = :green)(@sprintf("%-12.2f", p2)) : 
                          Crayon(foreground = :red)(@sprintf("%-12.2f", p2))
        delta_str = (p1 < 0 || p2 < 0) ? Crayon(foreground = :red)(@sprintf("%-6.4f", deltas[i])) : 
                                          @sprintf("%-6.4f", deltas[i])
        @printf("%s | %s | %s | %-12.2f | %-12.2f\n",
            delta_str, p1_str, p2_str, QSubs1[i], QSubs2[i])
    end
    
    # Print BFM-NL optimal row
    opt_p1 = modelDict[:P_1j] * kVA_B
    opt_p2 = modelDict[:P_2j] * kVA_B
    opt_q1 = modelDict[:Q_1j] * kVA_B
    opt_q2 = modelDict[:Q_2j] * kVA_B
    opt_crayon = Crayon(foreground = :white, bold = true)
    opt_delta_str = opt_crayon(@sprintf("%-6s", "'BFM-NL'"))
    opt_p1_str = opt_crayon(@sprintf("%-12.2f", opt_p1))
    opt_p2_str = opt_crayon(@sprintf("%-12.2f", opt_p2))
    opt_q1_str = opt_crayon(@sprintf("%-12.2f", opt_q1))
    opt_q2_str = opt_crayon(@sprintf("%-12.2f", opt_q2))
    @printf("%s | %s | %s | %s | %s\n",
        opt_delta_str, opt_p1_str, opt_p2_str, opt_q1_str, opt_q2_str)
    
    println("-"^length(header))
    
    # -------------------- GEN2 SIMULATION --------------------
    println("\nGenerator Simulation (Verification):")
    
    # Set up OpenDSS for gen2 simulation
    dss.Text.Command("Edit Vsource.grid2 enabled=no")
    dss.Text.Command("Edit Generator.gen2 enabled=yes")
    dss.Text.Command(@sprintf("Edit Generator.gen2 enabled=yes kW=%.6f kvar=%.6f", 
                             modelDict[:P_2j] * kVA_B, modelDict[:Q_2j] * kVA_B))
    dss.Text.Command("Solve")
    
    # Retrieve angles
    dss.Vsources.Name("grid1")
    dss.Circuit.SetActiveElement("Vsource.grid1")
    vmang_grid1 = dss.CktElement.VoltagesMagAng()
    angle_1s = vmang_grid1[2, 1]
    
    dss.Generators.Name("gen2")
    dss.Circuit.SetActiveElement("Generator.gen2")
    vmang_gen2 = dss.CktElement.VoltagesMagAng()
    angle_2s = vmang_gen2[2, 1]
    delta_12 = angle_1s - angle_2s
    
    # Retrieve power dispatches
    dss.Circuit.SetActiveElement("Vsource.grid1")
    powers_grid1 = dss.CktElement.Powers()
    p_grid1 = -real(powers_grid1[1])
    q_grid1 = -imag(powers_grid1[1])
    
    dss.Circuit.SetActiveElement("Generator.gen2")
    powers_gen2 = dss.CktElement.Powers()
    p_gen2 = -real(powers_gen2[1])
    q_gen2 = -imag(powers_gen2[1])
    
    # Print gen2 simulation results
    header_gen2 = @sprintf("%-6s | %-12s | %-12s | %-12s | %-12s", 
                          "δ [°]", "PSubs1 [kW]", "PGen2 [kW]", "QSubs1 [kVAr]", "QGen2 [kVAr]")
    println(header_gen2)
    println("-"^length(header_gen2))
    
    gen2_row_crayon = Crayon(foreground = :cyan, bold = true)
    gen2_delta_str = gen2_row_crayon(@sprintf("%-6.4f", delta_12))
    gen2_p1_str = gen2_row_crayon(@sprintf("%-12.2f", p_grid1))
    gen2_p2_str = gen2_row_crayon(@sprintf("%-12.2f", p_gen2))
    gen2_q1_str = gen2_row_crayon(@sprintf("%-12.2f", q_grid1))
    gen2_q2_str = gen2_row_crayon(@sprintf("%-12.2f", q_gen2))
    @printf("%s | %s | %s | %s | %s\n",
        gen2_delta_str, gen2_p1_str, gen2_p2_str, gen2_q1_str, gen2_q2_str)
    println(separator)
    
    # Save to file
    output_dir = "processedData/" * systemName * "/"
    isdir(output_dir) || mkpath(output_dir)
    output_file = output_dir * "substation-power-distribution-load_$(Int(round(P_L_kW)))_kW.txt"
    
    open(output_file, "w") do io
        println(io, "Total loading: $(round(P_L_kW, digits=2)) kW, $(round(Q_L_kVAr, digits=2)) kVAr")
        println(io, "Substation voltages (pu): " *
            join(["$(name)=$(@sprintf("%.2f", Vsource_pus[name]))" for name in Vsource_names], ", "))
        println(io, "Substation BasekV: " *
            join(["$(name)=$(@sprintf("%.4f", Vsource_basekv[name]))" for name in Vsource_names], ", "))
        println(io, header)
        println(io, separator)
        for i in eachindex(deltas)
            @printf(io, "%-6.4f | %-12.2f | %-12.2f | %-12.2f | %-12.2f\n",
                deltas[i], PSubs1[i], PSubs2[i], QSubs1[i], QSubs2[i])
        end
        @printf(io, "%-6s | %-12.2f | %-12.2f | %-12.2f | %-12.2f\n",
            "'BFM-NL'", opt_p1, opt_p2, opt_q1, opt_q2)
        println(io, "-"^length(header))
        println(io, header_gen2)
        println(io, "-"^length(header_gen2))
        @printf(io, "%-6.4f | %-12.2f | %-12.2f | %-12.2f | %-12.2f\n",
            delta_12, p_grid1, p_gen2, q_grid1, q_gen2)
        println(io, separator)
        println(io, "$(round(delta_12, digits=2)) deg")
    end
    
    println("\nResults saved to: $output_file")
    println("="^80)
    
    return Dict(
        :delta_12 => delta_12,
        :p_grid1 => p_grid1,
        :p_gen2 => p_gen2,
        :q_grid1 => q_grid1,
        :q_gen2 => q_gen2,
        :deltas => deltas,
        :PSubs1 => PSubs1,
        :PSubs2 => PSubs2,
        :QSubs1 => QSubs1,
        :QSubs2 => QSubs2
    )
end
