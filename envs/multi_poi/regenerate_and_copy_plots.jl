# Script to copy existing angle/voltage plots to LaTeX figures folder
# Run multi_poi_mpopf.jl first to generate the plots

println("="^80)
println("COPYING PLOTS TO LATEX FIGURES FOLDER")
println("="^80)

# Define source and destination directories using @__DIR__ for accurate paths
# @__DIR__ is envs/multi_poi/
plots_source = joinpath(@__DIR__, "processedData", "ieee123_5poi_1ph_T24", "plots")
# Go up to documents_general, then to documentsCreated
documents_general = joinpath(@__DIR__, "..", "..", "..")
latex_dest = joinpath(documents_general, "documentsCreated", "PESGM2026_Multi-Source-Multi-Period-OPF", "figures")

# Create destination if it doesn't exist
mkpath(latex_dest)

# Plot files to copy
plot_files = [
    "angle_voltage_slack1s.png",
    "angle_voltage_slack2s.png",
    "angle_voltage_slack3s.png",
    "angle_voltage_slack4s.png",
    "angle_voltage_slack5s.png"
]

println("\nSource: $plots_source")
println("Destination: $latex_dest")

for filename in plot_files
    src = joinpath(plots_source, filename)
    dst = joinpath(latex_dest, filename)
    
    if isfile(src)
        cp(src, dst, force=true)
        println("  ✓ Copied: $filename")
    else
        println("  ✗ Not found: $filename")
    end
end

println("\n" * "="^80)
println("DONE! All plots copied to LaTeX figures folder.")
println("="^80)
