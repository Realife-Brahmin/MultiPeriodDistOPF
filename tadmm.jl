#! Activate project environment
import Pkg; Pkg.activate(".")

using OpenDSSDirect
using Dates
using Printf
using Random

# Import the ODD parser
include("src/Parser/parseFromDSS.jl")
using .parseFromDSS

# --- USER DEFINED META VARIABLES ---
# System and simulation parameters
systemName = "ads10A_1ph"   # Change as needed
T = 24                      # Number of time steps

# Example: user can define LoadShape, PVShape, etc. here if needed
# For now, just placeholders
LoadShape = ones(T)         # Placeholder, replace with actual shape if needed
PVShape = ones(T)           # Placeholder

# --- PARSE SYSTEM ---
println("Parsing system '", systemName, "' with T=", T)
data = get_system_config_from_dss(systemName, T)

# --- INSPECT DATA ---
println("\n--- DATA DICT KEYS ---")
println(keys(data))

println("\n--- SAMPLE DATA CONTENTS ---")
for k in keys(data)
    println("[", k, "] => ", typeof(data[k]))
end

# Optionally, print a few actual values for inspection
println("\n--- EXAMPLE: Buses ---")
println(get(data, :Nset, "Nset not found"))

println("\n--- EXAMPLE: Branches ---")
println(get(data, :Lset, "Lset not found"))

println("\n--- EXAMPLE: Loads ---")
println(get(data, :NLset, "NLset not found"))

println("\n--- EXAMPLE: Batteries ---")
println(get(data, :Bset, "Bset not found"))

println("\n--- EXAMPLE: PVs ---")
println(get(data, :Dset, "Dset not found"))

println("\n--- EXAMPLE: Battery SOCs ---")
println(get(data, :B0, "B0 not found"))

# Add more inspection/printing as needed for your workflow

println("\n--- tadmm.jl skeleton complete ---")
