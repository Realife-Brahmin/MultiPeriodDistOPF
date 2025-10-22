import Pkg; Pkg.activate(".")
println("Testing tadmm.jl loading...")

# Try to load the file and check for errors
try
    include("tadmm.jl")
    println("✓ Script loaded successfully!")
catch e
    println("✗ Error loading script:")
    println(e)
    rethrow(e)
end
