"""
Test script to figure out how to extract Ipopt iteration count from JuMP/MOI
"""

using Pkg
Pkg.activate("envs/tadmm")  # Use tadmm environment

using JuMP
using Ipopt

# MOI is part of JuMP
const MOI = JuMP.MOI

println("="^70)
println("Testing Ipopt Iteration Extraction Methods")
println("="^70)

# Create a simple test problem
model = Model(Ipopt.Optimizer)
set_silent(model)  # Suppress Ipopt output for now

@variable(model, x >= 0)
@variable(model, y >= 0)
@objective(model, Min, x^2 + y^2)
@constraint(model, x + y >= 1)

println("\nSolving test problem...")
optimize!(model)
println("   Solution: x=$(value(x)), y=$(value(y)), obj=$(objective_value(model))")

println("\n1. Testing MOI.BarrierIterations() [Ipopt uses interior-point/barrier method]:")
try
    iters = MOI.get(model, MOI.BarrierIterations())
    println("   ✓ SUCCESS: BarrierIterations = $iters")
catch e
    println("   ✗ BarrierIterations failed: $(typeof(e))")
    println("      Error: $e")
end

println("\n2. Testing MOI.RawStatusString():")
try
    raw_status = MOI.get(model, MOI.RawStatusString())
    println("   Raw status string: '$raw_status'")
catch e
    println("   ERROR: $e")
end

println("\n3. Testing with more complex problem (to get >1 iteration):")
model2 = Model(Ipopt.Optimizer)
set_silent(model2)

n = 10
@variable(model2, x[1:n] >= 0)
@objective(model2, Min, sum(x[i]^2 for i=1:n))
@constraint(model2, sum(x) >= 10)
@constraint(model2, [i=1:n-1], x[i] + x[i+1] <= 3)

optimize!(model2)
println("   Objective: $(objective_value(model2))")

try
    iters = MOI.get(model2, MOI.BarrierIterations())
    println("   ✓ Complex problem: BarrierIterations = $iters")
catch e
    println("   ✗ BarrierIterations failed")
end

println("\n" * "="^70)
println("RECOMMENDATION:")
println("="^70)
if true  # We'll update this based on results
    println("Use: iters = MOI.get(model, MOI.BarrierIterations())")
    println("This should work for Ipopt since it uses interior-point algorithm")
else
    println("Need alternative method - MOI.BarrierIterations() doesn't work")
end
println("="^70)
