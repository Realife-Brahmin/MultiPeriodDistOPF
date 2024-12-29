using JuMP
using Ipopt
include("src/SolverArranger/SolverArranger.jl")
import .SolverArranger as SolverArranger
# Define constants
eta_C = 0.95
eta_D = 0.95
B0_1 = 1.0
B0_2 = 1.0

# Define the first model
model1 = Model(Ipopt.Optimizer)
@variable(model1, B1 >= 0)
@variable(model1, B2 >= 0)
@variable(model1, Pc1 >= 0)
@variable(model1, Pd1 >= 0)
@variable(model1, Pc2 >= 0)
@variable(model1, Pd2 >= 0)

@objective(model1, Min, (B1 - 2)^2 + (B2 - 4)^2)
@constraint(model1, B1 - (B0_1 + eta_C * Pc1 - 1 / eta_D * Pd1) == 0)
@constraint(model1, B2 - (B0_2 + eta_C * Pc2 - 1 / eta_D * Pd2) == 0)

optimize!(model1)

println("Model 1 Results:")
println("B1 = ", value(B1))
println("B2 = ", value(B2))

# Copy the first model to create the second model
model2, _ = JuMP.copy_model(model1)

# Modify the objective of the second model
@objective(model2, Min, (model2[:B1] - 2)^2 + (model2[:B2] - 6)^2)
SolverArranger.attach_solver(model2, "Ipopt")
optimize!(model2)

println("Model 2 Results:")
println("B1 = ", value(model2[:B1]))
println("B2 = ", value(model2[:B2]))