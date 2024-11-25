using JuMP
using Ipopt

# Create a model with the Ipopt solver
model = Model(Ipopt.Optimizer)

# Define the variable
@variable(model, x)

# Define the objective function (example: minimize x^2)
@objective(model, Min, x^2)

# Add the constraint x^2 - x = 0
@constraint(model, con, x^2 - x == 0)

# Optimize the model
optimize!(model)

# Retrieve the value of the variable
optimal_x = value(x)
println("Optimal value of x: ", optimal_x)

# Retrieve the dual variable of the constraint
dual_con = dual(con)
println("Dual variable of the constraint: ", dual_con)