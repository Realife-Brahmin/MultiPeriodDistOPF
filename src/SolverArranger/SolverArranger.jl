module SolverArranger

export 
    attach_solver,
    configure_solver

using JuMP
using Gurobi
using Ipopt

#region attach_solver
"""
    attach_solver(model, solver_name)

Attach the specified solver to the optimization model.

This function sets the optimizer for the given `model` based on the `solver_name`. 
It configures the optimizer attributes and sets the model to silent mode if needed.
"""
function attach_solver(model, solver_name)
    if solver_name == "Ipopt"
        optimizer = Ipopt.Optimizer
        optimizer_attributes = Dict("tol" => 1e-6, "max_iter" => 10000)
    elseif solver_name == "Gurobi"
        optimizer = Gurobi.Optimizer
        optimizer_attributes = Dict("TimeLimit" => 300)
    elseif solver_name == "Juniper"
        ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
        optimizer = optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt)
        optimizer_attributes = Dict()
    elseif solver_name == "EAGO"
        optimizer = EAGO.Optimizer
        optimizer_attributes = Dict()
    elseif solver_name == "MadNLP"
        optimizer = MadNLP.Optimizer
        optimizer_attributes = Dict()
    else
        error("Unsupported solver")
    end

    # Set the optimizer for the model
    set_optimizer(model, optimizer)

    # Set optimizer attributes
    for (attr, value) in optimizer_attributes
        set_optimizer_attribute(model, attr, value)
    end

    # Optionally, set the model to silent if needed
    set_silent(model)

    return model
end
#endregion

#region configure_solver
"""
    configure_solver(solver_name; gurobi_env=nothing)

Configure and return an optimization model with the specified solver.

This function creates and configures an optimization model based on the `solver_name`. 
It sets the optimizer attributes and returns the configured model.
For Gurobi, optionally accepts a gurobi_env to suppress license messages.
"""
function configure_solver(solver_name; gurobi_env=nothing)
    if solver_name == "Ipopt"
        model = Model(Ipopt.Optimizer)
        set_silent(model)
        set_optimizer_attribute(model, "tol", 1e-6)
        set_optimizer_attribute(model, "max_iter", 10000)
        # set_optimizer_attribute(model, "print_level", 5)
    elseif solver_name == "Gurobi"
        if gurobi_env !== nothing
            model = Model(() -> Gurobi.Optimizer(gurobi_env))
        else
            model = Model(Gurobi.Optimizer)
        end
        set_optimizer_attribute(model, "OutputFlag", 0)  # Suppress Gurobi output
        set_optimizer_attribute(model, "TimeLimit", 300)        # Limit time (in seconds)
    elseif solver == "Juniper"
        ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
        optimizer = optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt)
        model = Model(optimizer)
    elseif solver == "EAGO"
        model = Model(EAGO.Optimizer)
    elseif solver == "MadNLP"
        model = Model(MadNLP.Optimizer)
    else
        error("Unsupported solver")
    end

    return model
end
#endregion

end # module