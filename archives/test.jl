using JuMP

function display_all_constraints(model, filename::String)
    open(filename, "w") do file
        for (i, con) in enumerate(all_constraints(model, include_variable_in_set_constraints=true))
            println(file, "Constraint $i: $con")
        end
    end
end

# Example usage
@unpack model = modelDict
display_all_constraints(model, "constraints_output.txt")