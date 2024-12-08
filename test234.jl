using JuMP, Ipopt, Test
model = Model(Ipopt.Optimizer)
set_silent(model)
@variable(model, x >= 1)
@objective(model, Min, x + 0.5)
x_vals = Float64[]
function my_callback(
    alg_mod::Cint,
    iter_count::Cint,
    obj_value::Float64,
    inf_pr::Float64,
    inf_du::Float64,
    mu::Float64,
    d_norm::Float64,
    regularization_size::Float64,
    alpha_du::Float64,
    alpha_pr::Float64,
    ls_trials::Cint,
)
    push!(x_vals, callback_value(model, x))
    @test isapprox(obj_value, 1.0 * x_vals[end] + 0.5, atol=1e-1)
    # return `true` to keep going, or `false` to terminate the optimization.
    return iter_count < 1
end
MOI.set(model, Ipopt.CallbackFunction(), my_callback)
optimize!(model)
@test MOI.get(model, MOI.TerminationStatus()) == MOI.INTERRUPTED
@test length(x_vals) == 2