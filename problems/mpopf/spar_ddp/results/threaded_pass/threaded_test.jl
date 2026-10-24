# Development wrapper: the driver, with the threaded backward pass installed
# before it runs (the driver's own file is not touched while a sweep loads it).
const PS = joinpath(ENV["TP_REPO"], "ddp", "examples", "power_system")
include(joinpath(PS, "ieee123c_filterddp.jl"))      # defines main(); does not run it
include(joinpath(PS, "threaded_pass.jl"))
let system = ARGS[1], T = parse(Int, ARGS[2])
    ptag = haskey(ENV, "REDUCED_PROFILE") ? "_" * ENV["REDUCED_PROFILE"] : ""
    data = deserialize(joinpath(REPO, "ddp", "results", "network_filterddp", "network_data_$(system)_T$(T)$(ptag).jls"))
    idx, _ = control_layout(data)
    install_threaded_pass(data, idx)
end
main()
