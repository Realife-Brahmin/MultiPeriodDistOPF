# The driver with the paired backward pass installed (paired_pass.jl), until the
# driver has its own switch for it. Usage: as ieee123c_filterddp.jl.
const PS = @__DIR__
include(joinpath(PS, "ieee123c_filterddp.jl"))      # defines main(); does not run it
include(joinpath(PS, "paired_pass.jl"))
let system = ARGS[1], T = parse(Int, ARGS[2])
    ptag = haskey(ENV, "REDUCED_PROFILE") ? "_" * ENV["REDUCED_PROFILE"] : ""
    data = deserialize(joinpath(REPO, "spar_ddp", "results", "network_filterddp", "network_data_$(system)_T$(T)$(ptag).jls"))
    idx, _ = control_layout(data)
    install_paired_pass(data, idx)
end
main()
