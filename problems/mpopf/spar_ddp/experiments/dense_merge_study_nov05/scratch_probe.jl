using LinearAlgebra, Serialization
include(joinpath(ENV["REPLAY_REPO"],"ddp/examples/power_system/ieee123c_filterddp.jl"))
include(joinpath(ENV["REPLAY_REPO"],"ddp/examples/power_system/procs_pass.jl"))
w=PairWork(1020)
println("SCRATCH_COEFFICIENT_BYTES ",sum(sizeof,(w.K,w.Y1,w.Y2,w.g)))
println("SCRATCH_OBJECT_BYTES_AFTER_CONSTRUCTION ",Base.summarysize(w))
println("SCRATCH_FIELDS ",fieldnames(typeof(w)))
