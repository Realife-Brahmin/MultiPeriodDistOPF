include("sequential_dss.jl")
instance,solution,out,efficiency=ARGS
SequentialDSS.replay(instance,solution,out;efficiency=efficiency)
