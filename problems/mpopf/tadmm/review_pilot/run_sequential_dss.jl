include("sequential_dss.jl")
instance,solution,out,efficiency=ARGS
SequentialDSS.replay(instance,solution,out;efficiency=efficiency)
if get(ENV,"SEQUENTIAL_DSS_GZIP","0")=="1"
    for file in readdir(out;join=true)
        endswith(file,".csv") && run(`gzip $file`)
    end
end
