# Parse only: no includes, evaluation, dependencies, simulations or numerical tests.
function main()
    root = normpath(joinpath(@__DIR__, ".."))
    count = 0
    function inspect(ex, path)
        ex isa Expr || return
        if ex.head in (:error, :incomplete)
            error("Julia syntax error in $path: $(ex.args)")
        end
        for arg in ex.args
            inspect(arg, path)
        end
    end
    for folder in ("spar_ddp/DDP4OPF.jl/src", "spar_ddp/examples", "spar_ddp/experiments", "tools")
        for (dir, _, files) in walkdir(joinpath(root, folder))
            for file in files
                endswith(file, ".jl") || continue
                path = joinpath(dir, file)
                inspect(Meta.parseall(read(path, String); filename=path), path)
                count += 1
            end
        end
    end
    println("PASS: $count Julia files parsed; no code evaluated or numerical claims validated.")
end
main()
