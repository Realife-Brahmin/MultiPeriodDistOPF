struct Objective{nx, nu, O1, O2, O3, O4, O5, O6}
    l::O1
    lx::O2
    lu::O3
    lxx::O4
    lux::O5
    luu::O6
end

function Objective(l::F, nx::Int64, nu::Int64) where F<:Function
    x::Vector{Num} = Symbolics.variables(:x, 1:nx)
    u::Vector{Num} = Symbolics.variables(:u, 1:nu)

    l_ = Symbolics.simplify(l(x, u))
    lx = Symbolics.simplify(Symbolics.gradient(l_, x))
    lu = Symbolics.simplify(Symbolics.gradient(l_, u))
    lxx = Symbolics.simplify(Symbolics.jacobian(lx, x))
    lux = Symbolics.simplify(Symbolics.jacobian(lu, x))
    luu = Symbolics.simplify(Symbolics.jacobian(lu, u))

    # Compile ordinary sharded functions. RuntimeGeneratedFunction cannot shard
    # large array expressions and fails on network-scale dense Hessians.
    compile(a, args...) = Symbolics.build_function(a, args...;
        expression=Val(false), parallel=Symbolics.ShardedForm(), skipzeros=true)[1]
    l_func = compile([l_], x, u)
    lx_func = compile(lx, x, u)
    lu_func = compile(lu, x, u)
    lxx_func = compile(lxx, x, u)
    lux_func = compile(lux, x, u)
    luu_func = compile(luu, x, u)
    
    return Objective{nx, nu, typeof(l_func), typeof(lx_func), typeof(lu_func),
                     typeof(lxx_func), typeof(lux_func), typeof(luu_func)}(
                        l_func, lx_func, lu_func, lxx_func, lux_func, luu_func)
end
