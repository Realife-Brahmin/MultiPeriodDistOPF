struct EqualityConstraints{nx, nu, nc, F1, F2, F3, F4, F5, F6}
    c::F1
    cx::F2
    cu::F3
    cxx::F4
    cux::F5
    cuu::F6
end

function EqualityConstraints(c::F, nx::Int64, nu::Int64) where F<:Function
    x::Vector{Num} = Symbolics.variables(:x, 1:nx)
    u::Vector{Num} = Symbolics.variables(:u, 1:nu)

    c_ = Symbolics.simplify(c(x, u))
    cx = Symbolics.simplify(Symbolics.jacobian(c_, x))
    cu = Symbolics.simplify(Symbolics.jacobian(c_, u))
    nc = length(c_)

    if nc > 0
        compile(a, args...) = Symbolics.build_function(a, args...;
            expression=Val(false), parallel=Symbolics.ShardedForm(), skipzeros=true)[1]
        c_func = compile(c_, x, u)
        cx_func = compile(cx, x, u)
        cu_func = compile(cu, x, u)
        
        ϕ::Vector{Num} = Symbolics.variables(:ϕ, 1:nc)  # adjoint for second-order tensor contraction

        cxx = Symbolics.simplify(Symbolics.hessian(c_' * ϕ, x))
        cux = Symbolics.simplify(Symbolics.jacobian(cu' * ϕ, x))
        cuu = Symbolics.simplify(Symbolics.hessian(c_' * ϕ, u))
        cxx_func = compile(cxx, x, u, ϕ)
        cux_func = compile(cux, x, u, ϕ)
        cuu_func = compile(cuu, x, u, ϕ)
    else
        c_func = nothing
        cx_func = nothing
        cu_func = nothing
        cxx_func = nothing
        cux_func = nothing
        cuu_func = nothing
    end

    return EqualityConstraints{nx, nu, nc, typeof(c_func), typeof(cx_func), typeof(cu_func),
                    typeof(cxx_func), typeof(cux_func), typeof(cuu_func)}(
                        c_func, cx_func, cu_func, cxx_func, cux_func, cuu_func)
end

function EqualityConstraints(nx::Int64, nu::Int64)
    return EqualityConstraints((x, u) -> [], nx, nu)
end
