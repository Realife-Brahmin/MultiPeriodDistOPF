struct Dynamics{nx, nu, F1, F2, F3, F4, F5, F6}
    f::F1
    fx::F2
    fu::F3
    fxx::F4
    fux::F5
    fuu::F6
end

function Dynamics(f::F, nx::Int64, nu::Int64) where F<:Function
    x::Vector{Num} = Symbolics.variables(:x, 1:nx)
    u::Vector{Num} = Symbolics.variables(:u, 1:nu)

    y = Symbolics.simplify(f(x, u))
    λ::Vector{Num} = Symbolics.variables(:λ, 1:nx)  # vector variables for Hessian vector products
    
    fx = Symbolics.simplify(Symbolics.jacobian(y, x))
    fu = Symbolics.simplify(Symbolics.jacobian(y, u))
    compile(a, args...) = Symbolics.build_function(a, args...;
        expression=Val(false), parallel=Symbolics.ShardedForm(), skipzeros=true)[1]
    f_func = compile(y, x, u)
    fx_func = compile(fx, x, u)
    fu_func = compile(fu, x, u)
    
    fxx = Symbolics.simplify(Symbolics.hessian(λ' * y, x))
    fux = Symbolics.simplify(Symbolics.jacobian(Symbolics.gradient(λ' * y, u), x))
    fuu = Symbolics.simplify(Symbolics.hessian(λ' * y, u))
    fxx_func = compile(fxx, x, u, λ)
    fux_func = compile(fux, x, u, λ)
    fuu_func = compile(fuu, x, u, λ)

    return Dynamics{nx, nu, typeof(f_func), typeof(fx_func), typeof(fu_func),
                    typeof(fxx_func), typeof(fux_func), typeof(fuu_func)}(
                        f_func, fx_func, fu_func, fxx_func, fux_func, fuu_func)
end
