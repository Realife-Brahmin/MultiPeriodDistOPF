struct ControlLimits{T, nu}
    l::Vector{T}
    u::Vector{T}
    maskl::Vector{Bool}
    masku::Vector{Bool}
    nl::Int64
    nu::Int64
end

function ControlLimits(l::AbstractVector{T}, u::AbstractVector{T}) where T
    @assert length(l) == length(u)
    @assert all(u .>= l)
    maskl = (l .!= -floatmax(T)) .&& .!isinf.(l)
    masku = (u .!= floatmax(T)) .&& .!isinf.(u)
    nlo = sum(maskl)
    nup = sum(masku)
    nu = length(l)
    return ControlLimits{T, nu}(Vector{T}(l), Vector{T}(u),
        Vector{Bool}(maskl), Vector{Bool}(masku), nlo, nup)
end
