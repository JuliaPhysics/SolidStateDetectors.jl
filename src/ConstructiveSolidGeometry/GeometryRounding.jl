csg_default_tol(::Type{T}) where {T} = T(1e-12)
csg_default_tol(::Type{Float32}) = 1f-8

geom_sigdigits(::Type{Float32})::Int = 6
geom_sigdigits(::Type{Float64})::Int = 12

geom_atol_zero(::Type{Float32})::Float32 = 1e-6
geom_atol_zero(::Type{Float64})::Float64 = 1e-12

function geom_round(v::T)::T where {T <: Union{Float32, Float64}}
    return isapprox(v, 0, atol = geom_atol_zero(T)) ? T(0) : round(v, sigdigits = geom_sigdigits(T))
end

function geom_round(pt::CylindricalPoint{T})::CylindricalPoint{T} where {T <: Union{Float32, Float64}}
    return CylindricalPoint{T}( geom_round(pt.r), geom_round(pt.φ), geom_round(pt.z)  )
end

function geom_round(pt::CartesianPoint{T})::CartesianPoint{T} where {T <: Union{Float32, Float64}}
    return CartesianPoint{T}( geom_round(pt.x), geom_round(pt.y), geom_round(pt.z)  )
end

function geom_round(vpt::AbstractVector{<:AbstractCoordinatePoint})
    geom_round.(vpt)
end