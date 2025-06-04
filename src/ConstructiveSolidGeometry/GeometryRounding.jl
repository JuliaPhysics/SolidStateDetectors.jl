csg_round_lin(x::T; sigdigits::Int = 9) where {T} = round(x, sigdigits = sigdigits) # min. Δx = 1 nm 
csg_round_lin(x::Float32; sigdigits::Int = 5) = round(x, sigdigits = sigdigits) # min. Δx = 10 μm 
csg_round_rad(φ::T; digits::Int = 4) where {T} = round(φ, digits = digits) # min. Δφ ≈ 0.011 °

csg_default_tol(::Type{T}) where {T} = T(1e-12)
csg_default_tol(::Type{Float32}) = 1f-8

csg_isapprox(x::T, y::T; csgtol::T = T(1e-9)) where {T} = abs(x - y) < csgtol


geom_sigdigits(::Type{Int64})::Int = 12
geom_sigdigits(::Type{Float32})::Int = 6
geom_sigdigits(::Type{Float64})::Int = 12

geom_atol_zero(::Type{Int64})::Float64 = 1e-12
geom_atol_zero(::Type{Float32})::Float32 = 1e-6
geom_atol_zero(::Type{Float64})::Float64 = 1e-12

function geom_round(v::T)::T where {T}
    return isapprox(v, 0, atol = geom_atol_zero(T)) ? T(0) : round(v, sigdigits = geom_sigdigits(T))
end

function geom_round(pt::CylindricalPoint{T})::CylindricalPoint{T} where {T <: Real}
    return CylindricalPoint{T}( geom_round(pt.r), geom_round(pt.φ), geom_round(pt.z)  )
end

function geom_round(pt::CartesianPoint{T})::CartesianPoint{T} where {T <: Real}
    return CartesianPoint{T}( geom_round(pt.x), geom_round(pt.y), geom_round(pt.z)  )
end

function geom_round(vpt::Vector{CP}) where {CP <: AbstractCoordinatePoint}
    [geom_round(pt) for pt in vpt]
end