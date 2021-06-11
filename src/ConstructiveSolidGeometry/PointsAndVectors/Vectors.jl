
struct PlanarVector{T} <: AbstractPlanarVector{T}
    u::T
    v::T
end

"""
    struct CartesianVector{T} <: AbstractCoordinateVector{T, Cartesian}

* `x`: x-component in meter
* `y`: y-component in meter
* `z`: z-component in meter
"""
struct CartesianVector{T} <: AbstractCoordinateVector{T, Cartesian}
    x::T
    y::T
    z::T
end

zero(VT::AbstractCoordinateVector{T}) where {T} = VT(zero(T),zero(T),zero(T))

scale(v::CartesianVector{T}, s::SVector{3}) where {T} = CartesianVector{T}(v .* s)

@inline rotate(p::CartesianPoint{T}, r::RotMatrix{3,T,TT}) where {T, TT} = r.mat * p
@inline rotate(p::CylindricalPoint{T}, r::RotMatrix{3,T,TT}) where {T, TT} = CylindricalPoint(rotate(CartesianPoint(p), r))
@inline rotate!(vp::Vector{<:AbstractCoordinatePoint{T}}, r::RotMatrix{3,T,TT}) where {T, TT} = begin for i in eachindex(vp) vp[i] = rotate(vp[i], r) end; vp end
@inline rotate!(vvp::Vector{<:Vector{<:AbstractCoordinatePoint{T}}}, r::RotMatrix{3,T,TT}) where {T, TT} = begin for i in eachindex(vvp) rotate!(vvp[i], r) end; vvp end

@inline scale(p::CartesianPoint{T}, v::SVector{3,T}) where {T} = p .* v
@inline scale(p::CylindricalPoint{T}, v::SVector{3,T}) where {T} = CylindricalPoint(scale(CartesianPoint(p),v))
@inline scale!(vp::Vector{<:AbstractCoordinatePoint{T}}, v::SVector{3,T}) where {T} = begin for i in eachindex(vp) vp[i] = scale(vp[i], v) end; vp end
@inline scale!(vvp::Vector{<:Vector{<:AbstractCoordinatePoint{T}}}, v::SVector{3,T}) where {T} = begin for i in eachindex(vvp) scale!(vvp[i], v) end; vvp end

@inline translate(p::CartesianPoint{T}, v::CartesianVector{T}) where {T} = p + v
@inline translate(p::CylindricalPoint{T}, v::CartesianVector{T}) where {T} = CylindricalPoint(translate(CartesianPoint(p), v))
@inline translate!(vp::Vector{<:AbstractCoordinatePoint{T}}, v::CartesianVector{T}) where {T} = begin for i in eachindex(vp) vp[i] = translate(vp[i], v) end; vp end
@inline translate!(vvp::Vector{<:Vector{<:AbstractCoordinatePoint{T}}}, v::CartesianVector{T}) where {T} =  begin for i in eachindex(vvp) translate!(vvp[i], v) end; vvp end


"""
    struct CylindricalVector{T} <: AbstractCoordinateVector{T, Cylindrical}

* `r`: Radius in meter
* `φ`: Polar angle in radians. φ == 0 <=> Parallel to x-axis of cartesian coordinate system."
* `z`: z-coordinate in meter
"""
struct CylindricalVector{T} <: AbstractCoordinateVector{T, Cylindrical}
    r::T
    φ::T
    z::T
end

CartesianVector(v::CylindricalVector{T}) where {T} = 
    CartesianVector{T}(CartesianPoint(CylindricalPoint{T}(v)))