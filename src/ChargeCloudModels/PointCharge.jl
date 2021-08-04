"""
    struct NBodyChargeCloud{T} <: AbstractChargeCloud

Struct which defines a single point-like charge carrier.

## Fields
* `points::SVector{1, CartesianPoint{T}}`: Position of the charge carrier, saved as single entry of a `Vector`.
"""
struct PointCharge{T} <: AbstractChargeCloud
    points::SVector{1, CartesianPoint{T}}
end

function PointCharge(center::Vector{CartesianPoint{T}}, length::T = T(0)) where {T} 
    PointCharge{T}(center)
end

get_vertices(::Type{PointCharge{T}}) where {T} = 1

@recipe function f(t::PointCharge{T}) where {T}

    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "PointCharge"
        t.points
    end
end