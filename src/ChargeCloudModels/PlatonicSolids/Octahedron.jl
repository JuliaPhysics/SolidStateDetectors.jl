struct Octahedron{T} <: AbstractChargeCloud
    points::SVector{6, CartesianPoint{T}}
end

function Octahedron(center::CartesianPoint{T}, length::T = T(1)) where {T}
    points = CartesianPoint{T}[
        center + length * CartesianVector{T}(1,0,0),
        center + length * CartesianVector{T}(-1,0,0),
        center + length * CartesianVector{T}(0,1,0),
        center + length * CartesianVector{T}(0,-1,0),
        center + length * CartesianVector{T}(0,0,1),
        center + length * CartesianVector{T}(0,0,-1)
    ]
    Octahedron{T}( points )
end

get_vertices(::Type{Octahedron{T}}) where {T} = 6

@recipe function f(i::Octahedron{T}) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Octahedron"
        i.points
    end
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :line3d
        label --> ""
        linewidth --> 1
        i.points[vcat(1,3,2,4,1,5,2,6,3,5,4,6,1)]
    end
end
