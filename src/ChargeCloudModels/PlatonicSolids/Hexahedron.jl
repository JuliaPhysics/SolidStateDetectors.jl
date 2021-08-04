struct Hexahedron{T} <: AbstractChargeCloud
    points::SVector{8, CartesianPoint{T}}
end

function Hexahedron(center::CartesianPoint{T}, length::T = T(1)) where {T}
    points = CartesianPoint{T}[]
    for x in (-1,1)
        for y in (-1,1)
            for z in (-1,1)
                push!(points, center + length * sqrt(T(1/3)) * CartesianVector{T}(x,y,z))
            end
        end
    end
    Hexahedron{T}( points )
end

get_vertices(::Type{Hexahedron{T}}) where {T} = 8

@recipe function f(i::Hexahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Hexahedron"
        i.points
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(1,2,6,5,1,3,4,2)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(4,8,7,3)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(5,7)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(6,8)]
        end
    end
end