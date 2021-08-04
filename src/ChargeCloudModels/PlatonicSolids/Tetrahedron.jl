struct Tetrahedron{T} <: AbstractChargeCloud
    points::SVector{4, CartesianPoint{T}}
end

function Tetrahedron(center::CartesianPoint{T}, length::T = T(1)) where {T}
    points = CartesianPoint{T}[center + CartesianVector{T}(0,0,length)]
    for φ in (0,120,240)
        push!(points, center + length * CartesianVector{T}(sqrt(8)/3*cosd(φ), sqrt(8)/3*sind(φ), -1/3))
    end
    Tetrahedron{T}( points )
end

get_vertices(::Type{Tetrahedron{T}}) where {T} = 4

@recipe function f(t::Tetrahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Tetrahedron"
        t.points
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            t.points[vcat(1,2,3,4,2)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            t.points[vcat(3,1,4)]
        end
    end
end
