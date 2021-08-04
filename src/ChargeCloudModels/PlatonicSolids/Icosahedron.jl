struct Icosahedron{T} <: AbstractChargeCloud
    points::SVector{12, CartesianPoint{T}}
end

function Icosahedron(center::CartesianPoint{T}, length::T = T(1)) where {T}
    points = CartesianPoint{T}[]
    for i in (-1,1)
        for g in (-Base.MathConstants.golden, Base.MathConstants.golden)
            push!(points, center + length * CartesianVector{T}(i,0,g) / sqrt(T(2 + Base.MathConstants.golden)))
            push!(points, center + length * CartesianVector{T}(g,i,0) / sqrt(T(2 + Base.MathConstants.golden)))
            push!(points, center + length * CartesianVector{T}(0,g,i) / sqrt(T(2 + Base.MathConstants.golden)))
        end
    end
    Icosahedron{T}( points )
end

get_vertices(::Type{Icosahedron{T}}) where {T} = 12

@recipe function f(i::Icosahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Icosahedron"
        i.points
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(7,6,1,2,3,5,7,3,1,7,11,5,9,2,8,1)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(3,9,10,4,2,9,4,8,6,11,10,12,4)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(5,10,12,8)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(6,12,11)]
        end
    end
end


