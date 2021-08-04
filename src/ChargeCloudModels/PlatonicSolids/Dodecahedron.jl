struct Dodecahedron{T} <: AbstractChargeCloud
    points::SVector{20, CartesianPoint{T}}
end

function Dodecahedron(center::CartesianPoint{T}, length::T = T(1)) where {T}
    points = CartesianPoint{T}[
        center + sqrt(1/3) * length * CartesianVector{T}(-1,-1,-1),
        center + sqrt(1/3) * length * CartesianVector{T}(-1,-1,1),
        center + sqrt(1/3) * length * CartesianVector{T}(-1,1,-1),
        center + sqrt(1/3) * length * CartesianVector{T}(1,-1,-1),
        center + sqrt(1/3) * length * CartesianVector{T}(-1,1,1),
        center + sqrt(1/3) * length * CartesianVector{T}(1,-1,1),
        center + sqrt(1/3) * length * CartesianVector{T}(1,1,-1),
        center + sqrt(1/3) * length * CartesianVector{T}(1,1,1),
    ]
    for g in (-Base.MathConstants.golden, Base.MathConstants.golden)
        for invg in (-inv(Base.MathConstants.golden), inv(Base.MathConstants.golden))
            push!(points, center + sqrt(1/3) * length * CartesianVector{T}(invg, g, 0))
            push!(points, center + sqrt(1/3) * length * CartesianVector{T}(0, invg, g))
            push!(points, center + sqrt(1/3) * length * CartesianVector{T}(g, 0, invg))
        end
    end
    Dodecahedron{T}( points )
end

get_vertices(::Type{Dodecahedron{T}}) where {T} = 20

@recipe function f(i::Dodecahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Dodecahedron"
        i.points
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(1,9,2,14,5,15,3,11,1,10,4,12,6,16,2)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(5,19,16)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(19,8,18,15)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(3,13,10)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(4,17,7,13)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(7,18)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(8,20,17)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(6,20)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(9,12)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(11,14)]
        end
    end
end