#include("ConeMantle.jl")
#include("EllipsoidMantle.jl")
#include("EllipticalSurface.jl")
#include("Plane.jl")
#include("Polygon.jl")
#include("TorusMantle.jl")

@recipe function f(s::AbstractSurfacePrimitive; n_arc = 40, n_vert_lines = 2)
    seriestype --> :ssd
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :ssd 
        @series begin 
            label --> get_label_name(s)
            linecolor --> :black
            fillalpha := 1
            linewidth --> 1.5
            mesh(s, n_arc, n_vert_lines)
        end 
        @series begin 
            linewidth := 0
            linecolor --> :black
            mesh(s, n_arc)
        end
    elseif haskey(plotattributes, :seriestype) && plotattributes[:seriestype] in [:mesh3d, :mesh]
        label --> get_label_name(s)   
        linecolor --> :white
        mesh(s, n_arc)
    elseif haskey(plotattributes, :seriestype) && plotattributes[:seriestype] in [:path, :wireframe]
        label --> get_label_name(s)
        seriescolor --> :black
        fillalpha := 1
        linewidth --> 1
        mesh(s, n_arc, n_vert_lines)
    end    
end

@recipe function f(vp::AbstractVector{<:AbstractSurfacePrimitive})
    linecolor --> :black
    @series begin
        label --> "Faces"
        show_normal --> false
        vp[1]
    end
    if length(vp) > 1
        for p in vp[2:end]
            @series begin
                show_normal --> false
                label := nothing
                p
            end
        end
    end
end
