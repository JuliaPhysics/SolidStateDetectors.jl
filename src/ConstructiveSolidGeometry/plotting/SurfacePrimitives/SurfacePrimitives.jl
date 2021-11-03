#include("ConeMantle.jl")
#include("EllipsoidMantle.jl")
#include("EllipticalSurface.jl")
#include("Plane.jl")
#include("Polygon.jl")
#include("TorusMantle.jl")

@recipe function f(s::AbstractSurfacePrimitive{T}; n_arc = 40, n_vert_lines = 2, n_samples = 100) where {T}
    seriestype --> :csg
    l = get_label_name(s)
    if haskey(plotattributes, :seriestype) 
        if plotattributes[:seriestype] == :csg 
            @series begin 
                label := ""
                linewidth := 0
                linecolor --> :white
                mesh(s, n_arc)
            end
            @series begin 
                label --> l
                linecolor --> :black
                fillalpha := 1
                linewidth --> 1.5
                mesh(s, n_arc, n_vert_lines)
            end 
        elseif plotattributes[:seriestype] == :mesh3d
            label --> l   
            linecolor --> :white
            mesh(s, n_arc)
        elseif plotattributes[:seriestype] == :wireframe
            label --> l
            seriescolor --> :black
            fillalpha := 1
            linewidth --> 2
            mesh(s, n_arc, n_vert_lines)
        elseif plotattributes[:seriestype] == :samplesurface
            label --> l
            seriesalpha --> 0.4
            sample(s, extremum(s)/n_samples)
        else
            @warn "The only seriestypes wich will return a plot are :csg, :wireframe, :mesh3d, and :samplesurface"
        end
    end    
end

@recipe function f(::Type{Val{:samplesurface}}, x, y, z)
    seriescolor --> 1
    seriesalpha --> 0.2
    markerstrokewidth --> 0
    markersize --> 2
    seriestype := :scatter
    ()
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
