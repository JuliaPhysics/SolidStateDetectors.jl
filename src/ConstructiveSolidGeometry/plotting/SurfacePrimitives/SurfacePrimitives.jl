#include("ConeMantle.jl")
#include("EllipsoidMantle.jl")
#include("EllipticalSurface.jl")
#include("Plane.jl")
#include("Polygon.jl")
#include("TorusMantle.jl")

@recipe function f(s::Plane)
    throw(ArgumentError("No plot recipe defined for Plane."))
end

@recipe function f(s::AbstractSurfacePrimitive{T}; n_arc = 40, n_vert_lines = 2, n_samples = 40) where {T}
    seriestype --> :csg
    st = plotattributes[:seriestype]
    l = get_label_name(s)
    if st == :csg 
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
            linewidth --> 2
            mesh(s, n_arc, n_vert_lines)
        end 
    elseif st == :mesh3d
        label --> l   
        linecolor --> :white
        mesh(s, n_arc)
    elseif st == :wireframe
        label --> l
        seriescolor --> :black
        fillalpha := 1
        linewidth --> 2
        mesh(s, n_arc, n_vert_lines)
    elseif st == :samplesurface
        label --> l
        seriesalpha --> 0.4
        sample(s, extremum(s)/n_samples)
    elseif st == :slice
        @warn ":slice is not supported for AbstractSurfacePrimitive. Use :csg, :wireframe, :mesh3d, or :samplesurface."
    else
        @warn "The only seriestypes which will return a plot are :csg, :wireframe, :mesh3d, :samplesurface, and :slice."
    end  
end

@recipe function f(::Type{Val{:samplesurface}}, x, y, z)
    seriescolor --> 1
    seriesalpha --> 0.1
    markerstrokewidth --> 0
    markersize --> 3
    if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        aspect_ratio --> 1.0
    end 
    seriestype := :scatter
    ()
end

@recipe function f(::Type{Val{:slice}}, x, y, z)
    seriescolor --> 1
    seriesalpha --> 1
    markerstrokewidth --> 0
    markersize --> 1
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
