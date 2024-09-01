#include("ConeMantle.jl")
#include("EllipsoidMantle.jl")
#include("EllipticalSurface.jl")
#include("Plane.jl")
#include("Polygon.jl")
#include("TorusMantle.jl")

@recipe function f(s::Plane)
    throw(ArgumentError("No plot recipe defined for Plane."))
end

@recipe function f(s::AbstractSurfacePrimitive{T}; n_arc = 40, n_vert_lines = 2, n_samples = 40, slice_val = T(0)) where {T}
    seriestype --> :csg
    l = get_label_name(s)
    if haskey(plotattributes, :seriestype) 
        projections = [:x, :y, :z]
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
                linewidth --> 2
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
        elseif plotattributes[:seriestype] in projections
            label --> l
            spacing = extremum(s)/n_samples
            samples = filter(pt -> abs(getproperty(pt, plotattributes[:seriestype]) - slice_val) < spacing/2, sample(s, spacing))
            proj = filter(x -> x != plotattributes[:seriestype], projections)
            xguide --> string(proj[1])
            xunit --> internal_length_unit
            yguide --> string(proj[2])
            yunit --> internal_length_unit
            unitformat --> :slash
            internal_length_unit*getproperty.(samples, proj[1]), internal_length_unit*getproperty.(samples, proj[2])
        else
            @warn "The only seriestypes which will return a plot are :csg, :wireframe, :mesh3d, :samplesurface, and :x, :y, or :z."
        end
    end    
end

@recipe function f(::Type{Val{:samplesurface}}, x, y, z)
    seriescolor --> 1
    seriesalpha --> 0.2
    markerstrokewidth --> 0
    markersize --> 4
    if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        aspect_ratio --> 1.0
    end 
    seriestype := :scatter
    ()
end

@recipe function f(::Union{Type{Val{:x}}, Type{Val{:y}}, Type{Val{:z}}}, x, y, z)
    seriescolor --> 1
    seriesalpha --> 1
    markerstrokewidth --> 0
    markersize --> 1
    if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        aspect_ratio --> 1.0
    end 
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
