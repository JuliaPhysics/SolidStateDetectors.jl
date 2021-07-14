include("Polygon.jl")
include("EllipticalSurface.jl")
include("ConeMantle.jl")
include("EllipsoidMantle.jl")
include("TorusMantle.jl")

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
