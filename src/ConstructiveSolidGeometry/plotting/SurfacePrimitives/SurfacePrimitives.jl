include("Plane.jl")
include("ConalPlane.jl")
include("CylindricalAnnulus.jl")
include("ConeMantle.jl")
include("ToroidalAnnulus.jl")
include("TorusMantle.jl")
include("SphereMantle.jl")

@recipe function f(g::AbstractSurfacePrimitive; SSD_style = :wireframe, n = 30)
    if SSD_style == :wireframe #update to only plot real surfaces
        linewidth --> 2
        for points in get_plot_points(g, n = n)
            @series begin
                label := ""
                points
            end
        end
    elseif SSD_style == :surface
        linewidth --> 0.2
        mesh(g)
    end
end
