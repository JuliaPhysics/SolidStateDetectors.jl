function mesh(es::EllipticalSurface{T}; n = 30) where {T <: AbstractFloat}

    rMin::T, rMax::T = get_r_limits(es)
    φMin::T, φMax::T = get_φ_limits(es)
    f = (φMax - φMin)/(2π)
    n = Int(ceil(n*f))
    φ = range(φMin, φMax, length = n+1)
    r = range(rMin, rMax, length = 2)
    z = [0,0]

    X::Array{T,2} = [r_j*cos(φ_i) for φ_i in φ, r_j in r]
    Y::Array{T,2} = [r_j*sin(φ_i) for φ_i in φ, r_j in r]
    Z::Array{T,2} = [z_j for i in φ, z_j in z]

    Mesh{T}(X, Y, Z)
end

@recipe function f(es::EllipticalSurface; n = 40)
    colorbar := false
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :surface
            m = mesh(es, n = n)
            @series begin
                label --> "Elliptical Surface"
                occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend))) ? polymesh(m) : m
            end
    else
        ls = lines(es)
        linecolor --> :black
        @series begin
            label --> "Elliptical Surface"
            n := n
            ls[1]
        end
        if length(ls) > 1 
            for i in 2:length(ls)
                @series begin 
                    label := nothing
                    n := n
                    ls[i]
                end
            end
        end
    end
    if !haskey(plotattributes, :show_normal) || plotattributes[:show_normal]
        @series begin
            label := nothing
            seriestype := :vector
            _plt_get_start_point_for_normal(es), normalize(Plane(es).normal) * sqrt(_plt_area(es)) / 20 
        end
    end
end

_plt_area(es::CircularArea{T}) where {T} = π*es.r^2
_plt_area(es::PartialCircularArea{T}) where {T} = π*es.r^2 / ((es.φ[2]-es.φ[1])/2π)
_plt_area(es::Annulus{T}) where {T} = π*(es.r[2]^2 - es.r[1]^2)
_plt_area(es::PartialAnnulus{T}) where {T} = π*(es.r[2]^2 - es.r[1]^2) / ((es.φ[2]-es.φ[1])/2π)

_plt_get_start_point_for_normal(es::CircularArea{T}) where {T} = es.origin
function _plt_get_start_point_for_normal(es::PartialCircularArea{T}) where {T}
    cyl = CylindricalPoint{T}(es.r/2, (es.φ[2]+es.φ[1])/2, zero(T))
    _transform_into_global_coordinate_system(CartesianPoint(cyl), es)
end
function _plt_get_start_point_for_normal(es::Annulus{T}) where {T}
    cyl = CylindricalPoint{T}((es.r[2]+es.r[1])/2, zero(T), zero(T))
    _transform_into_global_coordinate_system(CartesianPoint(cyl), es)
end
function _plt_get_start_point_for_normal(es::PartialAnnulus{T}) where {T}
    cyl = CylindricalPoint{T}((es.r[2]+es.r[1])/2, (es.φ[2]+es.φ[1])/2, zero(T))
    _transform_into_global_coordinate_system(CartesianPoint(cyl), es)
end
