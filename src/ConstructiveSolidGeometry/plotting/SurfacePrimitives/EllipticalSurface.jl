@recipe function f(es::EllipticalSurface; n_arc = 40, n_vert_lines = 2)
    seriestype --> :mesh3d
    @series begin
        label --> "Elliptical Surface"
        mesh(es, n_arc = n_arc, n_vert_lines = n_vert_lines)
    end
    if haskey(plotattributes, :show_normal) && plotattributes[:show_normal]
        @series begin
            label := nothing
            seriestype := :vector
            _plt_get_start_point_for_normal(es), normalize(Plane(es).normal) * sqrt(_plt_area(es)) / 20 
        end
    end
end

_plt_area(es::CircularArea{T}) where {T} = π*es.r^2
_plt_area(es::PartialCircularArea{T}) where {T} = π*es.r^2 / (es.φ/2π)
_plt_area(es::Annulus{T}) where {T} = π*(es.r[2]^2 - es.r[1]^2)
_plt_area(es::PartialAnnulus{T}) where {T} = π*(es.r[2]^2 - es.r[1]^2) / (es.φ/2π)

_plt_get_start_point_for_normal(es::CircularArea{T}) where {T} = es.origin
function _plt_get_start_point_for_normal(es::PartialCircularArea{T}) where {T}
    cyl = CylindricalPoint{T}(es.r/2, es.φ/2, zero(T))
    _transform_into_global_coordinate_system(CartesianPoint(cyl), es)
end
function _plt_get_start_point_for_normal(es::Annulus{T}) where {T}
    cyl = CylindricalPoint{T}((es.r[2]+es.r[1])/2, zero(T), zero(T))
    _transform_into_global_coordinate_system(CartesianPoint(cyl), es)
end
function _plt_get_start_point_for_normal(es::PartialAnnulus{T}) where {T}
    cyl = CylindricalPoint{T}((es.r[2]+es.r[1])/2, es.φ/2, zero(T))
    _transform_into_global_coordinate_system(CartesianPoint(cyl), es)
end
