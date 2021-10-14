@recipe function f(em::EllipsoidMantle, n = 40; subn = 10)
    seriestype --> :mesh3d
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :mesh3d
        @series begin
            label --> "Ellipsoid Mantle"
            mesh(em, n = n)
        end
    else
        ls = lines(em)
        linecolor --> :black
        @series begin
            label --> "Ellipsoid Mantle"
            ls[1]
        end
        for i in 2:length(ls)
            @series begin
                label := nothing
                ls[i]
            end
        end
    end
    if (haskey(plotattributes, :show_normal) && plotattributes[:show_normal]) &&
            em.φ === nothing && em.θ === nothing
        @series begin
            label := nothing
            seriestype := :vector
            pts = _plt_points_for_normals(em)
            ns = broadcast(p -> normal(em, p) / 5, pts)
            [(pts[i], ns[i]) for i in eachindex(pts)]
        end
    end
end

function _plt_points_for_normals(em::EllipsoidMantle{T,NTuple{3,T}}) where {T}
    pts = [ CartesianPoint{T}( em.r[1], zero(T), zero(T)),
            CartesianPoint{T}(-em.r[1], zero(T), zero(T)),
            CartesianPoint{T}(zero(T),  em.r[2], zero(T)),
            CartesianPoint{T}(zero(T), -em.r[2], zero(T)),
            CartesianPoint{T}(zero(T), zero(T),  em.r[3]),
            CartesianPoint{T}(zero(T), zero(T), -em.r[3]) ]
    _transform_into_global_coordinate_system(pts, em)
end
function _plt_points_for_normals(em::EllipsoidMantle{T,T}) where {T}
    pts = [ CartesianPoint{T}( em.r, zero(T), zero(T)),
            CartesianPoint{T}(-em.r, zero(T), zero(T)),
            CartesianPoint{T}(zero(T),  em.r, zero(T)),
            CartesianPoint{T}(zero(T), -em.r, zero(T)),
            CartesianPoint{T}(zero(T), zero(T),  em.r),
            CartesianPoint{T}(zero(T), zero(T), -em.r) ]
    _transform_into_global_coordinate_system(pts, em)
end
