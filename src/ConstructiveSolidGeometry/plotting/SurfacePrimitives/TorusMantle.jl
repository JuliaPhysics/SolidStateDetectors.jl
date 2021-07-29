@recipe function f(tm::TorusMantle, n = 40; subn = 10)
    ls = lines(tm)
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
    if (!haskey(plotattributes, :show_normal) || plotattributes[:show_normal]) &&
            tm.φ === nothing && tm.θ === nothing
        @series begin
            label := nothing
            seriestype := :vector
            pts = _plt_points_for_normals(tm)
            ns = broadcast(p -> normal(tm, p) / 5, pts)
            [(pts[i], ns[i]) for i in eachindex(pts)]
        end
    end
end

function _plt_points_for_normals(em::TorusMantle{T}) where {T}
    pts = [ CartesianPoint{T}( em.r_torus+em.r_tube, zero(T), zero(T)),
            CartesianPoint{T}( em.r_torus-em.r_tube, zero(T), zero(T)),
            CartesianPoint{T}( em.r_torus, zero(T), em.r_tube ),
            CartesianPoint{T}( em.r_torus, zero(T),-em.r_tube ),
            CartesianPoint{T}(-em.r_torus+em.r_tube, zero(T), zero(T)),
            CartesianPoint{T}(-em.r_torus-em.r_tube, zero(T), zero(T)),
            CartesianPoint{T}(-em.r_torus, zero(T), em.r_tube ),
            CartesianPoint{T}(-em.r_torus, zero(T),-em.r_tube ),
            CartesianPoint{T}( zero(T), em.r_torus+em.r_tube, zero(T)),
            CartesianPoint{T}( zero(T), em.r_torus-em.r_tube, zero(T)),
            CartesianPoint{T}( zero(T), em.r_torus, em.r_tube ),
            CartesianPoint{T}( zero(T), em.r_torus,-em.r_tube ),
            CartesianPoint{T}( zero(T),-em.r_torus+em.r_tube, zero(T)),
            CartesianPoint{T}( zero(T),-em.r_torus-em.r_tube, zero(T)),
            CartesianPoint{T}( zero(T),-em.r_torus, em.r_tube ),
            CartesianPoint{T}( zero(T),-em.r_torus,-em.r_tube ) ]
    _transform_into_global_coordinate_system(pts, em)
end
