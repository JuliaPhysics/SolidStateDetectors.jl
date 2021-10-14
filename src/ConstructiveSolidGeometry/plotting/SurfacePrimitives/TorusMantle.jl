@recipe function f(tm::TorusMantle; n_arc = 40, subn = 10)
    seriestype --> :mesh3d
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :mesh3d
        @series begin
            label --> "Torus Mantle"
            mesh(tm, n_arc = n_arc)
        end
    else
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
    end
    if (haskey(plotattributes, :show_normal) && plotattributes[:show_normal]) &&
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

function _plt_points_for_normals(tm::TorusMantle{T}) where {T}
    pts = [ CartesianPoint{T}( tm.r_torus+tm.r_tube, zero(T), zero(T)),
            CartesianPoint{T}( tm.r_torus-tm.r_tube, zero(T), zero(T)),
            CartesianPoint{T}( tm.r_torus, zero(T), tm.r_tube ),
            CartesianPoint{T}( tm.r_torus, zero(T),-tm.r_tube ),
            CartesianPoint{T}(-tm.r_torus+tm.r_tube, zero(T), zero(T)),
            CartesianPoint{T}(-tm.r_torus-tm.r_tube, zero(T), zero(T)),
            CartesianPoint{T}(-tm.r_torus, zero(T), tm.r_tube ),
            CartesianPoint{T}(-tm.r_torus, zero(T),-tm.r_tube ),
            CartesianPoint{T}( zero(T), tm.r_torus+tm.r_tube, zero(T)),
            CartesianPoint{T}( zero(T), tm.r_torus-tm.r_tube, zero(T)),
            CartesianPoint{T}( zero(T), tm.r_torus, tm.r_tube ),
            CartesianPoint{T}( zero(T), tm.r_torus,-tm.r_tube ),
            CartesianPoint{T}( zero(T),-tm.r_torus+tm.r_tube, zero(T)),
            CartesianPoint{T}( zero(T),-tm.r_torus-tm.r_tube, zero(T)),
            CartesianPoint{T}( zero(T),-tm.r_torus, tm.r_tube ),
            CartesianPoint{T}( zero(T),-tm.r_torus,-tm.r_tube ) ]
    _transform_into_global_coordinate_system(pts, tm)
end
