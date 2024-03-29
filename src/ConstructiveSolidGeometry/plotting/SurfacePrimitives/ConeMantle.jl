@recipe function f(cm::ConeMantle; n_arc = 40, subn = 10, n_vert_lines = 2)
    seriestype --> :mesh3d
    @series begin
        label --> "Cone Mantle"
        if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :mesh3d   
            mesh(cm, n_arc)
        else
            mesh(cm, n_arc, n_vert_lines)
        end    
    end
    if haskey(plotattributes, :show_normal) && plotattributes[:show_normal]
        @series begin
            label := nothing
            seriestype := :vector
            nφ = isnothing(cm.φ) ? 0 : cm.φ/2 
            T = typeof(cm.hZ)
            npt_obj = CartesianPoint(CylindricalPoint{T}(radius_at_z(cm, zero(T)), nφ, zero(T)))
            npt = _transform_into_global_coordinate_system(npt_obj, cm)
            npt, normal(cm, npt)/3
        end
    end
end
