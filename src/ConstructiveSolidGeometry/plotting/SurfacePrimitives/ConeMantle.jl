@recipe function f(cm::ConeMantle, n = 40; subn = 10)
    seriestype --> :mesh3d
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :mesh3d
        @series begin
            label --> "Cone Mantle"
            mesh(cm, n = n)
        end
    else
        ls = lines(cm)
        linecolor --> :black
        @series begin
            label --> "Cone Mantle"
            ls[1]
        end
        for i in 2:length(ls)
            @series begin
                label := nothing
                ls[i]
            end
        end
    end
    if haskey(plotattributes, :show_normal) && plotattributes[:show_normal]
        @series begin
            label := nothing
            seriestype := :vector
            nφ = cm.φ == nothing ? 0 : (cm.φ[2] + cm.φ[1])/2 
            T = typeof(cm.hZ)
            npt_obj = CartesianPoint(CylindricalPoint{T}(radius_at_z(cm, zero(T)), nφ, zero(T)))
            npt = _transform_into_global_coordinate_system(npt_obj, cm)
            npt, normal(cm, npt)/3
        end
    end
end
