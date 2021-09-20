function mesh(cm::ConeMantle{T}; n = 30)::Mesh{T} where {T}
    φMin, φMax = get_φ_limits(cm)
    
    f = (φMax - φMin)/(2π)
    n = Int(ceil(n*f))
    φ = range(φMin, φMax, length = n+1)
    z = range(-cm.hZ, cm.hZ, length = 2)

    X::Array{T,2} = [radius_at_z(cm,z_j)*cos(φ_i) for φ_i in φ, z_j in z]
    Y::Array{T,2} = [radius_at_z(cm,z_j)*sin(φ_i) for φ_i in φ, z_j in z]
    Z::Array{T,2} = [z_j for i in φ, z_j in z]
    
    cm.rotation*Mesh{T}(X,Y,Z) + cm.origin
end

@recipe function f(cm::ConeMantle, n = 40; subn = 10)
    colorbar := false
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :surface
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
