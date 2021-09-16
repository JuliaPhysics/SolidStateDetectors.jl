function mesh(cm::ConeMantle{T}; n = 30)::Mesh{T} where {T}
    rbot, rtop = get_r_limits(cm)
    φMin, φMax = get_φ_limits(cm)
    
    f = (φMax - φMin)/(2π)
    n = Int(ceil(n*f))
    φ = range(φMin, φMax, length = n+1)
    z = range(-cm.hZ, cm.hZ, length = 2)

    X::Array{T,2} = [radius_at_z(cm,z[j])*cos(φ_i) for φ_i in φ, j in 1:length(z)]
    Y::Array{T,2} = [radius_at_z(cm,z[j])*sin(φ_i) for φ_i in φ, j in 1:length(z)]
    Z::Array{T,2} = [j for i in 1:length(φ), j in z]
    if cm.rotation != one(SMatrix{3, 3, T, 9})
        cm.rotation*Mesh{T}(X,Y,Z) + cm.origin
    else
        Mesh{T}(X,Y,Z) + cm.origin
    end
end

@recipe function f(cm::ConeMantle, n = 40; subn = 10)
    colorbar := false
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :surface
            m = mesh(cm, n = n)
            if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
                @series begin
                    label --> "Cone Mantle"
                    polymesh(m)
                end
            else
                @series begin
                    label --> "Cone Mantle"
                    m
                end
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
    if !haskey(plotattributes, :show_normal) || plotattributes[:show_normal]
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
