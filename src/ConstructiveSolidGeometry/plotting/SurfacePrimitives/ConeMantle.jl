function mesh(cm::ConeMantle{T}; n = 30)::Mesh{T} where {T}
    φMin, φMax = get_φ_limits(cm)
    
    f = (φMax - φMin)/(2π)
    n = Int(ceil(n*f))
    φ = range(φMin, φMax, length = n+1)
    rtop = radius_at_z(cm,cm.hZ)
    rbot = radius_at_z(cm,-cm.hZ)
    x = append!(rbot*cos.(φ), rtop*cos.(φ))
    y = append!(rbot*sin.(φ), rtop*sin.(φ))
    z = append!(rbot*ones(n+1), rtop*ones(n+1))
    connections = [[i,i+1,i+n+2,i+n+1] for i in 1:n]

    cm.rotation*Mesh{T}(x,y,z,connections) + cm.origin
end

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
