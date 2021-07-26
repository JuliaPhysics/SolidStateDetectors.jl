@recipe function f(cm::ConeMantle, n = 40; subn = 10)
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

    # RotZ(π) * -cm.rotation such that the normal vector points inside the cone (Convention)


# function mesh(c::ConeMantle{T}; n = 30) where {T <: AbstractFloat}

#     rbot::T, rtop::T = get_r_limits(c)
#     φMin::T, φMax::T, _ = get_φ_limits(c)
#     zMin::T, zMax::T = _linear_endpoints(c.hZ)
#     f = (φMax - φMin)/(2π)
#     n = Int(ceil(n*f))
#     m = (rtop-rbot)/(zMax-zMin)
#     φrange = range(φMin, φMax, length = n+1)
#     scφrange = sincos.(φrange)
#     z = (zMin, zMax)

#     X::Array{T,2} = [(m*(z_i-zMin)+rbot)*cφ for (_,cφ) in scφrange, z_i in z]
#     Y::Array{T,2} = [(m*(z_i-zMin)+rbot)*sφ for (sφ,_) in scφrange, z_i in z]
#     Z::Array{T,2} = [j for i in φrange, j in z]

#     Mesh(X, Y, Z)
# end
