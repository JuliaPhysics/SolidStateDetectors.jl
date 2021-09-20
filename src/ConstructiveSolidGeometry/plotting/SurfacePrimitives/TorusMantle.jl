function mesh(tm::TorusMantle{T}; n = 30)::Mesh{T} where {T}
    φMin::T, φMax::T = get_φ_limits(tm)
    θMin::T, θMax::T = get_θ_limits(tm)

    fφ = (φMax - φMin)/(2π)
    nφ = Int(ceil(n*fφ))

    fθ = (θMax - θMin)/(2π)
    nθ = Int(ceil(n*fθ))
    
    θrange = range(θMin, θMax, length = nθ + 1)
    sθrange = sin.(θrange)
    cθrange = cos.(θrange)
    φrange = range(φMin, φMax, length = nφ + 1)
    sφrange = sin.(φrange)
    cφrange = cos.(φrange)

    X::Array{T,2} = [(tm.r_torus + tm.r_tube*cθ)*cφ for cφ in cφrange, cθ in cθrange]
    Y::Array{T,2} = [(tm.r_torus + tm.r_tube*cθ)*sφ for sφ in sφrange, cθ in cθrange]
    Z::Array{T,2} = [tm.r_tube*sθ for i in φrange, sθ in sθrange]
    
    tm.rotation*Mesh{T}(X,Y,Z) + tm.origin
end

@recipe function f(tm::TorusMantle, n = 40; subn = 10)
    colorbar := false
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :surface
        @series begin
            label --> "Ellipsoid Mantle"
            mesh(tm, n = n)
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
