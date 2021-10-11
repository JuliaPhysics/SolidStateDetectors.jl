function mesh(tm::TorusMantle{T}; n = 30)::Mesh{T} where {T}
    φMin::T, φMax::T = get_φ_limits(tm)
    θMin::T, θMax::T = get_θ_limits(tm)

    fφ = (φMax - φMin)/(2π)
    nφ = Int(ceil(n*fφ))

    fθ = (θMax - θMin)/(2π)
    nθ = Int(ceil(n*fθ))
    
    θ = range(θMin, θMax, length = nθ + 1)
    sθ = sin.(θ)
    cθ = cos.(θ)
    φ = range(φMin, φMax, length = nφ + 1)
    sφ = sin.(φ)
    cφ = cos.(φ)
    
    x = [(tm.r_torus + tm.r_tube*cθ)*cφ for cθ in cθ for cφ in cφ]
    y = [(tm.r_torus + tm.r_tube*cθ)*sφ for cθ in cθ for sφ in sφ]
    z = [tm.r_tube*sθ for sθ in sθ for i in φ]
    connections = [[i+(nφ+1)*j,i+1+(nφ+1)*j,i+1+(nφ+1)*(j+1),i+(nφ+1)*(j+1)] for j in 0:nθ-1 for i in 1:nφ]
    
    tm.rotation*Mesh{T}(x,y,z,connections) + tm.origin
end

@recipe function f(tm::TorusMantle, n = 40; subn = 10)
    seriestype --> :mesh3d
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :mesh3d
        @series begin
            label --> "Torus Mantle"
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
