function mesh(em::EllipsoidMantle{T}; n = 30)::Mesh{T} where {T}
    rx, ry, rz = get_radii(em) 
    φMin::T, φMax::T = get_φ_limits(em)
    θMin::T, θMax::T = get_θ_limits(em)

    fφ = (φMax - φMin)/(2π)
    nφ = Int(ceil(n*fφ))

    fθ = (θMax - θMin)/(π)
    nθ = Int(ceil(n*fθ))
    
    θ = range(θMin, θMax, length = nθ + 1)
    sθ = sin.(θ)
    cθ = cos.(θ)
    φ = range(φMin, φMax, length = nφ + 1)
    sφ = sin.(φ)
    cφ = cos.(φ)

    x = [rx*cθ*cφ for cθ in cθ for cφ in cφ]
    y = [ry*cθ*sφ for cθ in cθ for sφ in sφ]
    z = [rz*sθ for sθ in sθ for i in φ]
    connections = [[i+(nφ+1)*j,i+1+(nφ+1)*j,i+1+(nφ+1)*(j+1),i+(nφ+1)*(j+1)] for j in 0:nθ-1 for i in 1:nφ]
    
    em.rotation*Mesh{T}(x,y,z,connections) + em.origin
end

@recipe function f(em::EllipsoidMantle, n = 40; subn = 10)
    seriestype --> :mesh3d
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :mesh3d
        @series begin
            label --> "Ellipsoid Mantle"
            mesh(em, n = n)
        end
    else
        ls = lines(em)
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
            em.φ === nothing && em.θ === nothing
        @series begin
            label := nothing
            seriestype := :vector
            pts = _plt_points_for_normals(em)
            ns = broadcast(p -> normal(em, p) / 5, pts)
            [(pts[i], ns[i]) for i in eachindex(pts)]
        end
    end
end

function _plt_points_for_normals(em::EllipsoidMantle{T,NTuple{3,T}}) where {T}
    pts = [ CartesianPoint{T}( em.r[1], zero(T), zero(T)),
            CartesianPoint{T}(-em.r[1], zero(T), zero(T)),
            CartesianPoint{T}(zero(T),  em.r[2], zero(T)),
            CartesianPoint{T}(zero(T), -em.r[2], zero(T)),
            CartesianPoint{T}(zero(T), zero(T),  em.r[3]),
            CartesianPoint{T}(zero(T), zero(T), -em.r[3]) ]
    _transform_into_global_coordinate_system(pts, em)
end
function _plt_points_for_normals(em::EllipsoidMantle{T,T}) where {T}
    pts = [ CartesianPoint{T}( em.r, zero(T), zero(T)),
            CartesianPoint{T}(-em.r, zero(T), zero(T)),
            CartesianPoint{T}(zero(T),  em.r, zero(T)),
            CartesianPoint{T}(zero(T), -em.r, zero(T)),
            CartesianPoint{T}(zero(T), zero(T),  em.r),
            CartesianPoint{T}(zero(T), zero(T), -em.r) ]
    _transform_into_global_coordinate_system(pts, em)
end
