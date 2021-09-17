function mesh(em::EllipsoidMantle{T}; n = 30) where {T}
    rx, ry, rz = get_radii(em) 
    φMin::T, φMax::T = get_φ_limits(em)
    θMin::T, θMax::T = get_θ_limits(em)

    fφ = (φMax - φMin)/(2π)
    nφ = Int(ceil(n*fφ))

    fθ = (θMax - θMin)/(2π)
    nθ = Int(ceil(n*fθ))
    
    φ = range(φMin, φMax, length = nφ + 1)
    θ = range(θMin, θMax, length = nθ + 1)

    X = [rx*cos(θ_j)*cos(φ_i) for φ_i in φ, θ_j in θ]
    Y = [ry*cos(θ_j)*sin(φ_i) for φ_i in φ, θ_j in θ]
    Z = [rz*sin(θ_j) for i in φ, θ_j in θ]
    
    em.rotation*Mesh{T}(X,Y,Z) + em.origin
end

@recipe function f(em::EllipsoidMantle, n = 40; subn = 10)
    colorbar := false
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :surface
            m = mesh(em, n = n)
            @series begin
                label --> "Ellipsoid Mantle"
                occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend))) ? polymesh(m) : m
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
    if (!haskey(plotattributes, :show_normal) || plotattributes[:show_normal]) &&
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
