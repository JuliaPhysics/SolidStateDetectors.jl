@recipe function f(tm::TorusMantle, n = 40; subn = 10)
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
    if (!haskey(plotattributes, :show_normal) || plotattributes[:show_normal]) &&
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

function _plt_points_for_normals(em::TorusMantle{T}) where {T}
    pts = [ CartesianPoint{T}( em.r_torus+em.r_tube, zero(T), zero(T)),
            CartesianPoint{T}( em.r_torus-em.r_tube, zero(T), zero(T)),
            CartesianPoint{T}( em.r_torus, zero(T), em.r_tube ),
            CartesianPoint{T}( em.r_torus, zero(T),-em.r_tube ),
            CartesianPoint{T}(-em.r_torus+em.r_tube, zero(T), zero(T)),
            CartesianPoint{T}(-em.r_torus-em.r_tube, zero(T), zero(T)),
            CartesianPoint{T}(-em.r_torus, zero(T), em.r_tube ),
            CartesianPoint{T}(-em.r_torus, zero(T),-em.r_tube ),
            CartesianPoint{T}( zero(T), em.r_torus+em.r_tube, zero(T)),
            CartesianPoint{T}( zero(T), em.r_torus-em.r_tube, zero(T)),
            CartesianPoint{T}( zero(T), em.r_torus, em.r_tube ),
            CartesianPoint{T}( zero(T), em.r_torus,-em.r_tube ),
            CartesianPoint{T}( zero(T),-em.r_torus+em.r_tube, zero(T)),
            CartesianPoint{T}( zero(T),-em.r_torus-em.r_tube, zero(T)),
            CartesianPoint{T}( zero(T),-em.r_torus, em.r_tube ),
            CartesianPoint{T}( zero(T),-em.r_torus,-em.r_tube ) ]
    _transform_into_global_coordinate_system(pts, em)
end

# function get_plot_points(t::TorusMantle{T}; n = 30) where {T <: AbstractFloat}
#     plot_points = Vector{CartesianPoint{T}}[]
#     φMin::T, φMax::T, φ_is_full_2π::Bool = get_φ_limits(t)
#     θMin::T, θMax::T, θ_is_full_2π::Bool = get_θ_limits(t)
#     sθMin, cθMin = sincos(θMin)

#     r1 = T(t.r_torus + t.r_tube*cθMin)
#     z1 = T(t.r_tube*sθMin) + t.z
#     θ2 = θ_is_full_2π ? π : θMax
#     sθ2, cθ2 = sincos(θ2)
#     r2 = T(t.r_torus + t.r_tube*cθ2)
#     z2 = T(t.r_tube*sθ2) + t.z

#     append!(plot_points, get_plot_points(CylindricalAnnulus(T,r1..r1,t.φ,z1), n = n))
#     append!(plot_points, get_plot_points(CylindricalAnnulus(T,r2..r2,t.φ,z2), n = n))
#     for φ in (φ_is_full_2π ? [φMin] : [φMin, φMax])
#         append!(plot_points, get_plot_points(ToroidalAnnulus(T,t.r_torus,t.r_tube..t.r_tube,φ,t.θ,t.z), n = n))
#     end
#     plot_points
# end

# function mesh(t::TorusMantle{T}; n = 30) where {T <: AbstractFloat}
#     φMin::T, φMax::T, _ = get_φ_limits(t)
#     θMin::T, θMax::T, _ = get_θ_limits(t)

#     fφ = (φMax - φMin)/(2π)
#     nφ = Int(ceil(n*fφ))

#     fθ = (θMax - θMin)/(2π)
#     nθ = Int(ceil(n*fθ))

#     θrange = range(θMin, θMax, length = nθ + 1)
#     sθrange = sin.(θrange)
#     cθrange = cos.(θrange)
#     φrange = range(φMin, φMax, length = nφ + 1)
#     sφrange = sin.(φrange)
#     cφrange = cos.(φrange)

#     X = [(t.r_torus + t.r_tube*cθ)*cφ for cφ in cφrange, cθ in cθrange]
#     Y = [(t.r_torus + t.r_tube*cθ)*sφ for sφ in sφrange, cθ in cθrange]
#     Z = [t.r_tube*sθ + t.z for i in 1:nφ+1, sθ in sθrange]
#     Mesh(X, Y, Z)
# end
