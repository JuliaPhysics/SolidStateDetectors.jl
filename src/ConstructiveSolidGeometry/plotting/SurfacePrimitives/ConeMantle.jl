@recipe function f(cm::ConeMantle, n = 40; subn = 10)
    e_top = get_top_ellipse(cm)
    e_bot = get_bot_ellipse(cm)
    e_top_edges = edges(e_top, n = n)
    e_bot_edges = edges(e_bot, n = n)
    linecolor --> :black
    @series begin
        label --> "Cone Mantle"
        e_top_edges
    end
    @series begin
        label := nothing
        e_bot_edges
    end
    @series begin
        label := nothing
        map(i -> Edge(e_top_edges[i].a, e_bot_edges[i].a), 1:subn:length(e_top_edges))
    end
end

    # RotZ(π) * -cm.rotation such that the normal vector points inside the cone (Convention)


# function get_plot_points(c::ConeMantle{T}; n = 30) where {T <: AbstractFloat}

#     plot_points = Vector{CartesianPoint{T}}[]

#     rbot::T, rtop::T = get_r_limits(c)
#     φMin::T, φMax::T, φ_is_full_2π::Bool = get_φ_limits(c)
#     zMin::T, zMax::T = get_z_limits(c)

#     φrange = range(φMin, φMax, length = n + 1)

#     #top and bottom circles
#     rbot != 0 ? push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rbot * cos(φ), rbot * sin(φ), zMin) for φ in φrange])) : nothing
#     rtop != 0 ? push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rtop * cos(φ), rtop * sin(φ), zMax) for φ in φrange])) : nothing

#     #side line(s)
#     for φ in (φ_is_full_2π ? T(0) : [φMin, φMax])
#         if rbot != 0 || rtop != 0
#             push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rbot * cos(φ), rbot * sin(φ), zMin), CartesianPoint{T}(rtop * cos(φ), rtop * sin(φ), zMax)]))
#         end
#     end
#     plot_points
# end

# function mesh(c::ConeMantle{T}; n = 30) where {T <: AbstractFloat}

#     rbot::T, rtop::T = get_r_limits(c)
#     φMin::T, φMax::T, _ = get_φ_limits(c)
#     zMin::T, zMax::T = get_z_limits(c)
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
