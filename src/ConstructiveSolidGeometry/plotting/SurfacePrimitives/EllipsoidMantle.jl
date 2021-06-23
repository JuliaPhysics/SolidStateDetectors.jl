@recipe function f(em::EllipsoidMantle, n = 40; subn = 10)
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
    # if !haskey(plotattributes, :show_normal) || plotattributes[:show_normal]
    #     @series begin
    #         label := nothing
    #         seriestype := :vector
    #         nφ = cm.φ == nothing ? 0 : (cm.φ[2] + cm.φ[1])/2 
    #         T = typeof(cm.hZ)
    #         npt_obj = CartesianPoint(CylindricalPoint{T}(radius_at_z(cm, zero(T)), nφ, zero(T)))
    #         npt = _transform_into_global_coordinate_system(npt_obj, cm)
    #         npt, normalize(normal(cm, npt))/10
    #     end
    # end
end

# function get_plot_points(s::SphereMantle{T}; n = 30) where {T <: AbstractFloat}
#     plot_points = Vector{CartesianPoint{T}}[]

#     φrange = range(0, 2π, length = n)

#     for φ in φrange
#         push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(s.r * sin(θ) * cos(φ), s.r * sin(θ) * sin(φ), s.r * cos(θ)) for θ in φrange]))
#     end
#     push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(s.r * cos(φ), s.r * sin(φ), 0) for φ in φrange]))

#     plot_points
# end

# function mesh(s::SphereMantle{T}; n = 30) where {T <: AbstractFloat}

#     θrange = range(-π/2, π/2, length = n)
#     sθrange = sin.(θrange)
#     cθrange = cos.(θrange)
#     φrange = range(0, 2π, length = n)
#     sφrange = sin.(φrange)
#     cφrange = cos.(φrange)

#     X = [s.r*cθ*cφ for cφ in cφrange, cθ in cθrange]
#     Y = [s.r*cθ*sφ for sφ in sφrange, cθ in cθrange]
#     Z = [s.r*sθ for i in 1:n, sθ in sθrange]
#     Mesh(X, Y, Z)
# end
