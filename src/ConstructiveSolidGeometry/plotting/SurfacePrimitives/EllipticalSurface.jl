@recipe function f(es::EllipticalSurface; n = 40)
    ls = lines(es)
    linecolor --> :black
    @series begin
        label --> "Elliptical Surface"
        n := n
        ls[1]
    end
    if !haskey(plotattributes, :show_normal) || plotattributes[:show_normal]
        @series begin
            label := nothing
            seriestype := :vector
            _plt_get_start_point_for_normal(es), normalize(Plane(es).normal) * sqrt(_plt_area(es)) / 20 
        end
    end
    if length(ls) > 1 
        for i in 2:length(ls)
            @series begin 
                label := nothing
                n := n
                ls[i]
            end
        end
    end
end

_plt_area(es::CircularArea{T}) where {T} = π*es.r^2
_plt_area(es::PartialCircularArea{T}) where {T} = π*es.r^2 / ((es.φ[2]-es.φ[1])/2π)
_plt_area(es::Annulus{T}) where {T} = π*(es.r[2]^2 - es.r[1]^2)
_plt_area(es::PartialAnnulus{T}) where {T} = π*(es.r[2]^2 - es.r[1]^2) / ((es.φ[2]-es.φ[1])/2π)

_plt_get_start_point_for_normal(es::CircularArea{T}) where {T} = es.origin
function _plt_get_start_point_for_normal(es::PartialCircularArea{T}) where {T}
    cyl = CylindricalPoint{T}(es.r/2, (es.φ[2]+es.φ[1])/2, zero(T))
    _transform_into_global_coordinate_system(CartesianPoint(cyl), es)
end
function _plt_get_start_point_for_normal(es::Annulus{T}) where {T}
    cyl = CylindricalPoint{T}((es.r[2]+es.r[1])/2, zero(T), zero(T))
    _transform_into_global_coordinate_system(CartesianPoint(cyl), es)
end
function _plt_get_start_point_for_normal(es::PartialAnnulus{T}) where {T}
    cyl = CylindricalPoint{T}((es.r[2]+es.r[1])/2, (es.φ[2]+es.φ[1])/2, zero(T))
    _transform_into_global_coordinate_system(CartesianPoint(cyl), es)
end

# function get_plot_points(a::CylindricalAnnulus{T}; n = 30) where {T <: AbstractFloat}

#     plot_points = Vector{CartesianPoint{T}}[]

#     rMin::T, rMax::T = get_r_limits(a)
#     φMin::T, φMax::T, φ_is_full_2π::Bool = get_φ_limits(a)
#     φrange = range(φMin, φMax, length = n + 1)

#     #circle(s)
#     for r in [rMin, rMax]
#         if r == 0 continue end
#         push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * cos(φ), r * sin(φ), a.z) for φ in φrange]))
#     end

#     #for incomplete φ: lines of cross-sections
#     if !φ_is_full_2π
#         for φ in [φMin, φMax]
#             push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rMin * cos(φ), rMin * sin(φ), a.z), CartesianPoint{T}(rMax * cos(φ), rMax * sin(φ), a.z)]))
#         end
#     end
#     plot_points
# end

# function mesh(a::CylindricalAnnulus{T}; n = 30) where {T <: AbstractFloat}

#     rMin::T, rMax::T = get_r_limits(a)
#     φMin::T, φMax::T, _ = get_φ_limits(a)
#     f = (φMax - φMin)/(2π)
#     n = Int(ceil(n*f))
#     φrange = range(φMin, φMax, length = n+1)
#     sφrange = sin.(φrange)
#     cφrange = cos.(φrange)
#     r = range(rMin, rMax, length = 2)
#     z = fill(a.z, length(r))

#     X::Array{T,2} = [r_j*cφ for cφ in cφrange, r_j in r]
#     Y::Array{T,2} = [r_j*sφ for sφ in sφrange, r_j in r]
#     Z::Array{T,2} = [z_j for i in 1:n+1, z_j in z]

#     Mesh(X, Y, Z)
# end
