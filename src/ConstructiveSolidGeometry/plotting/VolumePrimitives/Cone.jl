
_get_plot_r(c::Cone{T, <:Union{T, AbstractInterval{T}}, <:Any, <:Any}) where {T} =
    (_radial_endpoints(c.r)..., _radial_endpoints(c.r)...)
_get_plot_r(c::Cone{T, <:Tuple, <:Any, <:Any}) where {T} = 
    (_radial_endpoints(c.r[1])..., _radial_endpoints(c.r[2])...)

_get_plot_φ(c::Cone{T, <:Any, Nothing, <:Any}) where {T} = (T(0), T(2π), true)
_get_plot_φ(c::Cone{T, <:Any, <:AbstractInterval, <:Any}) where {T} = (c.φ.left, c.φ.right, false)

_get_plot_z(c::Cone{T}) where {T} = _linear_endpoints(c.z)


function get_plot_points(c::Cone{T}; n = 30) where {T <: AbstractFloat}
    
    plot_points = Vector{CartesianPoint{T}}[]
    
    rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = _get_plot_r(c)
    φMin::T, φMax::T, φ_is_full_2π::Bool = _get_plot_φ(c)
    zMin::T, zMax::T = _get_plot_z(c)
    
    φrange = range(φMin, φMax, length = n)
    
    #bottom circle(s)
    for r in [rbotMin, rbotMax]
        if r == 0 continue end
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * cos(φ), r * sin(φ), zMin) for φ in φrange]))
    end
    
    #top circle(s)
    for r in [rtopMin, rtopMax]
        if r == 0 continue end
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * cos(φ), r * sin(φ), zMax) for φ in φrange]))
    end
    
    #side line(s)
    for φ in (φ_is_full_2π ? T(0) : [φMin, φMax])    
        if rbotMin != 0 || rtopMin != 0
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rbotMin * cos(φ), rbotMin * sin(φ), zMin), CartesianPoint{T}(rtopMin * cos(φ), rtopMin * sin(φ), zMax)]))      
        end
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rbotMax * cos(φ), rbotMax * sin(φ), zMin), CartesianPoint{T}(rtopMax * cos(φ), rtopMax * sin(φ), zMax)]))       
    end
    
    #for incomplete φ: lines of cross-sections
    if !φ_is_full_2π
        for φ in [φMin, φMax]
            push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rbotMin * cos(φ), rbotMin * sin(φ), zMin), CartesianPoint{T}(rbotMax * cos(φ), rbotMax * sin(φ), zMin)]))
            push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rtopMin * cos(φ), rtopMin * sin(φ), zMax), CartesianPoint{T}(rtopMax * cos(φ), rtopMax * sin(φ), zMax)]))
        end
    end
    plot_points
end


# Old
# @recipe function f(t::Tube{T}; n = 30, seriescolor = :green) where {T}
#     linewidth --> 2
#     n --> n
#     @series begin
#         seriescolor --> seriescolor
#         label --> "Tube"
#         []
#     end
#     label := ""
#     seriescolor := seriescolor
#     LineSegments(t)
# end
# @recipe function f(c::Cone{T}; n = 30, seriescolor = :orange) where {T}
#     linewidth --> 2
#     n --> n
#     @series begin
#         seriescolor --> seriescolor
#         label --> "Cone"
#         []
#     end
#     seriescolor := seriescolor
#     label := ""
#     LineSegments(c)
# end