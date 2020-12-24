

# Old
# @recipe function f(l::AbstractLine{T, 3, :cartesian}) where {T}
#     x::Vector{T} = [l.org.x, l.org.x + l.dir.x]
#     y::Vector{T} = [l.org.y, l.org.y + l.dir.y]
#     z::Vector{T} = [l.org.z, l.org.z + l.dir.z]
#     x, y, z
# end

# @recipe function f(l::AbstractLine{T, 2, :cartesian}) where {T}
#     x::Vector{T} = [l.org.x, l.org.x + l.dir.x]
#     y::Vector{T} = [l.org.y, l.org.y + l.dir.y]
#     x, y
# end

# @recipe function f(pc::PartialCircle{T, 3, :cartesian}; n = 30) where {T}
#     phirange = range(pc.phiStart, pc.phiStop, length = round(Int, n + 1))
#     x::Vector{T} = pc.r .* cos.(phirange)
#     y::Vector{T} = pc.r .* sin.(phirange)
#     z::Vector{T} = map(phi -> 0.0, phirange)
#     points = map(p -> pc.Rotate*p + pc.Translate, CartesianPoint{T}.(x,y,z))
#     points
# end

# @recipe function f(ls::Array{<:AbstractLine{T, 3, :cartesian}, 1}) where {T}
#     seriescolor --> 1
#     for l in ls
#         @series begin
#             label := ""
#             l
#         end
#     end
#     @series begin
#         label --> ""
#         [], []
#     end
# end
