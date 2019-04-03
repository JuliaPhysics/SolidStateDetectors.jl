# function in(pt::AnyPoint{T}, vg::Vector{AbstractGeometry{T}})::Bool where {T <: AbstractFloat}
#     is_inside::Bool = false
#     for g in vg
#         if in(pt, g)
#             is_inside = true
#             break
#         end
#     end
#     return is_inside
# end
