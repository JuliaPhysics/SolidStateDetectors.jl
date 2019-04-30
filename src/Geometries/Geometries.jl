# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

include("Units.jl")

include("Points.jl")
include("Vectors.jl")

abstract type AbstractGeometry{T, N} end
abstract type AbstractVolumePrimitive{T, N} <: AbstractGeometry{T, N} end
abstract type AbstractSurfacePrimitive{T, N} <: AbstractGeometry{T, N} end
abstract type AbstractSet{T, N} <: AbstractGeometry{T, N} end

function Geometry(T, geom_dict::Union{Dict{String,Any},Dict{Any,Any}}, inputunit_dict::Dict{String, Unitful.Units} )
    Geometry(T, Val{Symbol(geom_dict["type"])}(), Dict{Any,Any}(geom_dict), inputunit_dict)
end
# function Geometry(T, geom_sub_dict::Dict, geom_unit)
#     haskey(geom_sub_dict, "positive") ? geometry_positive = Geometries(T, geom_sub_dict["positive"], geom_unit) : geometry_positive = AbstractGeometry[]
#     haskey(geom_sub_dict, "negative") ? geometry_negative = Geometries(T, geom_sub_dict["negative"], geom_unit) : geometry_negative = AbstractGeometry[]
#     geometry = geometry_positive[1]
#     for g in sort!(vcat(geometry_positive,geometry_negative))
#         g in geometry_negative ? geometry -= g : nothing
#         g in geometry_positive ? geometry += g : nothing
#     end
#     return geometry, geometry_positive, geometry_negative
# end

function Geometry(T, geom_sub_dict::Vector, geom_unit)
    geometry_positive = Geometries(T, geom_sub_dict, geom_unit)
    geometry_negative = AbstractGeometry[]
    geometry = geometry_positive[1]
    for g in sort!(vcat(geometry_positive,geometry_negative))
        g in geometry_negative ? geometry -= g : nothing
        g in geometry_positive ? geometry += g : nothing
    end
    return geometry, geometry_positive, geometry_negative
end

function in(pt::AbstractCoordinatePoint{T}, vg::Vector{AbstractGeometry{T}})::Bool where {T <: SSDFloat}
    inside::Bool = false
    for g in vg
        if in(pt, g)
            inside = true
            break
        end
    end
    return inside
end

include("VolumePrimitives/VolumePrimitives.jl")
include("SurfacePrimitives/SurfacePrimitives.jl")
include("LinePrimitives/LinePrimitives.jl")

include("Sets.jl")


include("IO.jl")
function get_decomposed_volumes(vol::AbstractGeometry{T})::Tuple{Vector{<:AbstractGeometry},Vector{<:AbstractGeometry}} where T
    positive_volumes = AbstractVolumePrimitive[]
    negative_volumes = AbstractVolumePrimitive[]
    translate = CartesianVector{T}(0.0,0.0,0.0)
    decompose_volume(vol, positive_volumes, negative_volumes, translate,  1)
    positive_volumes, negative_volumes
end

function decompose_volume(vol::Difference, pos, neg, translate, flag = 1)::Nothing
    ismissing(vol.translate) ? nothing : translate += vol.translate
    decompose_volume(vol.a, pos, neg, translate, flag)
    decompose_volume(vol.b, pos, neg, translate, flag *= -1)
    nothing
end

function decompose_volume(vol::GeometryUnion, pos, neg, translate, flag = 1)::Nothing
    ismissing(vol.translate) ? nothing : translate += vol.translate
    decompose_volume(vol.a ,pos, neg, translate, flag)
    decompose_volume(vol.b, pos, neg, translate, flag)
    nothing
end

function decompose_volume(vol::Intersection, pos, neg, translate, flag = 1)::Nothing
    ismissing(vol.translate) ? nothing : translate += vol.translate
    decompose_volume(vol.a ,pos, neg, translate, flag)
    decompose_volume(vol.b, pos, neg, translate, flag)
    nothing
end

function decompose_volume(vol::AbstractVolumePrimitive, pos, neg, translate, flag = 1)::Nothing
    # ismissing(vol.translate) ? vol.translate = translate : vol.translate += translate
    flag == 1 ? push!(pos, vol+translate) : push!(neg, vol+translate)
    nothing
end

function sample(geometry::AbstractGeometry{T}, stepsize::Vector{T}) where T
    return filter(x-> x in geometry, unique!(vcat([sample(g, stepsize) for g in get_decomposed_volumes(geometry)[1] ]...)))
end
