function get_decomposed_volumes(vol::AbstractGeometry{T})::Tuple{Vector{<:AbstractGeometry},Vector{<:AbstractGeometry}} where {T}
    positive_volumes = AbstractGeometry[]
    negative_volumes = AbstractGeometry[]
    transformation = CSGTransformation[]
    decompose!(positive_volumes, negative_volumes, vol, transformation, Val{:PositiveVolume}, AbstractVolumePrimitive)
    positive_volumes, negative_volumes
end

function get_decomposed_surfaces(vol::AbstractGeometry{T})::Vector{<:AbstractGeometry} where {T}
    positive_surfaces = AbstractGeometry[]
    negative_surfaces = AbstractGeometry[]
    transformation = CSGTransformation[]
    decompose!(positive_surfaces, negative_surfaces, vol, transformation, Val{:PositiveVolume}, AbstractSurfacePrimitive)
    positive_surfaces
end

@inline invert_volume_type(::Type{Val{:PositiveVolume}}) = Val{:NegativeVolume}
@inline invert_volume_type(::Type{Val{:NegativeVolume}}) = Val{:PositiveVolume}

function decompose!(pos, neg, vol::CSGDifference, transformation, ::Type{VT}, ::Type{P})::Nothing where {VT <:Union{Val{:PositiveVolume}, Val{:NegativeVolume}}, P <: Union{AbstractVolumePrimitive, AbstractSurfacePrimitive}}
    decompose!(pos, neg, vol.a, transformation, VT, P)
    decompose!(pos, neg, vol.b, transformation, invert_volume_type(VT), P)
    nothing
end

function decompose!(pos, neg, vol::Union{CSGUnion, CSGIntersection}, transformation, ::Type{VT}, ::Type{P})::Nothing where {VT <:Union{Val{:PositiveVolume}, Val{:NegativeVolume}}, P <: Union{AbstractVolumePrimitive, AbstractSurfacePrimitive}}
    decompose!(pos, neg, vol.a, transformation, VT, P)
    decompose!(pos, neg, vol.b, transformation, VT, P)
    nothing
end

function decompose!(pos, neg, vol::TranslatedGeometry, transformation, ::Type{VT}, ::Type{P})::Nothing where {VT <:Union{Val{:PositiveVolume}, Val{:NegativeVolume}}, P <: Union{AbstractVolumePrimitive, AbstractSurfacePrimitive}}
    push!(transformation, vol.t)
    decompose!(pos, neg, vol.p, transformation, VT, P)
    pop!(transformation)
    nothing
end

function decompose!(pos, neg, vol::RotatedGeometry, transformation, ::Type{VT}, ::Type{P})::Nothing where {VT <:Union{Val{:PositiveVolume}, Val{:NegativeVolume}}, P <: Union{AbstractVolumePrimitive, AbstractSurfacePrimitive}}
    push!(transformation, inv(vol.inv_r))
    decompose!(pos, neg, vol.p, transformation, VT, P)
    pop!(transformation)
    nothing
end

function decompose!(pos, neg, vol::ScaledGeometry, transformation, ::Type{VT}, ::Type{P})::Nothing where {VT <:Union{Val{:PositiveVolume}, Val{:NegativeVolume}}, P <: Union{AbstractVolumePrimitive, AbstractSurfacePrimitive}}
    push!(transformation, inv.(vol.inv_s))
    decompose!(pos, neg, vol.p, transformation, VT, P)
    pop!(transformation)
    nothing
end

decompose!(pos, neg, vol::AbstractVolumePrimitive, transformation, ::Type{Val{:PositiveVolume}}, ::Type{AbstractVolumePrimitive}) = push!(pos, transform(vol, transformation))
decompose!(pos, neg, vol::AbstractVolumePrimitive, transformation, ::Type{Val{:NegativeVolume}}, ::Type{AbstractVolumePrimitive}) = push!(neg, transform(vol, transformation))
decompose!(pos, neg, vol::AbstractVolumePrimitive, transformation, ::Type{Val{:PositiveVolume}}, ::Type{AbstractSurfacePrimitive}) = append!(pos, transform.(get_decomposed_surfaces(vol), transformation))
#decompose!(pos, neg, surf::AbstractSurfacePrimitive, transformation, ::Type{Val{:PositiveVolume}}) = push!(pos, surf)
