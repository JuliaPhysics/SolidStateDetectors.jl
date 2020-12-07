function get_decomposed_volumes(vol::AbstractGeometry{T})::Tuple{Vector{<:AbstractGeometry},Vector{<:AbstractGeometry}} where {T}
    positive_volumes = AbstractGeometry[]
    negative_volumes = AbstractGeometry[]
    translate = CartesianVector{T}(0.0,0.0,0.0)
    decompose_volume!(positive_volumes, negative_volumes, vol, translate, Val{:PositiveVolume})
    positive_volumes, negative_volumes
end

@inline invert_volume_type(::Type{Val{:PositiveVolume}}) = Val{:NegativeVolume}
@inline invert_volume_type(::Type{Val{:NegativeVolume}}) = Val{:PositiveVolume}

function decompose_volume!(pos, neg, vol::CSGDifference, translate, ::Type{VT})::Nothing where {VT <:Union{Val{:PositiveVolume}, Val{:NegativeVolume}}}
    decompose_volume!(pos, neg, vol.a, translate, VT)
    decompose_volume!(pos, neg, vol.b, translate, invert_volume_type(VT))
    nothing
end

function decompose_volume!(pos, neg, vol::Union{CSGUnion, CSGIntersection}, translate, ::Type{VT})::Nothing where {VT <:Union{Val{:PositiveVolume}, Val{:NegativeVolume}}}
    decompose_volume!(pos, neg, vol.a, translate, VT)
    decompose_volume!(pos, neg, vol.b, translate, VT)
    nothing
end

function decompose_volume!(pos, neg, vol::TranslatedGeometry, translate, ::Type{VT})::Nothing where {VT <:Union{Val{:PositiveVolume}, Val{:NegativeVolume}}}
    translate += vol.t
    decompose_volume!(pos, neg, vol.p, translate, VT)
    nothing
end

decompose_volume!(pos, neg, vol::Union{AbstractVolumePrimitive, RotatedGeometry, ScaledGeometry}, translate, ::Type{Val{:PositiveVolume}}) = push!(pos, vol + translate)
decompose_volume!(pos, neg, vol::Union{AbstractVolumePrimitive, RotatedGeometry, ScaledGeometry}, translate, ::Type{Val{:NegativeVolume}}) = push!(neg, vol + translate)