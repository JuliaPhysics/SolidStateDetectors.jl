function get_decomposed_volumes(vol::AbstractGeometry{T})::Tuple{Vector{<:AbstractGeometry},Vector{<:AbstractGeometry}} where {T}
    positive_volumes = AbstractGeometry[]
    negative_volumes = AbstractGeometry[]
    transformation = CSGTransformation[]
    decompose_volume!(positive_volumes, negative_volumes, vol, transformation, Val{:PositiveVolume})
    positive_volumes, negative_volumes
end

@inline invert_volume_type(::Type{Val{:PositiveVolume}}) = Val{:NegativeVolume}
@inline invert_volume_type(::Type{Val{:NegativeVolume}}) = Val{:PositiveVolume}

function decompose_volume!(pos, neg, vol::CSGDifference, transformation, ::Type{VT})::Nothing where {VT <:Union{Val{:PositiveVolume}, Val{:NegativeVolume}}}
    decompose_volume!(pos, neg, vol.a, transformation, VT)
    decompose_volume!(pos, neg, vol.b, transformation, invert_volume_type(VT))
    nothing
end

function decompose_volume!(pos, neg, vol::Union{CSGUnion, CSGIntersection}, transformation, ::Type{VT})::Nothing where {VT <:Union{Val{:PositiveVolume}, Val{:NegativeVolume}}}
    decompose_volume!(pos, neg, vol.a, transformation, VT)
    decompose_volume!(pos, neg, vol.b, transformation, VT)
    nothing
end

function decompose_volume!(pos, neg, vol::TranslatedGeometry, transformation, ::Type{VT})::Nothing where {VT <:Union{Val{:PositiveVolume}, Val{:NegativeVolume}}}
    push!(transformation, vol.t)
    decompose_volume!(pos, neg, vol.p, transformation, VT)
    pop!(transformation)
    nothing
end

function decompose_volume!(pos, neg, vol::RotatedGeometry, transformation, ::Type{VT})::Nothing where {VT <:Union{Val{:PositiveVolume}, Val{:NegativeVolume}}}
    push!(transformation, inv(vol.inv_r))
    decompose_volume!(pos, neg, vol.p, transformation, VT)
    pop!(transformation)
    nothing
end

function decompose_volume!(pos, neg, vol::ScaledGeometry, transformation, ::Type{VT})::Nothing where {VT <:Union{Val{:PositiveVolume}, Val{:NegativeVolume}}}
    push!(transformation, inv.(vol.inv_s))
    decompose_volume!(pos, neg, vol.p, transformation, VT)
    pop!(transformation)
    nothing
end

decompose_volume!(pos, neg, vol::AbstractVolumePrimitive, transformation, ::Type{Val{:PositiveVolume}}) = push!(pos, transform(vol, transformation))
decompose_volume!(pos, neg, vol::AbstractVolumePrimitive, transformation, ::Type{Val{:NegativeVolume}}) = push!(neg, transform(vol, transformation))