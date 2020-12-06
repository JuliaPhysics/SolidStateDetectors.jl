function get_decomposed_volumes(vol::AbstractGeometry{T})::Tuple{Vector{<:AbstractGeometry},Vector{<:AbstractGeometry}} where T
    positive_volumes = AbstractGeometry[]
    negative_volumes = AbstractGeometry[]
    translate = CartesianVector{T}(0.0,0.0,0.0)
    decompose_volume(vol, positive_volumes, negative_volumes, translate, 1)
    positive_volumes, negative_volumes
end

function decompose_volume(vol::CSGDifference, pos, neg, translate, flag = 1)::Nothing
    decompose_volume(vol.a, pos, neg, translate, flag)
    decompose_volume(vol.b, pos, neg, translate, flag *= -1)
    nothing
end

function decompose_volume(vol::CSGUnion, pos, neg, translate, flag = 1)::Nothing
    decompose_volume(vol.a ,pos, neg, translate, flag)
    decompose_volume(vol.b, pos, neg, translate, flag)
    nothing
end

function decompose_volume(vol::CSGIntersection, pos, neg, translate, flag = 1)::Nothing
    decompose_volume(vol.a ,pos, neg, translate, flag)
    decompose_volume(vol.b, pos, neg, translate, flag)
    nothing
end

function decompose_volume(vol::TranslatedGeometry, pos, neg, translate, flag = 1)::Nothing
    translate += vol.t
    decompose_volume(vol.p, pos, neg, translate, flag)
    nothing
end

function decompose_volume(vol::AbstractGeometry, pos, neg, translate, flag = 1)::Nothing
    flag == 1 ? push!(pos, vol + translate) : push!(neg, vol + translate)
    nothing
end