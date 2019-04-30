"""
    struct GeometryUnion{T, N, A, B} <: AbstractSet{T, N}

a || b
"""
struct GeometryUnion{T <: Real, N, A, B} <: AbstractSet{T, N}
    # just `Union` seems to conflict with Core.Union... :(
    a::A
    b::B
    translate::Union{CartesianVector{T},Missing}
end

function GeometryUnion{T, N}(a::A, b::B) where {T, N, A, B}
    GeometryUnion{T, N, A, B}(a, b, missing)
end
function GeometryUnion{T, N}(a::A, b::B, translate::CartesianVector{T}) where {T, N, A, B}
    GeometryUnion{T, N, A, B}(a, b, translate)
end

function in(pt::CylindricalPoint{T}, set::GeometryUnion{T, N})::Bool where {T <: Real, N}
    inside::Bool = false
    ismissing(set.translate) ? nothing : pt = CylindricalPoint(CartesianPoint(pt) - set.translate)
    if in(pt, set.a) || in(pt, set.b)
        inside = true
    end
    return inside
end

function in(pt::CartesianPoint{T}, set::GeometryUnion{T, N})::Bool where {T <: Real, N}
    inside::Bool = false
    ismissing(set.translate) ? nothing : pt -= set.translate
    if in(pt, set.a) || in(pt, set.b)
        inside = true
    end
    return inside
end

function Geometry(T::DataType, t::Val{:union}, dict::Dict{Union{Any,String}, Any}, inputunit_dict::Dict{String,Unitful.Units})
    union_of_parts = sum( map(x-> Geometry(T, x, inputunit_dict), dict["parts"]) )
    if haskey(dict,"translate")
        translate = CartesianVector{T}( haskey(dict["translate"],"x") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["x"]) * inputunit_dict["length"] ))) : 0.0,
                                        haskey(dict["translate"],"y") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["y"]) * inputunit_dict["length"] ))) : 0.0,
                                        haskey(dict["translate"],"z") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["z"]) * inputunit_dict["length"] ))) : 0.0)
        return GeometryUnion{T,3}(union_of_parts.a, union_of_parts.b, translate )
    else
        return union_of_parts
    end
end

"""
    struct Intersection{T, N, A, B} <: AbstractSet{T, N}

a && b
"""
struct Intersection{T, N, A, B} <: AbstractSet{T, N}
    a::A
    b::B
    translate::Union{CartesianVector{T},Missing}
end
function Intersection{T, N}(a::A, b::B) where {T, N, A, B}
    Intersection{T, N, A, B}(a, b, missing)
end

function Intersection{T, N}(a::A, b::B, translate::CartesianVector{T}) where {T, N, A, B}
    Intersection{T, N, A, B}(a, b, translate)
end

function in(pt::CartesianPoint{T}, set::Intersection{T, N})::Bool where {T <: Real, N}
    inside::Bool = false
    ismissing(set.translate) ? nothing : pt -= set.translate
    if in(pt, set.a) && in(pt, set.b)
        inside = true
    end
    return inside
end
function in(pt::CylindricalPoint{T}, set::Intersection{T, N})::Bool where {T <: Real, N}
    inside::Bool = false
    ismissing(set.translate) ? nothing : pt = CylindricalPoint(CartesianPoint(pt) - set.translate)
    if in(pt, set.a) && in(pt, set.b)
        inside = true
    end
    return inside
end

function Geometry(T::DataType, t::Val{:intersection}, dict::Dict{Union{Any,String}, Any}, inputunit_dict::Dict{String,Unitful.Units})
    parts = map(x-> Geometry(T, x, inputunit_dict), dict["parts"]) 
    intersection = parts[1]
    for part in parts[2:end]
        intersection = Intersection{T,3}(intersection, part)
    end
    if haskey(dict,"translate")
        translate = CartesianVector{T}( haskey(dict["translate"],"x") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["x"]) * inputunit_dict["length"] ))) : 0.0,
                                        haskey(dict["translate"],"y") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["y"]) * inputunit_dict["length"] ))) : 0.0,
                                        haskey(dict["translate"],"z") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["z"]) * inputunit_dict["length"] ))) : 0.0)
        return Intersection{T,3}(intersection.a, intersection.b, translate )
    else
        return intersection
    end
end


"""
    struct Difference{T, N, A, B} <: AbstractSet{T, N}

a && !b
"""
struct Difference{T, N, A, B} <: AbstractSet{T, N}
    a::A
    b::B
    translate::Union{CartesianVector{T},Missing}
end

function Difference{T, N}(a::A, b::B) where {T, N, A, B}
    Difference{T, N, A, B}(a, b, missing)
end
function Difference{T, N}(a::A, b::B, translate::Union{CartesianVector{T},Missing}) where {T, N, A, B}
    Difference{T, N, A, B}(a, b, translate)
end
function in(pt::CartesianPoint{T}, set::Difference{T, N})::Bool where {T <: Real, N}
    inside::Bool = false
    ismissing(set.translate) ? nothing : pt -= set.translate
    if in(pt, set.a) && !in(pt, set.b)
        inside = true
    end
    return inside
end

function in(pt::CylindricalPoint{T}, set::Difference{T, N})::Bool where {T <: Real, N}
    inside::Bool = false
    ismissing(set.translate) ? nothing : pt = CylindricalPoint(CartesianPoint(pt) - set.translate)
    if in(pt, set.a) && !in(pt, set.b)
        inside = true
    end
    return inside
end

function Geometry(T::DataType, t::Val{:difference}, dict::Dict{Union{Any,String}, Any}, inputunit_dict::Dict{String,Unitful.Units})
    difference = Geometry(T, dict["parts"][1], inputunit_dict) - sum( map(x-> Geometry(T, x, inputunit_dict), dict["parts"][2:end]) )
    if haskey(dict,"translate")
        translate = CartesianVector{T}( haskey(dict["translate"],"x") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["x"]) * inputunit_dict["length"] ))) : 0.0,
                                        haskey(dict["translate"],"y") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["y"]) * inputunit_dict["length"] ))) : 0.0,
                                        haskey(dict["translate"],"z") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["z"]) * inputunit_dict["length"] ))) : 0.0)
        return Difference{T,3}(difference.a, difference.b, translate )
    else
        return difference
    end
end

# function Geometry(T::DataType, t::Val{:difference}, dict::Dict{Union{Any,String}, Any}, inputunit_dict::Dict{String,Unitful.Units})
#     return Geometry(T, dict["parts"][1], inputunit_dict) - sum( map(x-> Geometry(T, x, inputunit_dict), dict["parts"][2:end]) )
# end

function (+)(a::A, b::B)::AbstractSet{T, N} where {T, N, A <: AbstractGeometry{T, N}, B <: AbstractGeometry{T, N}}
    return GeometryUnion{T, N}(a, b)
end

function (-)(a::A, b::B)::AbstractSet{T, N} where {T, N, A <: AbstractGeometry{T, N}, B <: AbstractGeometry{T, N}}
    return Difference{T, N}(a, b)
end

function (&)(a::A, b::B)::AbstractSet{T, N} where {T, N, A <: AbstractGeometry{T, N}, B <: AbstractGeometry{T, N}}
    return Intersection{T, N}(a, b)
end
