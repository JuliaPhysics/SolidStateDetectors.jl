struct Line{T<:AbstractFloat} <: AbstractLinePrimitive{T}
    origin::CartesianPoint{T}
    direction::CartesianVector{T}
end

#Type promotion happens here
function Line(origin::PT, direction::DIR) where {PT<:CartesianPoint, DIR<:CartesianVector}
    eltypes = _csg_get_promoted_eltype.((PT, DIR))
    T = float(promote_type(eltypes...))
    Line{T}(origin, direction)
end

function Line(;
    origin = zero(CartesianPoint{Int}), 
    direction = CartesianVector{Int}(0,0,1)
) 
    Line(origin, direction)
end

function Line{T}(;
    origin = zero(CartesianPoint{Float64}), 
    direction = CartesianVector{Float64}(0,0,1)
) where {T}
    Line{T}(origin, direction)
end

distance(pt::CartesianPoint, l::Line) = norm((pt - l.origin) × l.direction) / norm(l.direction)

function _transform_into_object_coordinate_system(l::Line{T}, p::AbstractPrimitive) where {T}
    origin = _transform_into_object_coordinate_system(l.origin, p) 
    direction = inv(rotation(p)) * l.direction
    Line( origin, CartesianVector(direction) )  
end