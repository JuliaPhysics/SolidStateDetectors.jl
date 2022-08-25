struct Line{T} <: AbstractLinePrimitive{T}
    origin::CartesianPoint{T}
    direction::CartesianVector{T}
    
    Line{T}(origin, direction) where T = new{T}(origin,direction)
    function Line(origin::PT, direction::DIR) where {PT<:CartesianPoint, DIR<:CartesianVector}
        #Type promotion happens here
        eltypes = _csg_get_promoted_eltype.((PT, DIR))
        T = float(promote_type(eltypes...))
        Line{T}(origin, direction)
    end
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
    direction = _transform_into_object_coordinate_system(l.direction, p)
    Line( origin, direction )  
end