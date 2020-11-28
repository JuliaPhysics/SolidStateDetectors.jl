@inline _left(r::T) where {T <: Real} = -r
@inline _left(r::AbstractInterval) = r.left

@inline _right(r::T) where {T <: Real} = r
@inline _right(r::AbstractInterval) = r.right

@inline _width(z::T) where {T <: Real} = T(2) * z
@inline _width(z::AbstractInterval) = width(z)