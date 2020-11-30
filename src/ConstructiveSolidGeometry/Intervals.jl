@inline _left_linear_interval(z::T) where {T <: Real} = -z
@inline _left_linear_interval(z::AbstractInterval) = z.left

@inline _right_linear_interval(z::T) where {T <: Real} = z
@inline _right_linear_interval(z::AbstractInterval) = z.right

@inline _width_linear_interval(z::T) where {T <: Real} = 2z
@inline _width_linear_interval(z::AbstractInterval) = width(z)