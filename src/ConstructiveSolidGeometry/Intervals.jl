@inline _left_linear_interval(z::Real) = -z
@inline _left_linear_interval(z::AbstractInterval) = z.left

@inline _right_linear_interval(z::Real) = z
@inline _right_linear_interval(z::AbstractInterval) = z.right

@inline _width_linear_interval(z::Real) = 2z
@inline _width_linear_interval(z::AbstractInterval) = width(z)


@inline _left_radial_interval(r::T) where {T <: Real} = T(0)
@inline _left_radial_interval(r::AbstractInterval) = r.left

@inline _right_radial_interval(r::Real) = r
@inline _right_radial_interval(r::AbstractInterval) = r.right