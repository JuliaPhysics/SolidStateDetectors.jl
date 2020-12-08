const internal_unit_length     = u"m"
const internal_unit_angle      = u"rad"

to_internal_units(x::Real) = x
to_internal_units(x::Quantity{<:Real, Unitful.𝐋}) = ustrip(uconvert(internal_unit_length, x))
to_internal_units(x::Quantity{<:Real, Unitful.NoDims}) = ustrip(uconvert(internal_unit_angle, x))
to_internal_units(x::AbstractArray{<:Quantity}) = to_internal_units.(x)