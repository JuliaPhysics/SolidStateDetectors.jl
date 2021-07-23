const internal_unit_length     = u"m"
const internal_unit_angle      = u"rad"

_get_TDU(x::Quantity{T,D,U}) where {T,D,U} = T, D, U
const length_unit = Quantity{<:Real, Unitful.𝐋}
const angle_unit = Quantity{<:Real, NoDims, <:Union{_get_TDU(1u"rad")[3], _get_TDU(1u"°")[3]}}

to_internal_units(x::Real) = x
to_internal_units(x::length_unit) = ustrip(uconvert(internal_unit_length, x))
to_internal_units(x::angle_unit) = ustrip(uconvert(internal_unit_angle, x))
to_internal_units(x::Quantity{<:Real, Unitful.𝐋^(-3)}) = ustrip(uconvert(internal_unit_length^(-3), x)) # densities
to_internal_units(x::Quantity{<:Real, Unitful.𝐋^(-4)}) = ustrip(uconvert(internal_unit_length^(-4), x)) # density gradients
to_internal_units(x::AbstractArray) = to_internal_units.(x)