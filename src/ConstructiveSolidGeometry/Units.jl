const internal_length_unit    = u"m"
const internal_angle_unit     = u"rad"

_get_TDU(x::Quantity{T,D,U}) where {T,D,U} = T, D, U
const LengthQuantity = Quantity{<:Real, Unitful.𝐋}
const AngleQuantity = Quantity{<:Real, NoDims, <:Union{_get_TDU(1u"rad")[3], _get_TDU(1u"°")[3]}}

to_internal_units(x::Real) = x
to_internal_units(x::LengthQuantity) = ustrip(internal_length_unit, x)
to_internal_units(x::AngleQuantity)  = ustrip(internal_angle_unit, x)
to_internal_units(x::Quantity{<:Real, Unitful.𝐋^(-1)}) = ustrip(internal_length_unit^(-1), x) # charge trapping
to_internal_units(x::Quantity{<:Real, Unitful.𝐋^(-3)}) = ustrip(internal_length_unit^(-3), x) # densities
to_internal_units(x::Quantity{<:Real, Unitful.𝐋^(-4)}) = ustrip(internal_length_unit^(-4), x) # density gradients
to_internal_units(x::Quantity{<:Real, Unitful.𝐋^(-5)}) = ustrip(internal_length_unit^(-5), x) # density quadratic components
to_internal_units(x::AbstractArray) = to_internal_units.(x)

from_internal_units(x::Real, unit::Unitful.Units{<:Any, Unitful.𝐋}) = uconvert(unit, x * internal_length_unit)
from_internal_units(x::Real, unit::Unitful.Units{<:Any, Unitful.NoDims}) = uconvert(unit, x * internal_angle_unit)
from_internal_units(x::AbstractArray, unit::Unitful.Units) = broadcast(x -> from_internal_units(x, unit), x)