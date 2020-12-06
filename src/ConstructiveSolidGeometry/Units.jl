const internal_length_unit      = u"m"
const internal_angle_unit       = u"rad"
const internal_time_unit        = u"s"
const internal_voltage_unit     = u"V"
const internal_efield_unit      = u"V / m"
const internal_energy_unit      = u"eV"
const internal_temperature_unit = u"K"

const MaybeWithUnits{T} = Union{T, Quantity{<:T}}
const RealQuantity = MaybeWithUnits{<:Real}

to_internal_units(u_internal::Unitful.Units, x::Real) = x
to_internal_units(u_internal::Unitful.Units, x::Quantity) = ustrip(uconvert(u_internal, x))

to_internal_units(u_internal::Unitful.Units, x::AbstractArray{<:Real}) = x
to_internal_units(u_internal::Unitful.Units, x::AbstractArray{<:Quantity}) = ustrip.(uconvert.(u_internal, x))

from_internal_units(u_external::typeof(Unitful.NoUnits), u_internal::Unitful.Units, x::Real) where {T<:Real} = x
from_internal_units(u_external::Unitful.Units, u_internal::Unitful.Units, x::Real) where {T<:Quantity} = uconvert(u_external, x * u_internal)

from_internal_units(u_external::typeof(Unitful.NoUnits), u_internal::Unitful.Units, x::AbstractArray{<:Real}) where {T<:Real} = x
from_internal_units(u_external::Unitful.Units, u_internal::Unitful.Units, x::AbstractArray{<:Real}) where {T<:Quantity} = uconvert.(u_external, x * u_internal)