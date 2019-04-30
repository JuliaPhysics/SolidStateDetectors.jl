# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

# Note: Most (if not all) of these types are temporary, to be replaced later
# on by generic types from RadiationDetectorSignals.jl (when the latter is
# mature enough).

# Internal units should be SI units
const internal_length_unit  = u"m"
const internal_angle_unit   = u"rad"
const internal_time_unit    = u"s"
const internal_voltage_unit = u"V"
const internal_efield_unit  = u"V / m"
const internal_energy_unit  = u"eV"

to_internal_units(u_internal::Unitful.Units, x::Real) = x
to_internal_units(u_internal::Unitful.Units, x::Quantity) = ustrip(uconvert(u_internal, x))

to_internal_units(u_internal::Unitful.Units, x::AbstractArray{<:Real}) = x
to_internal_units(u_internal::Unitful.Units, x::AbstractArray{<:Quantity}) = ustrip.(uconvert.(u_internal, x))

from_internal_units(u_external::typeof(Unitful.NoUnits), u_internal::Unitful.Units, x::Real) where {T<:Real} = x
from_internal_units(u_external::Unitful.Units, u_internal::Unitful.Units, x::Real) where {T<:Quantity} = uconvert(u_external, x * u_internal)

from_internal_units(u_external::typeof(Unitful.NoUnits), u_internal::Unitful.Units, x::AbstractArray{<:Real}) where {T<:Real} = x
from_internal_units(u_external::Unitful.Units, u_internal::Unitful.Units, x::AbstractArray{<:Real}) where {T<:Quantity} = uconvert.(u_external, x * u_internal)


const MaybeWithUnits{T} = Union{T, Quantity{<:T}}
const RealQuantity = MaybeWithUnits{<:Real}


DetectorHitEvents = TypedTables.Table{
    <:NamedTuple{
        (:evtno, :detno, :thit, :edep, :pos),
        <:Tuple{
            Union{Integer, AbstractVector{<:Integer}},
            Union{Integer, AbstractVector{<:Integer}},
            AbstractVector{<:RealQuantity},
            AbstractVector{<:RealQuantity},
            AbstractVector{<:AbstractVector{<:RealQuantity}}
        }
    }
}