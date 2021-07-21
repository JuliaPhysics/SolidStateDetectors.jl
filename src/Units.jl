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

_get_TDU(x::Quantity{T,D,U}) where {T,D,U} = T, D, U
const length_unit = Quantity{<:Real, Unitful.𝐋}
const angle_unit = Quantity{<:Real, NoDims, <:Union{_get_TDU(1u"rad")[3], _get_TDU(1u"°")[3]}}

unit_conversion = Dict{String, Unitful.Units}(
    "nm" => u"nm", "um" => u"μm", "mm" => u"mm", "cm" => u"cm", "m" => u"m", #length
    "deg" => u"°","rad" => u"rad", #angle
    "V" => u"V", "kV" => u"kV", #potential
    "K" => u"K", "Kelvin" => u"K", "C" => u"°C", "Celsius" => u"°C", #temperature
)

function default_unit_tuple()::NamedTuple{<:Any, <:NTuple{4, Unitful.Units}}
    return (
        length = u"m", # change this to u"m" ? SI Units
        potential = u"V",
        angle = u"°",
        temperature = u"K"
    )
end

function construct_units(config_file_dict::AbstractDict)
    dunits::NamedTuple = default_unit_tuple()
    if haskey(config_file_dict, "units")
        d = config_file_dict["units"]
        dunits = (
            length = haskey(d, "length") ? unit_conversion[d["length"]] : dunits.length, 
            angle  = haskey(d, "angle") ? unit_conversion[d["angle"]] : dunits.angle,
            potential = haskey(d, "potential") ? unit_conversion[d["potential"]] : dunits.potential,
            temperature = haskey(d, "temperature") ? unit_conversion[d["temperature"]] : dunits.temperature
        )
    end
    dunits
end



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