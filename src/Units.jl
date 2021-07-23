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

_get_TDU(x::Quantity{T,D,U}) where {T,D,U} = T, D, U
to_internal_units(x::Quantity) = ustrip(uconvert(u_internal, x))
to_internal_units(x::Quantity{<:Real, Unitful.𝐓}) = ustrip(uconvert(internal_time_unit, x))
to_internal_units(x::Quantity{<:Real, _get_TDU(1u"V")[2]}) = ustrip(uconvert(internal_voltage_unit, x))
to_internal_units(x::Quantity{<:Real, _get_TDU(1u"V / m")[2]}) = ustrip(uconvert(internal_efield_unit, x))
to_internal_units(x::Quantity{<:Real, _get_TDU(1u"eV")[2]}) = ustrip(uconvert(internal_energy_unit, x))

from_internal_units(u_external::typeof(Unitful.NoUnits), u_internal::Unitful.Units, x::Real) where {T<:Real} = x
from_internal_units(u_external::Unitful.Units, u_internal::Unitful.Units, x::Real) where {T<:Quantity} = uconvert(u_external, x * u_internal)

from_internal_units(u_external::typeof(Unitful.NoUnits), u_internal::Unitful.Units, x::AbstractArray{<:Real}) where {T<:Real} = x
from_internal_units(u_external::Unitful.Units, u_internal::Unitful.Units, x::AbstractArray{<:Real}) where {T<:Quantity} = uconvert.(u_external, x * u_internal)

unit_conversion = Dict{String, Unitful.Units}(
    "nm" => u"nm", "um" => u"μm", "mm" => u"mm", "cm" => u"cm", "m" => u"m", #length
    "deg" => u"°","rad" => u"rad", #angle
    "V" => u"V", "kV" => u"kV", #potential
    "K" => u"K", "Kelvin" => u"K", "C" => u"°C", "Celsius" => u"°C", #temperature
)

const UnitTuple = NamedTuple{(:length, :angle, :potential, :temperature), Tuple{
                    Unitful.FreeUnits{<:Any, Unitful.𝐋}, 
                    Unitful.FreeUnits{<:Any, NoDims}, 
                    Unitful.FreeUnits{<:Any, Unitful.𝐋^2 * Unitful.𝐌 * Unitful.𝐈^-1 * Unitful.𝐓^-3}, 
                    Unitful.FreeUnits{<:Any, Unitful.𝚯}}}

function default_unit_tuple()::UnitTuple
    return (
        length = u"m",
        angle = u"°",
        potential = u"V",
        temperature = u"K"
    )
end

function construct_unit(ustring::String)::Unitful.Units
    try
        parsed_unit = uparse(ustring)
        @assert parsed_unit isa Unitful.Units
        return parsed_unit
    catch 
        @assert ustring in keys(unit_conversion) "Unit string $(ustring) cannot be interpreted..."
        return unit_conversion[ustring]
    end
end

function construct_units(config_file_dict::AbstractDict)::UnitTuple
    dunits::NamedTuple = default_unit_tuple()
    if haskey(config_file_dict, "units")
        d = config_file_dict["units"]
        dunits = (
            length = haskey(d, "length") ? construct_unit(d["length"]) : dunits.length, 
            angle  = haskey(d, "angle") ? construct_unit(d["angle"]) : dunits.angle,
            potential = haskey(d, "potential") ? construct_unit(d["potential"]) : dunits.potential,
            temperature = haskey(d, "temperature") ? construct_unit(d["temperature"]) : dunits.temperature
        )
    end
    matching_units::Bool = dimension.(values(dunits)) == (Unitful.𝐋, Unitful.NoDims, Unitful.𝐋^2 * Unitful.𝐌 * Unitful.𝐈^-1 * Unitful.𝐓^-3, Unitful.𝚯) 
    @assert matching_units "Some of the units are not given in the right format:"*
        (dimension(dunits.length)      == Unitful.𝐋       ? "" : "\n  - units.length ($(dunits.length)) does not describe a length.")*
        (dimension(dunits.angle)       == Unitful.NoDims  ? "" : "\n  - units.angle ($(dunits.angle)) does not describe an angle.")*
        (dimension(dunits.potential)   == Unitful.𝐋^2 * Unitful.𝐌 * Unitful.𝐈^-1 * Unitful.𝐓^-3 ? "" : "\n  - units.potential ($(dunits.potential)) does not describe a potential.")*
        (dimension(dunits.temperature) == Unitful.𝚯       ? "" : "\n  - units.temperature ($(dunits.temperature)) does not describe a temperature.")
        
    return UnitTuple(dunits)
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