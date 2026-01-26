# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

# Note: Most (if not all) of these types are temporary, to be replaced later
# on by generic types from RadiationDetectorSignals.jl (when the latter is
# mature enough).

# Internal units should be SI units
const internal_length_unit  = ConstructiveSolidGeometry.internal_length_unit
const internal_angle_unit   = ConstructiveSolidGeometry.internal_angle_unit
const internal_time_unit    = u"s"
const internal_voltage_unit = u"V"
const internal_energy_unit  = u"eV"
const internal_efield_unit  = internal_voltage_unit / internal_length_unit
const internal_charge_unit  = u"C"
const internal_mobility_unit = u"m^2/(V*s)"
const internal_diffusion_unit = internal_length_unit ^ 2 / internal_time_unit
const internal_charge_density_unit = internal_charge_unit / internal_length_unit ^ 3
const internal_charge_density_gradient_unit = internal_charge_unit / internal_length_unit ^ 4
const internal_temperature_unit = u"K"
const internal_activation_energy_unit = u"cal/mol"

const external_charge_unit  = u"e_au" # elementary charge - from UnitfulAtomic.jl

to_internal_units(x::Quantity{<:Real}) = throw(ArgumentError("Unit $(unit(x)) unknown to SolidStateDetectors.jl"))
to_internal_units(x::Quantity{<:Real, dimension(internal_time_unit)})    = ustrip(internal_time_unit,    x)
to_internal_units(x::Quantity{<:Real, dimension(internal_voltage_unit)}) = ustrip(internal_voltage_unit, x)
to_internal_units(x::Quantity{<:Real, dimension(internal_efield_unit)})  = ustrip(internal_efield_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_energy_unit)})  = ustrip(internal_energy_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_charge_unit)})  = ustrip(internal_charge_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_mobility_unit)})  = ustrip(internal_mobility_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_diffusion_unit)})  = ustrip(internal_diffusion_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_charge_density_unit)})  = ustrip(internal_charge_density_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_charge_density_gradient_unit)}) = ustrip(internal_charge_density_gradient_unit, x)
to_internal_units(x::Quantity{<:Real, dimension(internal_temperature_unit)})     = ustrip(internal_temperature_unit,  x)
to_internal_units(x::Quantity{<:Real, dimension(internal_activation_energy_unit)}) = ustrip(internal_activation_energy_unit,  x)

from_internal_units(x::Real, unit::Unitful.Units{<:Any, dimension(internal_time_unit)})    = uconvert(unit, x * internal_time_unit)
from_internal_units(x::Real, unit::Unitful.Units{<:Any, dimension(internal_voltage_unit)}) = uconvert(unit, x * internal_voltage_unit)
from_internal_units(x::Real, unit::Unitful.Units{<:Any, dimension(internal_efield_unit)})  = uconvert(unit, x * internal_efield_unit)
from_internal_units(x::Real, unit::Unitful.Units{<:Any, dimension(internal_energy_unit)})  = uconvert(unit, x * internal_energy_unit)
from_internal_units(x::Real, unit::Unitful.Units{<:Any, dimension(internal_charge_unit)})  = uconvert(unit, x * internal_charge_unit)
from_internal_units(x::Real, unit::Unitful.Units{<:Any, dimension(internal_mobility_unit)})  = uconvert(unit, x * internal_mobility_unit)
from_internal_units(x::Real, unit::Unitful.Units{<:Any, dimension(internal_diffusion_unit)})  = uconvert(unit, x * internal_diffusion_unit)
from_internal_units(x::Real, unit::Unitful.Units{<:Any, dimension(internal_temperature_unit)}) = uconvert(unit, x * internal_temperature_unit)

# Internal function for now: We should also fano noise here (optionally)
_convert_internal_energy_to_external_charge(material) = inv(to_internal_units(material.E_ionisation)) * external_charge_unit
_convert_internal_energy_to_external_unit(unit::Unitful.Units{<:Any, dimension(internal_charge_unit)}, material) = uconvert(unit, _convert_internal_energy_to_external_charge(material))
_convert_internal_energy_to_external_unit(unit::Unitful.Units{<:Any, dimension(internal_energy_unit)}, material) = uconvert(unit, 1*internal_energy_unit)
_convert_internal_energy_to_external_unit(unit::Unitful.Units, material) = throw(ArgumentError("The unit $(unit) is neither charge nor energy."))

unit_conversion = Dict{String, Unitful.Units}(
    "um" => u"μm", # length
    "deg" => u"°", # angle
    "Kelvin" => u"K", "Celsius" => u"°C", #temperature
    "e" => u"e_au" # elementary charge
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
        parsed_unit = uparse(ustring, unit_context=[Unitful, UnitfulAtomic])
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

# If this is modified --> update the docstring of `simulate_waveforms`
const DHE_column_names = (:evtno, :detno, :thit, :edep, :pos)
const DHE_column_types = (
    Union{<:Integer, <:AbstractVector{<:Integer}},
    Union{<:Integer, <:AbstractVector{<:Integer}},
    AbstractVector{<:Union{<:Real, Unitful.Time}},
    AbstractVector{<:Union{<:Real, Unitful.Energy}},
    AbstractVector{<:Union{<:AbstractVector{<:RealQuantity}, <:AbstractCoordinatePoint}}
)
DetectorHitEvents = TypedTables.Table{<:NamedTuple{DHE_column_names, <:Tuple{DHE_column_types...}}}

@inline get_detector_hits_table(t::DetectorHitEvents) = t
@inline get_detector_hits_table(t::TypedTables.Table) = TypedTables.getproperties(t, DHE_column_names)

@inline is_detector_hits_table(t::TypedTables.Table) = begin
    for (name, type) in zip(DHE_column_names, DHE_column_types)
        hasproperty(t, name) || throw(ArgumentError("Expected detector hit events table to have column named `$(name)`"))
        eltype(getproperty(t, name)) <: type || throw(ArgumentError("Expected detector hit events table column `$(name)` entries to be of type `$(type)`"))
    end
    return true
end