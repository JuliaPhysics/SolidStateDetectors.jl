"""
    struct ConstantChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}

Returns always a fixed charge density.
"""
struct ConstantChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T} 
    ρ::T
end

function get_charge_density(cdm::ConstantChargeDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    return cdm.ρ
end


function ChargeDensity(T::DataType, t::Val{:constant}, dict::Union{Dict{String, Any}, Dict{Any, Any}}, input_units::NamedTuple)
    ρ::T = haskey(dict, "value") ? _parse_value(T, dict["value"], input_units.length^(-3)) : T(0)
    ConstantChargeDensity{T}( ρ )
end
