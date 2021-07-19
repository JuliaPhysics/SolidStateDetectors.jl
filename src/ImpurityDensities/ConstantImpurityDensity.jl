"""
    struct ConstantImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

Returns always a fixed impurity density.
"""
struct ConstantImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T} 
    ρ::T
end

function get_impurity_density(cdm::ConstantImpurityDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    return cdm.ρ
end

function ImpurityDensity(T::DataType, t::Val{:constant}, dict::Union{Dict{String, Any}, Dict{Any, Any}}, input_units::NamedTuple)
    ρ::T = haskey(dict, "value") ? _parse_value(T, dict["value"], input_units.length^(-3)) : T(0)
    ConstantImpurityDensity{T}( ρ )
end