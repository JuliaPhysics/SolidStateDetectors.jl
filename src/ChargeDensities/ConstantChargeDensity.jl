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
    unit_factor::T = 1
    if haskey(input_units, :length) 
        lunit = input_units.length
        unit_factor = inv(ustrip(uconvert( internal_length_unit^3, 1 * lunit^3 )))
    end
    return ConstantChargeDensity{T}( dict, unit_factor )
end


function ConstantChargeDensity{T}(dict::Union{Dict{String, Any}, Dict{Any, Any}}, unit_factor::T)::ConstantChargeDensity{T} where {T <: SSDFloat}
    ρ::T = if haskey(dict, "charge_density")   
        # geom_round(unit_factor * T(dict["charge_density"]))
        unit_factor * T(dict["charge_density"])
    else
        T(0)
    end
    ConstantChargeDensity{T}( ρ )
end
