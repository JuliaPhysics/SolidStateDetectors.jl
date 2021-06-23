"""
    struct ConstantImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

Returns always a fixed Impurity density.
"""
struct ConstantImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T} 
    ρ::T
end

function get_impurity_density(cdm::ConstantImpurityDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    return cdm.ρ
end


function ImpurityDensity(T::DataType, t::Val{:constant}, dict::Union{Dict{String, Any}, Dict{Any, Any}}, input_units::NamedTuple)
    unit_factor::T = 1
    if haskey(input_units, :length)
        lunit = input_units.length
        unit_factor = inv(ustrip(uconvert( internal_length_unit^3, 1 * lunit^3 )))
    end
    return ConstantImpurityDensity{T}( dict, unit_factor )
end


function ConstantImpurityDensity{T}(dict::Union{Dict{String, Any}, Dict{Any, Any}}, unit_factor::T)::ConstantImpurityDensity{T} where {T <: SSDFloat}
    ρ::T = if haskey(dict, "Impurity_density")   
        # geom_round(unit_factor * T(dict["Impurity_density"]))
        unit_factor * T(dict["Impurity_density"])
    else
        T(0)
    end
    ConstantImpurityDensity{T}( ρ )
end
