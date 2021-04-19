
"""
    struct CylindricalChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}

Simple charge density model which assumes a linear gradient in charge density in each spatial dimension of a cylindrical coordinate system.
`offsets::NTuple{3, T}` are the charge densities at 0 and `gradients::NTuple{3, T}` are the linear slopes in `r` and `z` direction.
"""
struct CylindricalChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}
    offsets::NTuple{3, T}
    gradients::NTuple{3, T}
end

function ChargeDensity(T::DataType, t::Val{:cylindrical}, dict::Union{Dict{String, Any}, Dict{Any, Any}}, input_units::NamedTuple)
    unit_factor::T = 1
    gradient_unit_factor::T = 1
    if haskey(input_units, :length)
        lunit = input_units.length
        unit_factor = inv(ustrip(uconvert( internal_length_unit^3, 1 * lunit^3 )))
        gradient_unit_factor = inv(ustrip(uconvert( internal_length_unit^4, 1 * lunit^4 )))
    end
    return CylindricalChargeDensity{T}( dict, unit_factor, gradient_unit_factor )
end

function CylindricalChargeDensity{T}(dict::Union{Dict{String, Any}, Dict{Any, Any}}, unit_factor::T, gradient_unit_factor::T)::CylindricalChargeDensity{T} where {T <: SSDFloat}
    offsets, gradients = zeros(T,3), zeros(T,3)
    if prod(map(k -> k in ["r","z"], collect(keys(dict)))) @warn "Only r and z are supported in the cylindrical charge density model.\nChange the charge density model in the config file or remove all other entries." end
    if haskey(dict, "r")     offsets[1] = geom_round(unit_factor * T(dict["r"]["init"]));     gradients[1] = geom_round(gradient_unit_factor * T(dict["r"]["gradient"]))    end
    if haskey(dict, "z")     offsets[3] = geom_round(unit_factor * T(dict["z"]["init"]));     gradients[3] = geom_round(gradient_unit_factor * T(dict["z"]["gradient"]))    end
    CylindricalChargeDensity{T}( NTuple{3, T}(offsets), NTuple{3, T}(gradients) )
end

function get_charge_density(lcdm::CylindricalChargeDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    pt::CylindricalPoint{T} = CylindricalPoint(pt)
    ρ::T = 0
    for i in eachindex(lcdm.offsets)
        ρ += (lcdm.offsets[i] + pt[i] * lcdm.gradients[i]) #* T(1e16) # * T(1e10) * T(1e6) -> 1/cm^3 -> 1/m^3
    end
    return ρ
end
