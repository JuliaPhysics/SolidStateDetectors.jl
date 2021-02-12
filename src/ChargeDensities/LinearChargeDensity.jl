
"""
    struct LinearChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}

Simple charge density model which assumes a linear gradient in charge density in each dimension of a Cartesian coordinate system.
`offsets::NTuple{3, T}` are the charge densities at 0 and `gradients::NTuple{3, T}` are the linear slopes in each dimension, `x`, `y` and `z`.
"""
struct LinearChargeDensity{T <: SSDFloat} <: AbstractChargeDensity{T}
    offsets::NTuple{3, T}
    gradients::NTuple{3, T}
end

function ChargeDensity(T::DataType, t::Val{:linear}, dict::Union{Dict{String, Any}, Dict{Any, Any}}, inputunit_dict::Dict)
    unit_factor::T = 1
    gradient_unit_factor::T = 1
    if haskey(inputunit_dict, "length")
        lunit = inputunit_dict["length"]
        unit_factor = inv(ustrip(uconvert( internal_length_unit^3, 1 * lunit^3 )))
        gradient_unit_factor = inv(ustrip(uconvert( internal_length_unit^4, 1 * lunit^4 )))
    end
    return LinearChargeDensity{T}( dict, unit_factor, gradient_unit_factor )
end

function LinearChargeDensity{T}(dict::Union{Dict{String, Any}, Dict{Any, Any}}, unit_factor::T, gradient_unit_factor::T)::LinearChargeDensity{T} where {T <: SSDFloat}
    offsets, gradients = zeros(T,3), zeros(T,3)
    if prod(map(k -> k in ["x","y","z"], collect(keys(dict)))) @warn "Only x, y and z are supported in the linear charge density model.\nChange the charge density model in the config file or remove all other entries." end
    if haskey(dict, "x")     offsets[1] = geom_round(unit_factor * T(dict["x"]["init"]));     gradients[1] = geom_round(gradient_unit_factor * T(dict["x"]["gradient"]))    end
    if haskey(dict, "y")     offsets[2] = geom_round(unit_factor * T(dict["y"]["init"]));     gradients[2] = geom_round(gradient_unit_factor * T(dict["y"]["gradient"]))    end
    if haskey(dict, "z")     offsets[3] = geom_round(unit_factor * T(dict["z"]["init"]));     gradients[3] = geom_round(gradient_unit_factor * T(dict["z"]["gradient"]))    end
    LinearChargeDensity{T}( NTuple{3, T}(offsets), NTuple{3, T}(gradients) )
end

function get_charge_density(lcdm::LinearChargeDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    pt::CartesianPoint{T} = CartesianPoint(pt)
    ρ::T = 0
    for i in eachindex(lcdm.offsets)
        ρ += (lcdm.offsets[i] + pt[i] * lcdm.gradients[i]) #* T(1e16) # * T(1e10) * T(1e6) -> 1/cm^3 -> 1/m^3
    end
    return ρ
end
