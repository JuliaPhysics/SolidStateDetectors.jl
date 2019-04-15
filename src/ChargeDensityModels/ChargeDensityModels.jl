abstract type AbstractChargeDensityModel{T <: SSDFloat} end

@inline function ChargeDensityModel(T::DataType, dict::Dict{String, Any})
    return ChargeDensityModel(T, Val{Symbol(dict["name"])}(), dict)
end


"""
    struct LinearChargeDensityModel{T <: SSDFloat} <: AbstractChargeDensityModel{T}

Simple charge density model which assumes a linear gradient in charge density in each dimension.
`offsets::NTuple{3, T}` are the charge densities at 0 and `gradients::NTuple{3, T}` are the linear
slopes in each dimension. 
"""
struct LinearChargeDensityModel{T <: SSDFloat} <: AbstractChargeDensityModel{T}
    offsets::NTuple{3, T}
    gradients::NTuple{3, T}
end

function get_charge_density(lcdm::LinearChargeDensityModel{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    ρ::T = 0
    for i in eachindex(lcdm.offsets)
        ρ += (lcdm.offsets[i] + pt[i] * lcdm.gradients[i]) * T(1e16) # * T(1e10) * T(1e6) -> 1/cm^3 -> 1/m^3
    end
    return ρ
end


function ChargeDensityModel(T::DataType, t::Val{:Linear}, dict::Dict{String, Any})
    return LinearChargeDensityModel{T}( dict )
end

function LinearChargeDensityModel{T}(dict::Dict{String, Any})::LinearChargeDensityModel{T} where {T <: SSDFloat}
    offsets, gradients = if "r_0" in keys(dict) # :Cylindrical
        (dict["r_0"], dict["phi_0"], dict["z_0"]), (dict["r_gradient"], dict["phi_gradient"], dict["z_gradient"])
    else # :Cartesian
        (dict["x_0"], dict["y_0"], dict["z_0"]), (dict["x_gradient"], dict["y_gradient"], dict["z_gradient"])
    end
    LinearChargeDensityModel{T}( offsets, gradients )
end

