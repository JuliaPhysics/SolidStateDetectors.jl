
"""
    struct LinearImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

Simple Impurity density model which assumes a linear gradient in Impurity density in each dimension of a Cartesian coordinate system.
`offsets::NTuple{3, T}` are the Impurity densities at 0 and `gradients::NTuple{3, T}` are the linear slopes in each dimension, `x`, `y` and `z`.
"""
struct LinearImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}
    offsets::NTuple{3, T}
    gradients::NTuple{3, T}
end

function ImpurityDensity(T::DataType, t::Val{:linear}, dict::Union{Dict{String, Any}, Dict{Any, Any}}, input_units::NamedTuple)
    offsets, gradients = zeros(T,3), zeros(T,3)
    density_unit = input_units.length^(-3)
    density_gradient_unit = input_units.length^(-4)
    if prod(map(k -> k in ["x","y","z"], collect(keys(dict)))) @warn "Only x, y and z are supported in the linear impurity density model.\nChange the impurity density model in the config file or remove all other entries." end
    if haskey(dict, "x")     
        if haskey(dict["x"], "init")     offsets[1]   = _parse_value(T, dict["x"]["init"], density_unit) end
        if haskey(dict["x"], "gradient") gradients[1] = _parse_value(T, dict["x"]["gradient"], density_gradient_unit) end
    end
    if haskey(dict, "y")     
        if haskey(dict["y"], "init")     offsets[2]   = _parse_value(T, dict["y"]["init"], density_unit) end
        if haskey(dict["y"], "gradient") gradients[2] = _parse_value(T, dict["y"]["gradient"], density_gradient_unit) end
    end
    if haskey(dict, "z")     
        if haskey(dict["z"], "init")     offsets[3]   = _parse_value(T, dict["z"]["init"], density_unit) end
        if haskey(dict["z"], "gradient") gradients[3] = _parse_value(T, dict["z"]["gradient"], density_gradient_unit) end
    end
    LinearImpurityDensity{T}( NTuple{3, T}(offsets), NTuple{3, T}(gradients) )
end

function get_impurity_density(lcdm::LinearImpurityDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    pt::CartesianPoint{T} = CartesianPoint(pt)
    ρ::T = 0
    for i in eachindex(lcdm.offsets)
        ρ += (lcdm.offsets[i] + pt[i] * lcdm.gradients[i]) #* T(1e16) # * T(1e10) * T(1e6) -> 1/cm^3 -> 1/m^3
    end
    return ρ
end
