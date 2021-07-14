
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
    unit_factor::T = 1
    gradient_unit_factor::T = 1
    if haskey(input_units, :length)
        lunit = input_units.length
        unit_factor = inv(ustrip(uconvert( internal_length_unit^3, 1 * lunit^3 )))
        gradient_unit_factor = inv(ustrip(uconvert( internal_length_unit^4, 1 * lunit^4 )))
    end
    return LinearImpurityDensity{T}( dict, unit_factor, gradient_unit_factor )
end

function LinearImpurityDensity{T}(dict::Union{Dict{String, Any}, Dict{Any, Any}}, unit_factor::T, gradient_unit_factor::T)::LinearImpurityDensity{T} where {T <: SSDFloat}
    offsets, gradients = zeros(T,3), zeros(T,3)
    if prod(map(k -> k in ["x","y","z"], collect(keys(dict)))) @warn "Only x, y and z are supported in the linear Impurity density model.\nChange the Impurity density model in the config file or remove all other entries." end
    # if haskey(dict, "x")     offsets[1] = geom_round(unit_factor * T(dict["x"]["init"]));     gradients[1] = geom_round(gradient_unit_factor * T(dict["x"]["gradient"]))    end
    # if haskey(dict, "y")     offsets[2] = geom_round(unit_factor * T(dict["y"]["init"]));     gradients[2] = geom_round(gradient_unit_factor * T(dict["y"]["gradient"]))    end
    # if haskey(dict, "z")     offsets[3] = geom_round(unit_factor * T(dict["z"]["init"]));     gradients[3] = geom_round(gradient_unit_factor * T(dict["z"]["gradient"]))    end
    if haskey(dict, "x")     offsets[1] = unit_factor * T(dict["x"]["init"]);     gradients[1] = gradient_unit_factor * T(dict["x"]["gradient"])    end
    if haskey(dict, "y")     offsets[2] = unit_factor * T(dict["y"]["init"]);     gradients[2] = gradient_unit_factor * T(dict["y"]["gradient"])    end
    if haskey(dict, "z")     offsets[3] = unit_factor * T(dict["z"]["init"]);     gradients[3] = gradient_unit_factor * T(dict["z"]["gradient"])    end
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
