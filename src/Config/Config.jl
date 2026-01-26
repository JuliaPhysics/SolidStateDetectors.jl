# """
#     abstract type AbstractConfig{T <: SSDFloat} end
# 
# Supertype of all detector/world/object configs.
# 
# User defined geometries must be subtype of `AbstractConfig{T}`.
# 
# There are a few functions which must be defined for a user config, e.g. `struct UserConfig{T} <: AbstractConfig{T}`:
# 
# For cylindrical grids:
# 
# * in(pt::CylindricalPoint{T}, config::UserConfig{T})::Bool where {T <: SSDFloat}
# * Grid(config::UserConfig{T})::CylindricalGrid{T} where {T <: SSDFloat}
# * get\\_ρ\\_and\\_ϵ(pt::CylindricalPoint{T}, config::UserConfig{T})::Tuple{T, T} where {T <: SSDFloat} 
# * set\\_point\\_types\\_and\\_fixed\\_potentials!(point_types::Array{PointType, 3}, potential::Array{T, 3}, 
#         grid::CylindricalGrid{T}, config::UserConfig{T}; weighting\\_potential\\_contact\\_id::Union{Missing, Int} = missing)::Nothing where {T <: SSDFloat}
# 
# For cartesian grids:
# 
# * in(pt::CartesianPoint{3, T}, config::UserConfig{T})::Bool 
# * Grid(config::UserConfig{T})::CartesianGrid3D{T} where {T <: SSDFloat}
# * get\\_ρ\\_and\\_ϵ(pt::CartesianPoint{3, T}, config::UserConfig{T})::Tuple{T, T} where {T <: SSDFloat} 
# * set\\_point\\_types\\_and\\_fixed\\_potentials!(point_types::Array{PointType, 3}, potential::Array{T, 3}, 
#         grid::CartesianGrid3D{T}, config::UserConfig{T}; weighting\\_potential\\_contact\\_id::Union{Missing, Int} = missing)::Nothing where {T <: SSDFloat}
# 
# """
abstract type AbstractConfig{T <: SSDFloat} end

function _update_density_config!(cf::AbstractDict, density::AbstractString, units::UnitTuple = default_unit_tuple())
    
    OLD_KEYS = Dict("linear" => ("x", "y", "z"), "cylindrical" => ("r", "z"))[cf["name"]]
    NEW_KEYS = ("offset", "gradient")
    UNIT = Dict("impurity" => units.length^(-3), "charge" => internal_charge_unit * units.length^(-3))[density]
    UNIT_STRING = Dict("impurity" => string(units.length^(-3)), "charge" => string(internal_charge_unit) * "*" * string(units.length^(-3)))[density]

    if any(haskey.(Ref(cf), OLD_KEYS)) # still in the old format
        msg = """Deprecation warning: Since v0.11.0, the configuration syntax for the $(cf["name"]) $(density) density has been updated.
        The following syntax is no longer valid (encountered unexpected fields $(filter(k -> haskey(cf, k), OLD_KEYS))):\n
        $(YAML.write(Dict("$(density)_density" => cf)))\n"""

        # convert config dictionary into the correct format
        for k in OLD_KEYS
            if haskey(cf, k)
                v = pop!(cf, k)
                # combine all inits together
                if haskey(v, "init")
                    value = v["init"] isa Quantity ? v["init"] : _parse_value(Float64, v["init"], UNIT)
                    if !iszero(value) 
                        if !haskey(cf, "offset"); cf["offset"] = zero(value) end
                        # propagate the units of the first offset value with units
                        if (value isa Quantity) && !(cf["offset"] isa Quantity) cf["offset"] = float(uconvert(unit(value), cf["offset"] * UNIT)) end
                        if !(value isa Quantity) && (cf["offset"] isa Quantity) value = float(uconvert(unit(cf["offset"]), value * UNIT)) end
                        cf["offset"] += value
                    end
                end
                # rearrange the gradient, but no need to parse/add anything, so just take the value
                if haskey(v, "gradient")
                    value = _parse_value(Float64, v["gradient"], UNIT / internal_length_unit)
                    if !iszero(value)
                        if !haskey(cf, "gradient"); cf["gradient"] = OrderedCollections.OrderedDict() end
                        cf["gradient"][k] = replace(string(v["gradient"]), " " => "*")
                    end
                end

            end
        end

        # manually handle the units added to the offset in the config file (if existent)
        merge!(cf, Dict("offset" => "$(replace(string(cf["offset"]), " " => "*"))"))

        # if there are no gradients: consider using a constant density
        if !haskey(cf, "gradient")
            cf["name"] = "constant"
            cf["value"] = replace(string(pop!(cf, "offset", 0)), " " => "*")
        end

        @warn msg * "Updating the config dictionary to the updated syntax:\n\n$(YAML.write(Dict("$(density)_density" => cf)))\n"
    end
    return cf
end

@inline update_impurity_density_config!(cf::AbstractDict, units::UnitTuple) = get(cf, "name", "") in ("linear", "cylindrical") ? _update_density_config!(cf, "impurity", units) : cf
@inline update_charge_density_config!(cf::AbstractDict, units::UnitTuple)   = get(cf, "name", "") in ("linear", "cylindrical") ? _update_density_config!(cf, "charge", units) : cf