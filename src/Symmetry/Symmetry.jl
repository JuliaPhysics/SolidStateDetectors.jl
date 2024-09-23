include("Mirror_Symmetry.jl")

function Symmetry(::Type{T}, symmetry_type::AbstractString, config_dict::AbstractDict, units::NamedTuple) where T
    if symmetry_type == "mirror"
        axis = collect(keys(config_dict))
        @assert length(axis) == 1 "Mirror Symmetry needs only one axis. $(length(axis)) axes were given: $(axis)"
        return MirrorSymmetry(axis[1],T(config_dict[axis[1]]), units)
    end
end

function Symmetry(::Type{T}, symmetry_types::AbstractArray, config_dicts::AbstractArray, units::NamedTuple) where T
    amt_symmetries = length(symmetry_types)
    @assert amt_symmetries == 1 "More than one symmetry per potential not yet implemented. 
        The following symmetries were specified $(symmetry_types)"
    Symmetry(T, symmetry_types[1], config_dicts[1], units)
end

function Symmetry(::Type{T}, config_file_dict::AbstractDict, units::NamedTuple)::NamedTuple where T
    nt = NamedTuple()
    if haskey(config_file_dict, "symmetry")
        sym_keys = collect(keys(config_file_dict["symmetry"])) # Vector with all components of specified symmetry (e.g. contact_1, electric_potential)
        sym_values = collect(values(config_file_dict["symmetry"]))
        sym_values = Tuple(Symmetry(T, collect(keys(sym_values[i])), collect(values(sym_values[i])), units) for i in 1:length(sym_values))
        sym_keys = Tuple(Symbol(sym_keys[i]) for i in 1:length(sym_keys))
        nt = NamedTuple{sym_keys}(sym_values)
    end
    return nt
end

    
 
