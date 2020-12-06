# converts "r" : {"bottom": {"from": ..., "to": ...}, "top": {"from": ..., "to": ...}} to the respective AbstractFloat/Interval/Tuple for Cone
function _get_r_of_primitive(T::DataType, dictr::Union{Dict{String,Any}, Dict{Any,Any}}, length_unit::Unitful.Units, t::Val{:cone})
    dictrbot::Dict{Any,Any} = dictr["bottom"]
    rbot = _get_r_of_primitive(T, dictrbot, length_unit)
    dictrtop::Dict{Any,Any} = dictr["top"]
    rtop = _get_r_of_primitive(T, dictrtop, length_unit)
    rbot, rtop
end

# converts "r" :  {"from": ... , "to": ...} to the respective AbstractFloat/Interval
function _get_r_of_primitive(T::DataType, dictr::Union{Dict{String,Any}, Dict{Any,Any}}, length_unit::Unitful.Units, t::Val{:tube} = Val{:tube}())
    rTo::T = geom_round(T(to_internal_units(internal_length_unit, dictr["to"] * length_unit)))
    if haskey(dictr, "from")
        rFrom::T = geom_round(T(to_internal_units(internal_length_unit, dictr["from"] * length_unit)))
        if rFrom == 0
            rTo
        else
            rFrom..rTo
        end
    else
        rTo
    end
end

# converts "r" : ... to the respective AbstractFloat
function _get_r_of_primitive(T::DataType, dictr::Real, length_unit::Unitful.Units, t::Val{:tube} = Val{:tube}())
    T(dictr)
end

# converts "phi" : {"from": ..., "to": ...} to the respective Nothing/Interval
# if no φ is given, then a full 360° interval (i.e. φ = nothing) is assumed
function _get_φ_of_primitive(T::DataType, dict::Union{Dict{String,Any}, Dict{Any,Any}}, angle_unit::Unitful.Units)
    φ = if !haskey(dict,"phi")
        @info "No 'phi' is specified in the '$(dict["type"])'. Assuming 'phi' to go from 0 to 360°."
        nothing
    else #haskey(dict, "phi")
        dictφ::Dict{Any,Any} = dict["phi"]
        φFrom::T = geom_round(to_internal_units(internal_angle_unit, T(dictφ["from"]) * angle_unit))
        φTo::T = geom_round(to_internal_units(internal_angle_unit, T(dictφ["to"]) * angle_unit))
        if abs(rem2pi(φFrom - φTo, RoundNearest)) > 2*eps(T)
            φFrom..φTo
            #else nothing 
        end
    end
end

# converts either "h" or "z" to the respective AbstractFloat/Interval
function _get_h_or_z_of_primitive(T::DataType, dict::Union{Dict{String,Any}, Dict{Any,Any}}, length_unit::Unitful.Units)
    @assert haskey(dict,"h") || haskey(dict,"z") "Please specify 'h' or 'z' of the '$(dict["type"])'."
    z = if haskey(dict,"h")
        _get_h_of_primitive(T, dict, length_unit)
    else #haskey(dict, "z")
        _get_linear_interval_of_primitive(T, "z", dict, length_unit)
    end
end

# converts "h" : ... to the respective AbstractFloat
function _get_h_of_primitive(T::DataType, dict::Union{Dict{String,Any}, Dict{Any,Any}}, length_unit::Unitful.Units)
    @assert haskey(dict,"h") "Please specify 'h' of the '$(dict["type"])'."
    h::T = geom_round(T(to_internal_units(internal_length_unit, dict["h"] * length_unit)))
    h/2
end


# converts any linear interval "...": {"from": ..., "to":... } to the respective AbstractFloat/Interval
function _get_linear_interval_of_primitive(T::DataType, s::String, dict::Union{Dict{String,Any}, Dict{Any,Any}}, length_unit::Unitful.Units)
    @assert haskey(dict, s) "Please specify '$s' of the '$(dict["type"])'."
    dicts::Dict{Any,Any} = dict[s]
    From::T = geom_round(T(to_internal_units(internal_length_unit, dicts["from"] * length_unit)))
    To::T = geom_round(T(to_internal_units(internal_length_unit, dicts["to"] * length_unit)))
    if To == -From 
        To 
    else
        From..To
    end
end

# converts "translate" : { ... } to the respective CartesianVector 
function _get_translate_vector(T::DataType, dict::Dict{Any,Any}, iud::Dict{String, Unitful.Units})::CartesianVector{T}
    conversion_factor::T = ustrip(uconvert(internal_length_unit, one(T) * iud["length"]))
    x::T = haskey(dict, "x") ? geom_round(conversion_factor * T(dict["x"])) : T(0)
    y::T = haskey(dict, "y") ? geom_round(conversion_factor * T(dict["y"])) : T(0)
    z::T = haskey(dict, "z") ? geom_round(conversion_factor * T(dict["z"])) : T(0)
    CartesianVector{T}(x,y,z)
end


# parses a geometry (with possible translate vector) into the respective AbstractGeometry
function Geometry(T::DataType, dict::Union{Dict{String,Any}, Dict{Any,Any}}, iud::Dict{String,Unitful.Units})
    if haskey(dict, "translate")
        tdict::Dict{Any,Any} = dict["translate"] 
        t::CartesianVector{T} = _get_translate_vector(T, tdict, iud)
        gdict::Dict{Any,Any} = filter(p -> first(p) != "translate", dict)
        return t == CartesianVector{T}(0,0,0) ? Geometry(T, gdict, iud) : translate(Geometry(T, gdict, iud), t)
    end
    Geometry(T, Val{Symbol(dict["type"])}(), dict, iud)
end

