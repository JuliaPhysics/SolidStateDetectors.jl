volume_primitive_dict = Dict{String, Any}(
    "tube" => Cone,
    "cone" => Cone,
    "sphere" => Sphere,
    "box" => Box,
    "HexagonalPrism" => HexagonalPrism,
    "union" => CSGUnion,
    "difference" => CSGDifference,
    "intersection" => CSGIntersection
)


#### INTERNAL PARSE FUNCTIONS

# parses dictionary entries of type Real or String to their value in internal units
@inline _parse_value(::Type{T}, x::Real, unit::Unitful.Units) where {T} = to_internal_units(T(x) * unit)
@inline _parse_value(::Type{T}, s::String, ::Unitful.Units) where {T} = to_internal_units(T(uparse(s)))

# parses dictionary entries of type {"from": ..., "to": ... } to a Tuple of the interval boundaries
function _parse_interval_from_to(::Type{T}, dict::Union{Dict{String,Any}, Dict{Any,Any}}, unit::Unitful.Units)::Tuple{T,T} where {T}
    To::T = _parse_value(T, dict["to"], unit)
    From::T = _parse_value(T, dict["from"], unit)
    @assert From <= To "Entry 'from' should be smaller than entry 'to' in $(dict)."
    (From, To)
end

# parses dictionary entries of type Real, String or {"from": ..., "to": ... } to respective AbstractFloat/Interval
@inline _parse_radial_interval(::Type{T}, x::Union{Real, String}, unit::Unitful.Units) where {T} = _parse_value(T, x, unit)
function _parse_radial_interval(::Type{T}, dict::Union{Dict{String,Any}, Dict{Any,Any}}, unit::Unitful.Units) where {T}
    @assert haskey(dict, "from") && haskey(dict, "to") "Please specify 'from' and 'to' in $(dict)."
    From::T, To::T = _parse_interval_from_to(T, dict, unit)
    @assert From >= 0 && To >= 0 "Entries 'from' and 'to' of radial $(dict) should be non-zero."
    From == 0 ? To : From..To
end

# parses dictionary entries of type Real, String or {"from": ..., "to": ... } to respective AbstractFloat/Interval
@inline _parse_linear_interval(::Type{T}, x::Union{Real, String}, unit::Unitful.Units) where {T} = _parse_value(T, x, unit)/2
function _parse_linear_interval(::Type{T}, dict::Union{Dict{String,Any}, Dict{Any,Any}}, unit::Unitful.Units) where {T}
    @assert haskey(dict, "from") && haskey(dict, "to") "Please specify 'from' and 'to' in $(dict)."
    From::T, To::T = _parse_interval_from_to(T, dict, unit)
    To == -From ? To : From..To
end


# parses dictionary entry for φ-interval that has {"from" ..., "to": ...} to the respective Nothing/Interval
function _parse_angular_interval(::Type{T}, dict::Union{Dict{String,Any}, Dict{Any,Any}}, unit::Unitful.Units) where {T}
    φTo::T = _parse_value(T, dict["to"], unit)
    φFrom::T = _parse_value(T, dict["from"], unit)
    if abs(rem2pi(φFrom - φTo, RoundNearest)) > 2*eps(T)
        φFrom..φTo
        #else nothing 
    end
end



### ADAPTED FOR PRIMITIVES (should throw Errors if something is not defined!)

# converts "r" to the respective AbstractFloat/Interval/Tuple for Cone
function parse_r_of_primitive(::Type{T}, dict::Union{Dict{String,Any}, Dict{Any,Any}}, unit::Unitful.Units) where {T}
    @assert haskey(dict, "r") "Please specify 'r' of the '$(dict["type"])'."
    dictr = dict["r"]
    # "r" : {"bottom": {"from": ..., "to": ...}, "top": {"from": ..., "to": ...}}
    if haskey(dictr, "bottom") && haskey(dictr, "top")
        bottom = _parse_radial_interval(T, dictr["bottom"], unit)
        top = _parse_radial_interval(T, dictr["top"], unit)
        # bottom and top need to be same type for Cone
        typeof(bottom) == typeof(top) ? (bottom, top) : _extend_number_to_zero_interval.((bottom, top))
    # "r" : {"from": ..., "to": ...}
    else
        _parse_radial_interval(T, dictr, unit)
    end
end

# converts "phi" : {...} to the respective Nothing/Interval
# if no φ is given, then a full 360° interval (i.e. φ = nothing) is assumed
function parse_φ_of_primitive(::Type{T}, dict::Union{Dict{String,Any}, Dict{Any,Any}}, unit::Unitful.Units) where {T}
    φ = if !haskey(dict,"phi")
        @info "No 'phi' is specified in the '$(dict["type"])'. Assuming 'phi' to go from 0 to 360°."
        nothing
    else #haskey(dict, "phi")
        _parse_angular_interval(T, dict["phi"], unit)
    end
end

# converts either "h" or "z" to the respective AbstractFloat/Interval
function parse_height_of_primitive(::Type{T}, dict::Union{Dict{String,Any}, Dict{Any,Any}}, unit::Unitful.Units) where {T}
    @assert haskey(dict,"h") || haskey(dict,"z") "Please specify 'h' or 'z' of the '$(dict["type"])'."
    haskey(dict,"h") ? _parse_linear_interval(T, dict["h"], unit) : _parse_linear_interval(T, dict["z"], unit)
end

# converts the content of the dictionary to respective AbstractFloat/Interval.
function parse_interval_of_primitive(::Type{T}, s::String, dict::Union{Dict{String,Any}, Dict{Any,Any}}, unit::Unitful.Units) where {T}
    @assert haskey(dict, s) "Please specify '$(s)' of the '$(dict["type"])'."
    _parse_linear_interval(T, dict[s], unit)
end


# converts {"x": ..., "y": ..., "z": ... } to the respective CartesianVector
function parse_translate_vector(::Type{T}, dict::Union{Dict{String,Any}, Dict{Any,Any}}, unit::Unitful.Units)::CartesianVector{T} where {T}
    x::T = haskey(dict, "x") ? _parse_value(T, dict["x"], unit) : T(0)
    y::T = haskey(dict, "y") ? _parse_value(T, dict["y"], unit) : T(0)
    z::T = haskey(dict, "z") ? _parse_value(T, dict["z"], unit) : T(0)
    CartesianVector{T}(x,y,z)
end


# parses a geometry (with possible translate vector) into the respective AbstractGeometry
function Geometry(::Type{T}, dict::Union{Dict{String,Any}, Dict{Any,Any}}, input_units::NamedTuple) where {T}
    if haskey(dict, "translate")
        length_unit = input_units.length
        t::CartesianVector{T} = parse_translate_vector(T, dict["translate"], length_unit)
        gdict::Dict{Any,Any} = filter(p -> first(p) != "translate", dict)
        return Geometry(T, gdict, input_units) + t
    end
    Geometry(T, volume_primitive_dict[dict["type"]], dict, input_units)
end

