const CSG_dict = Dict{String, Any}(
    "tube" => Cone,
    "cone" => Cone,
    "sphere" => Sphere,
    "box" => Box,
    "torus" => Torus,
    "TriangularPrism" => TriangularPrism,
    "SquarePrism"  => SquarePrism,
    "PentagonalPrism" => PentagonalPrism,
    "HexagonalPrism"  => HexagonalPrism,
    "union" => CSGUnion,
    "difference" => CSGDifference,
    "intersection" => CSGIntersection,
    "translate" => TranslatedGeometry,
    "rotate" => RotatedGeometry,
    "scale" => ScaledGeometry
)

function get_geometry_key(::Type{T}, dict::AbstractDict, input_units::NamedTuple)::Tuple{String, Vector{CSGTransformation}} where {T}
    dict_keys = filter(k -> k in keys(CSG_dict), keys(dict))
    transformations = sort!(filter(k -> 
                                (k == "translate" && !any(broadcast(key -> key in keys(CSG_dict), keys(dict["translate"])))) || 
                                (k == "rotate" && !any(broadcast(key -> key in keys(CSG_dict), keys(dict["rotate"])))), 
                            collect(dict_keys)), 
                      rev = true) #this will ensure that the rotation is parsed first and the translation is parsed second
    primitives = setdiff(dict_keys, transformations)
    @assert length(primitives) <= 1 "Too many geometry entries in dictionary: $(length(dict_keys))."
    @assert length(primitives) >= 1 "None of the entries $(keys(dict)) describes a Geometry."
    first(primitives), broadcast(t -> parse_CSG_transformation(T, dict, CSG_dict[t], input_units), transformations)
end


#### INTERNAL PARSE FUNCTIONS

# parses dictionary entries of type Real or String to their value in internal units
@inline _parse_value(::Type{T}, x::Real, unit::Unitful.Units) where {T} = to_internal_units(T(x) * unit)
@inline _parse_value(::Type{T}, s::String, ::Unitful.Units) where {T} = to_internal_units(T(uparse(s)))
@inline _parse_value(::Type{T}, a::Vector, unit::Unitful.Units) where {T} = _parse_value.(T, a, unit)

# parses dictionary entries of type {"from": ..., "to": ... } to a Tuple of the interval boundaries
function _parse_interval_from_to(::Type{T}, dict::AbstractDict, unit::Unitful.Units)::Tuple{T,T} where {T}
    To::T = _parse_value(T, dict["to"], unit)
    From::T = _parse_value(T, dict["from"], unit)
    @assert From <= To "Entry 'from' should be smaller than entry 'to' in $(dict)."
    (From, To)
end

# parses dictionary entries of type Real, String or {"from": ..., "to": ... } to respective AbstractFloat/Interval
@inline _parse_radial_interval(::Type{T}, x::Union{Real, String}, unit::Unitful.Units) where {T} = _parse_value(T, x, unit)
function _parse_radial_interval(::Type{T}, dict::AbstractDict, unit::Unitful.Units) where {T}
    @assert haskey(dict, "from") && haskey(dict, "to") "Please specify 'from' and 'to' in $(dict)."
    From::T, To::T = _parse_interval_from_to(T, dict, unit)
    @assert From >= 0 && To >= 0 "Entries 'from' and 'to' of radial $(dict) should be non-zero."
    From == 0 ? To : From..To
end

# parses dictionary entries of type Real, String or {"from": ..., "to": ... } to respective AbstractFloat/Interval
@inline _parse_linear_interval(::Type{T}, x::Union{Real, String}, unit::Unitful.Units) where {T} = _parse_value(T, x, unit)/2
function _parse_linear_interval(::Type{T}, dict::AbstractDict, unit::Unitful.Units) where {T}
    @assert haskey(dict, "from") && haskey(dict, "to") "Please specify 'from' and 'to' in $(dict)."
    From::T, To::T = _parse_interval_from_to(T, dict, unit)
    To == -From ? To : From..To
end


# parses dictionary entry for φ-interval that has {"from" ..., "to": ...} to the respective Nothing/Interval
function _parse_angular_interval(::Type{T}, dict::AbstractDict, unit::Unitful.Units) where {T}
    φTo::T = _parse_value(T, dict["to"], unit)
    φFrom::T = _parse_value(T, dict["from"], unit)
    if abs(rem2pi(φFrom - φTo, RoundNearest)) > 2*eps(T)
        φFrom..φTo
        #else nothing 
    end
end



### ADAPTED FOR PRIMITIVES (should throw Errors if something is not defined!)

# converts "r" to the respective AbstractFloat/Interval/Tuple for Cone
function parse_r_of_primitive(::Type{T}, dict::AbstractDict, unit::Unitful.Units) where {T}
    @assert haskey(dict, "r") "Please specify 'r'."
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
function parse_φ_of_primitive(::Type{T}, dict::AbstractDict, unit::Unitful.Units) where {T}
    φ = if !haskey(dict,"phi")
        #@info "No 'phi' is specified. Assuming 'phi' to go from 0 to 360°."
        nothing
    else #haskey(dict, "phi")
        _parse_angular_interval(T, dict["phi"], unit)
    end
end

function parse_θ_of_primitive(::Type{T}, dict::AbstractDict, unit::Unitful.Units) where {T}
    θ = if !haskey(dict,"theta")
        #@info "No 'theta' is specified. Assuming 'theta' to go from 0 to 360°."
        nothing
    else #haskey(dict, "theta")
        _parse_angular_interval(T, dict["theta"], unit)
    end
end

# converts either "h" or "z" to the respective AbstractFloat/Interval
function parse_height_of_primitive(::Type{T}, dict::AbstractDict, unit::Unitful.Units) where {T}
    @assert haskey(dict,"h") || haskey(dict,"z") "Please specify 'h' or 'z'."
    haskey(dict,"h") ? _parse_linear_interval(T, dict["h"], unit) : _parse_linear_interval(T, dict["z"], unit)
end

# converts the content of the dictionary to respective AbstractFloat/Interval.
function parse_interval_of_primitive(::Type{T}, s::String, dict::AbstractDict, unit::Unitful.Units) where {T}
    @assert haskey(dict, s) "Please specify '$(s)'."
    _parse_linear_interval(T, dict[s], unit)
end


# converts {"x": ..., "y": ..., "z": ... } to the respective CartesianVector
function parse_translate_vector(::Type{T}, dict::AbstractDict, unit::Unitful.Units)::CartesianVector{T} where {T}
    x::T = haskey(dict, "x") ? _parse_value(T, dict["x"], unit) : T(0)
    y::T = haskey(dict, "y") ? _parse_value(T, dict["y"], unit) : T(0)
    z::T = haskey(dict, "z") ? _parse_value(T, dict["z"], unit) : T(0)
    CartesianVector{T}(x,y,z)
end

function parse_CSG_transformation(::Type{T}, dict::AbstractDict, ::Type{TranslatedGeometry}, input_units::NamedTuple)::CSGTransformation where {T}
    parse_translate_vector(T, dict["translate"], input_units.length)
end

function Geometry(::Type{T}, ::Type{TranslatedGeometry}, dict::AbstractDict, input_units::NamedTuple) where {T}
    length_unit = input_units.length
    translate_vector::CartesianVector{T} = parse_translate_vector(T, dict, length_unit)
    key::String, transformations::Vector{CSGTransformation} = get_geometry_key(T, dict, input_units)
    translate(transform(Geometry(T, CSG_dict[key], dict[key], input_units), transformations), translate_vector)
end

function parse_scale_vector(::Type{T}, dict::AbstractDict, unit::Unitful.Units)::SVector{3,T} where {T}
    x::T = haskey(dict, "x") ? _parse_value(T, dict["x"], unit) : T(1)
    y::T = haskey(dict, "y") ? _parse_value(T, dict["y"], unit) : T(1)
    z::T = haskey(dict, "z") ? _parse_value(T, dict["z"], unit) : T(1)
    SVector{3,T}(x,y,z)
end

function Geometry(::Type{T}, ::Type{ScaledGeometry}, dict::AbstractDict, input_units::NamedTuple) where {T}
    length_unit = input_units.length
    scale_vector::SVector{3,T} = parse_scale_vector(T, dict, length_unit)
    key::String, transformations::Vector{CSGTransformation} = get_geometry_key(T, dict, input_units)
    scale(transform(Geometry(T, CSG_dict[key], dict[key], input_units), transformations), scale_vector)
end

_rot_keys = ["X","Y","Z","XY","XZ","YX","YZ","ZX","ZY","XYX","XYZ","XZX","XZY","YXY","YXZ","YZX","YZY","ZXY","ZXZ","ZYX","ZYZ"]
function parse_rotation_matrix(::Type{T}, dict::AbstractDict, unit::Unitful.Units)::RotMatrix3{T} where {T}
    dict_keys = keys(dict)
    if "M" in dict_keys
        R = RotMatrix3{T}(dict["M"])
        @assert isrotation(R) "R is not a Rotation Matrix."
        return R
    else 
        rotation = filter(k -> k in _rot_keys, uppercase.(dict_keys))
        @assert length(rotation) == 1 "Rotations must be defined using exactly one of the following keywords:\n$(vcat("M",_rot_keys))"
        key = first(rotation)
        R = getfield(Rotations, Symbol("Rot"*key))
        RotMatrix3{T}(R{Float64}(_parse_value(Float64,dict[key],unit)...)) # parse rotations matrix values as Float64 to minimize rounding errors
    end
end

function parse_CSG_transformation(::Type{T}, dict::AbstractDict, ::Type{RotatedGeometry}, input_units::NamedTuple)::CSGTransformation where {T}
    parse_rotation_matrix(T, dict["rotate"], input_units.angle)
end

function Geometry(::Type{T}, ::Type{RotatedGeometry}, dict::AbstractDict, input_units::NamedTuple) where {T}
    angle_unit = input_units.angle
    rotation_matrix = parse_rotation_matrix(T, dict, angle_unit)
    key::String, transformations::Vector{CSGTransformation} = get_geometry_key(T, dict, input_units)
    rotate(transform(Geometry(T, CSG_dict[key], dict[key], input_units), transformations), rotation_matrix)
end


function Geometry(::Type{T}, dict::AbstractDict, input_units::NamedTuple) where {T}
    key::String, transformations::Vector{CSGTransformation} = get_geometry_key(T, dict, input_units)
    transform(Geometry(T, CSG_dict[key], dict[key], input_units), transformations)
end


function show_CSG_tree(csgtree; start = "", tab = "", CSG = false)
    CT = typeof(csgtree)
    if CT <: AbstractConstructiveGeometry
        if !CSG 
            println(start*" $(CT.name.name){$(CT.parameters[1])}")
            tab *= "\t"
        end
        next_CSG::Bool = typeof(csgtree.a).name == CT.name
        show_CSG_tree(csgtree.a, start = tab*"├"*(CT <: CSGDifference ? " +" : "──"), tab = tab*(next_CSG ? "" : "│"), CSG = next_CSG)
        show_CSG_tree(csgtree.b, start = tab*(CSG ? "├" : "└")*(CT <: CSGDifference ? " ─" : "──"), tab = tab*(CSG ? "│" : " "))
    end
    if CT <: AbstractTransformedGeometry
        println(start*" $(CT.name.name){$(CT.parameters[1])}")
        tab *= "\t"
        show_CSG_tree(csgtree.p, start = tab*"└──", tab = tab)
    end
    if CT <: AbstractPrimitive
        println(start*" $(CT.name.name){$(CT.parameters[1])}")
    end
end
