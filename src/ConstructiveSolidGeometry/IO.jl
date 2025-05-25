struct ConfigFileError <: Exception
    msg::AbstractString
end
Base.showerror(io::IO, e::ConfigFileError) = print(io, "ConfigFileError: ", e.msg)


const CSG_dict = Dict{String, Any}(
    "tube" => Cone,
    "cone" => Cone,
    "sphere" => Ellipsoid,
    "box" => Box,
    "torus" => Torus,
    "polycone" => Polycone,
    "TriangularPrism" => TriangularPrism,
    "QuadranglePrism" => QuadranglePrism,
    "PentagonalPrism" => PentagonalPrism,
    "HexagonalPrism"  => HexagonalPrism,
    "union" => CSGUnion,
    "difference" => CSGDifference,
    "intersection" => CSGIntersection,
    "translate" => CartesianVector, # we just need some type to dispatch on
    "rotate" => Rotations.Rotation  # we just need some type to dispatch on
)

function get_geometry_key(::Type{T}, dict::AbstractDict, input_units::NamedTuple) where {T}
    g = collect(filter(k -> k in keys(CSG_dict), keys(dict)))
    filter!(k -> !(
        (k == "translate" && !any(broadcast(key -> key in keys(CSG_dict), keys(dict["translate"])))) ||
        (k == "rotate" && !any(broadcast(key -> key in keys(CSG_dict), keys(dict["rotate"]))))
    )  , g) # exclude translate and rotate that have no primitives inside
    @assert !(length(g) > 1) "Too many ($(length(g))) geometry entries in dictionary: $(g)."
    @assert !(length(g) == 0) "None of the entries $(keys(dict)) describes a Geometry."
    g[1]
end

# This function will only be called when constructing Semiconductor, Contacts or Passives
# and if "translate" and "rotate" are given OUTSIDE of the "geometry" entry
function parse_CSG_transformation(::Type{T}, dict::AbstractDict, input_units::NamedTuple) where {T}
    r = if haskey(dict, "rotate")
        SMatrix{3, 3, T, 9}(parse_rotation_matrix(T, dict["rotate"], input_units.angle))
        # We might want to improve this to avoid numerical precision errors. 
        # For now this is fine. 
    else
        one(SMatrix{3, 3, T, 9})
    end
    t = if haskey(dict, "translate")
        parse_translate_vector(T, dict["translate"], input_units.length)
    else
        zero(CartesianVector{T})
    end
    (rotation = r, translation = t)
end

#### INTERNAL PARSE FUNCTIONS

# parses dictionary entries of type Real or String to their value in internal units
@inline _parse_value(::Type{T}, x::Real, unit::Unitful.Units) where {T} = T(to_internal_units(float(x) * unit))
@inline _parse_value(::Type{T}, x::Quantity, unit::Unitful.Units) where {T} = begin
    if dimension(x) != dimension(unit) throw(ConfigFileError("The quantity '$(x)' cannot be converted to unit '$(unit)'")) end
    T(to_internal_units(float(x)))
end
@inline _parse_value(::Type{T}, s::String, unit::Unitful.Units) where {T} = _parse_value(T, uparse(s), unit)
@inline _parse_value(::Type{T}, a::Union{<:AbstractVector, <:Tuple}, unit::Unitful.Units) where {T} = _parse_value.(T, a, unit)

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
    @assert From >= 0 && To >= 0 "Entries 'from' and 'to' of radial $(dict) should `>= 0`."
    (From, To)
end

# parses dictionary entries of type Real, String or {"from": ..., "to": ... } to respective AbstractFloat/Interval
@inline _parse_linear_interval(::Type{T}, x::Union{Real, String}, unit::Unitful.Units) where {T} = _parse_value(T, x, unit)/2
function _parse_linear_interval(::Type{T}, dict::AbstractDict, unit::Unitful.Units) where {T}
    @assert haskey(dict, "from") && haskey(dict, "to") "Please specify 'from' and 'to' in $(dict)."
    From::T, To::T = _parse_interval_from_to(T, dict, unit)
    To == -From == zero(T) ? To : (From, To) # if != 0 is influences the origin 
end

function _parse_linear_interval(::Type{T}, tuple::Tuple{Real,Real}, unit::Unitful.Units) where {T}
    tuple[2] == -tuple[1] == zero(T) ? tuple[2] : (tuple[1], tuple[2]) # if != 0 is influences the origin 
end

# parses dictionary entry for φ-interval that has {"from" ..., "to": ...} to the respective Nothing/Interval
function _parse_angular_interval(::Type{T}, dict::AbstractDict, unit::Unitful.Units) where {T}
    φTo::T = _parse_value(T, dict["to"], unit)
    φFrom::T = _parse_value(T, dict["from"], unit)
    if abs(rem2pi(φFrom - φTo, RoundNearest)) > 2*eps(T)
        φFrom, φTo
    else 
        nothing 
    end
end



### ADAPTED FOR PRIMITIVES (should throw Errors if something is not defined!)

# converts "r" to the respective AbstractFloat/Interval/Tuple for Cone
function parse_r_of_primitive(::Type{T}, dict::AbstractDict, unit::Unitful.Units) where {T}
    @assert haskey(dict, "r") "Please specify 'r'."
    dictr = dict["r"]
    if haskey(dictr, "bottom") && haskey(dictr, "top")
        # "r" : {"bottom": {"from": ..., "to": ...}, "top": {"from": ..., "to": ...}}
        bottom = _parse_radial_interval(T, dictr["bottom"], unit)
        top = _parse_radial_interval(T, dictr["top"], unit)
        bottom, top
    else
        # "r" : {"from": ... , "to": ...}
        _parse_radial_interval(T, dictr, unit)
    end
end
# converts "r" to the respective AbstractFloat/Interval/Tuple for Cone
function parse_r_of_primitive(::Type{T}, dict::AbstractDict, unit::Unitful.Units, ::Type{Cone}) where {T}
    @assert haskey(dict, "r") "Please specify 'r'."
    dictr = dict["r"]
    r = if haskey(dictr, "bottom") && haskey(dictr, "top")
        # "r" : {"bottom": {"from": ..., "to": ...}, "top": {"from": ..., "to": ...}}
        bottom = _parse_radial_interval(T, dictr["bottom"], unit)
        top = _parse_radial_interval(T, dictr["top"], unit)
        if bottom isa Real bottom = (zero(T), bottom) end
        if top isa Real top = (zero(T), top) end
        (bottom, top)
    elseif haskey(dictr, "from") && haskey(dictr, "to")
        # "r" : {"from": ... , "to": ...}
        r = _parse_radial_interval(T, dictr, unit)
        if r isa Real r = (zero(T), r) end
        (r, r)
    else
        r = _parse_radial_interval(T, dictr, unit)
        ((zero(T),r),(zero(T),r))
    end
    # Not all cases implemented yet 
    r_bot_in, r_bot_out, r_top_in, r_top_out = r[1][1], r[1][2], r[2][1], r[2][2]
    if r_bot_in == r_top_in == zero(T) 
        if r_bot_out == r_top_out
            return r_bot_out # Cylinder
        else
            return ((r_bot_out,), (r_top_out,)) # VaryingCylinder
        end
    elseif r_top_in == r_top_out && r_bot_in != r_bot_out
        if r_top_in == 0
            return ((r_bot_in, r_bot_out), nothing)
        else
            return r # return ((r_bot_in, r_bot_out), r_top_in)
        end
    elseif r_top_in != r_top_out && r_bot_in == r_bot_out
        if r_bot_in == 0
            return (nothing, (r_top_in, r_top_out))
        else
            return r # return (r_bot_in, (r_top_in, r_top_out))
        end
    else
        r # VaryingTube
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

function Geometry(::Type{T}, ::Type{CartesianVector}, dict::AbstractDict, input_units::NamedTuple, outer_transformations::Transformations{T}) where {T}
    translate_vector = parse_translate_vector(T, dict, input_units.length)
    key = get_geometry_key(T, dict, input_units)
    inner_transformations = (rotation = one(SMatrix{3, 3, T, 9}), translation = translate_vector)
    transformations = combine_transformations(inner_transformations, outer_transformations)
    Geometry(T, CSG_dict[key], dict[key], input_units, transformations)
end

# function parse_scale_vector(::Type{T}, dict::AbstractDict)::SVector{3,T} where {T}
#     x::T = haskey(dict, "x") ? _parse_value(T, dict["x"], Unitful.FreeUnits{(), NoDims, nothing}()) : T(1)
#     y::T = haskey(dict, "y") ? _parse_value(T, dict["y"], Unitful.FreeUnits{(), NoDims, nothing}()) : T(1)
#     z::T = haskey(dict, "z") ? _parse_value(T, dict["z"], Unitful.FreeUnits{(), NoDims, nothing}()) : T(1)
#     SVector{3,T}(x,y,z)
# end

function get_origin(::Type{T}, dict::AbstractDict, unit::Unitful.Units) where {T}
    origin_names = ["origin", "translate", "translation"]
    haskeys = map(n -> haskey(dict, n), origin_names)
    i = findfirst(i -> i, haskeys)
    haskeyssum = sum(haskeys)
    if haskeyssum == 1
        o = dict[origin_names[i]]
        if o isa AbstractVector
            CartesianPoint{T}(_parse_value(T, o, unit))
        else
            x, y, z = zero(T), zero(T), zero(T)
            if haskey(o, "x") x = _parse_value(T, o["x"], unit) end
            if haskey(o, "y") y = _parse_value(T, o["y"], unit) end
            if haskey(o, "z") z = _parse_value(T, o["z"], unit) end
            CartesianPoint{T}(x, y, z)
        end
    elseif haskeyssum > 1
        throw(ConfigFileError("Multiple fields for origin of primitive detected. Use only `origin`, `translate` or `translation`."))
    else
        zero(CartesianPoint{T})
    end
end


# function Geometry(::Type{T}, ::Type{ScaledGeometry}, dict::AbstractDict, input_units::NamedTuple) where {T}
#     length_unit = input_units.length
#     scale_vector::SVector{3,T} = parse_scale_vector(T, dict, length_unit)
#     key::String, transformations = get_geometry_key(T, dict, input_units)
#     scale(transform(Geometry(T, CSG_dict[key], dict[key], input_units), transformations), scale_vector)
# end

_rot_keys = ["X","Y","Z","XY","XZ","YX","YZ","ZX","ZY","XYX","XYZ","XZX","XZY","YXY","YXZ","YZX","YZY","ZXY","ZXZ","ZYX","ZYZ"]
function parse_rotation_matrix(::Type{T}, dict::AbstractDict, unit::Unitful.Units)::RotMatrix3{T} where {T}
    dict_keys = keys(dict)
    if "M" in dict_keys
        R = RotMatrix3{T}(dict["M"])
        @assert isrotation(R) "R is not a Rotation Matrix."
        return R
    else 
        rotation = filter(k -> uppercase(k) in _rot_keys, dict_keys)
        @assert length(rotation) == 1 "Rotations must be defined using exactly one of the following keywords:\n$(vcat("M",_rot_keys))"
        key = first(rotation)
        R = getfield(Rotations, Symbol("Rot"*uppercase(key)))
        RotMatrix3{T}(R{Float64}(_parse_value(Float64,dict[key],unit)...)) # parse rotations matrix values as Float64 to minimize rounding errors
    end
end

function get_rotation(::Type{T}, dict::AbstractDict, unit::Unitful.Units) where {T}
    rot_names = ["rotate", "rotation"]
    haskeys = map(n -> haskey(dict, n), rot_names)
    i = findfirst(i -> i, haskeys)
    haskeyssum = sum(haskeys)
    if haskeyssum == 1
        parse_rotation_matrix(T, dict[rot_names[i]], unit)
    elseif haskeyssum > 1
        throw(ConfigFileError("Multiple fields for rotation of primitive detected. Use only `rotation` or `rotate`."))
    else
        one(SMatrix{3, 3, T, 9})
    end
end

# function Geometry(::Type{T}, ::Type{RotatedGeometry}, dict::AbstractDict, input_units::NamedTuple) where {T}
#     angle_unit = input_units.angle
#     rotation_matrix = parse_rotation_matrix(T, dict, angle_unit)
#     key::String, transformations = get_geometry_key(T, dict, input_units)
#     rotate(transform(Geometry(T, CSG_dict[key], dict[key], input_units), transformations), rotation_matrix)
# end

function Geometry(::Type{T}, ::Type{Rotations.Rotation}, dict::AbstractDict, input_units::NamedTuple, outer_transformations::Transformations{T}) where {T}
    rotation_matrix = SMatrix{3, 3, T, 9}(parse_rotation_matrix(T, dict, input_units.angle))
    key = get_geometry_key(T, dict, input_units)
    inner_transformations = (rotation = rotation_matrix, translation = CartesianVector{T}(0,0,0))
    transformations = combine_transformations(inner_transformations, outer_transformations)
    Geometry(T, CSG_dict[key], dict[key], input_units, transformations)
end


function Geometry(::Type{T}, dict::AbstractDict, input_units::NamedTuple, outer_transformations::Transformations{T}) where {T}
    key = get_geometry_key(T, dict, input_units)
    Geometry(T, CSG_dict[key], dict[key], input_units, outer_transformations)
end

function Geometry(::Type{T}, filename::String, input_units::NamedTuple = (length = internal_length_unit, angle = internal_angle_unit)) where {T}
    @assert isfile(filename) "The given filename '$(filename)' does not lead to a valid file."
    dict = if endswith(filename, ".json")
        JSON.parsefile(filename)
    elseif endswith(filename, ".yaml")
        YAML.load_file(filename)
    else
        @error "Only JSON and YAML formats are supported at the moment."
    end
    transformation = (rotation = one(SMatrix{3, 3, T, 9}), translation = zero(CartesianVector{T}))
    Geometry(T, dict, input_units, transformation)
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
    # if CT <: AbstractTransformedGeometry
    #     println(start*" $(CT.name.name){$(CT.parameters[1])}")
    #     tab *= "\t"
    #     show_CSG_tree(csgtree.p, start = tab*"└──", tab = tab)
    # end
    if CT <: AbstractPrimitive
        println(start*" $(CT.name.name){$(CT.parameters[1])}")
    end
end
