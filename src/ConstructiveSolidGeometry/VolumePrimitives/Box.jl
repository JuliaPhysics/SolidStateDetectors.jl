struct Box{T,TX,TY,TZ} <: AbstractVolumePrimitive{T}
    x::TX
    y::TY
    z::TZ
    
    function Box( ::Type{T},
                  x::Union{T, <:AbstractInterval{T}},
                  y::Union{T, <:AbstractInterval{T}},
                  z::Union{T, <:AbstractInterval{T}}) where {T}
        new{T,typeof(x),typeof(y),typeof(z)}(x, y, z)
    end
end

#Constructors
function Box(;xMin = -1, xMax = 1, yMin = -1, yMax = 1, zMin = -1, zMax = 1)
    T = float(promote_type(typeof.((xMin, xMax, yMin, yMax, zMin, zMax))...))
    x = xMax == -xMin ? T(xMax) : T(xMin)..T(xMax)
    y = yMax == -yMin ? T(yMax) : T(yMin)..T(yMax)
    z = zMax == -zMin ? T(zMax) : T(zMin)..T(zMax)
    Box( T, x, y, z)
end
Box(xMin, xMax, yMin, yMax, zMin, zMax) = Box(; xMin = xMin, xMax = xMax, yMin = yMin, yMax = yMax, zMin = zMin, zMax = zMax)

function Box(x::X, y::Y, z::Z) where {X<:Real, Y<:Real, Z<:Real}
    T = float(promote_type(X,Y,Z))
    Box( T, T(x/2), T(y/2), T(z/2))
end

in(p::AbstractCoordinatePoint, b::Box) =
    _in_x(p, b.x) && _in_y(p, b.y) && _in_z(p, b.z)
    
    
# read-in
function Geometry(::Type{T}, ::Type{Box}, dict::AbstractDict, input_units::NamedTuple) where {T}
    length_unit = input_units.length
    x = parse_interval_of_primitive(T, "x", dict, length_unit)
    y = parse_interval_of_primitive(T, "y", dict, length_unit)
    z = parse_interval_of_primitive(T, "z", dict, length_unit)
    return Box(T, x, y, z)
end

function Dictionary(b::Box{T}) where {T}
    dict = OrderedDict{String,Any}()
    dict["x"] = typeof(b.x) == T ? b.x * 2 : OrderedDict{String,Any}("from" => b.x.left, "to" => b.x.right)
    dict["y"] = typeof(b.y) == T ? b.y * 2 : OrderedDict{String,Any}("from" => b.y.left, "to" => b.y.right)
    dict["z"] = typeof(b.z) == T ? b.z * 2 : OrderedDict{String,Any}("from" => b.z.left, "to" => b.z.right)
    OrderedDict{String,Any}("box" => dict)
end

get_x_limits(b::Box{T}) where {T} = (_left_linear_interval(b.x), _right_linear_interval(b.x))
get_y_limits(b::Box{T}) where {T} = (_left_linear_interval(b.y), _right_linear_interval(b.y))
get_z_limits(b::Box{T}) where {T} = (_left_linear_interval(b.z), _right_linear_interval(b.z))

function get_decomposed_surfaces(b::Box{T}) where {T}
    xMin::T, xMax::T = get_x_limits(b)
    yMin::T, yMax::T = get_y_limits(b)
    zMin::T, zMax::T = get_z_limits(b)
    tol = geom_atol_zero(T)
    if isapprox(xMin, xMax, atol = tol)
        return AbstractSurfacePrimitive[RectangleX(b, xMin)]
    elseif isapprox(yMin, yMax, atol = tol)
        return AbstractSurfacePrimitive[RectangleY(b, yMin)]
    elseif isapprox(zMin, zMax, atol = tol)
        return AbstractSurfacePrimitive[RectangleZ(b, zMin)]
    else
        return AbstractSurfacePrimitive[
                                            RectangleX(b, xMin), RectangleX(b, xMax), 
                                            RectangleY(b, yMin), RectangleY(b, yMax), 
                                            RectangleZ(b, zMin), RectangleZ(b, zMax)
                                        ]
    end
end

function sample(b::Box{T}, Nsamps::NTuple{3,Int} = (2,2,2)) where {T}
    xMin::T, xMax::T = get_x_limits(b)
    yMin::T, yMax::T = get_y_limits(b)
    zMin::T, zMax::T = get_z_limits(b)
    samples = [
        CartesianPoint{T}(x,y,z)
        for x in (Nsamps[1] ≤ 1 ? xMin : range(xMin, xMax, length = Nsamps[1]))
        for y in (Nsamps[2] ≤ 1 ? yMin : range(yMin, yMax, length = Nsamps[2]))
        for z in (Nsamps[3] ≤ 1 ? zMin : range(zMin, zMax, length = Nsamps[3]))
    ]
end
