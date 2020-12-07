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
Box(xMin, xMax, yMin, yMax, zMin, zMax) = Box(;xMin, xMax, yMin, yMax, zMin, zMax)

function Box(x::X, y::Y, z::Z) where {X<:Real, Y<:Real, Z<:Real}
    T = float(promote_type(X,Y,Z))
    Box( T, T(x/2), T(y/2), T(z/2))
end

in(p::AbstractCoordinatePoint, b::Box{<:Any, <:Any, <:Any, <:Any}) =
    _in_x(p, b.x) && _in_y(p, b.y) && _in_z(p, b.z)
    
    
# read-in
function Geometry(::Type{T}, ::Type{Box}, dict::Union{Dict{String,Any}, Dict{Any,Any}}, input_units::NamedTuple) where {T}
    length_unit = input_units.length
    x = parse_interval_of_primitive(T, "x", dict, length_unit)
    y = parse_interval_of_primitive(T, "y", dict, length_unit)
    z = parse_interval_of_primitive(T, "z", dict, length_unit)
    return Box(T, x, y, z)
end
        
    
# plotting
function get_plot_points(b::Box{T}) where {T <: AbstractFloat}
    
    plot_points = Vector{CartesianPoint{T}}[]
    
    xMin::T = _left_linear_interval(b.x)
    xMax::T = _right_linear_interval(b.x)
    yMin::T = _left_linear_interval(b.y)
    yMax::T = _right_linear_interval(b.y)
    zMin::T = _left_linear_interval(b.z)
    zMax::T = _right_linear_interval(b.z)
    
    #horizontal squares
    for z in [zMin, zMax]
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(xMin, yMin, z), CartesianPoint{T}(xMin, yMax, z),
                            CartesianPoint{T}(xMax, yMax, z), CartesianPoint{T}(xMax, yMin, z), 
                            CartesianPoint{T}(xMin, yMin, z)]))
    end
    
    #vertical lines
    for x in [xMin, xMax]
        for y in [yMin, yMax]
            push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(x,y,zMin), CartesianPoint{T}(x,y,zMax)]))
        end
    end
        
    plot_points
end

    