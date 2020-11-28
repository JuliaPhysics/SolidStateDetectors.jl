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
    T = promote_type(typeof.((xMin, xMax, yMin, yMax, zMin, zMax))...)
    x = xMax == -xMin ? T(xMax) : T(xMin)..T(xMax)
    y = yMax == -yMin ? T(yMax) : T(yMin)..T(yMax)
    z = zMax == -zMin ? T(zMax) : T(zMin)..T(zMax)
    Box( T, x, y, z)
end
Box(xMin, xMax, yMin, yMax, zMin, zMax) = Box(;xMin, xMax, yMin, yMax, zMin, zMax)

function Box(x::X, y::Y, z::Z) where {X,Y,Z}
    T = promote_type(X,Y,Z)
    Box( T, T(x/2), T(y/2), T(z/2))
end

in(p::AbstractCoordinatePoint, b::Box{<:Any, <:Any, <:Any, <:Any}) =
    _in_x(p, b.x) && _in_y(p, b.y) && _in_z(p, b.z)

    