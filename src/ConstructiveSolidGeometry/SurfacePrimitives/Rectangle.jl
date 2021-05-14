struct Rectangle{TN,T,TL,TW} <: AbstractSurfacePrimitive{T}
    l::TL
    w::TW
    loc::T
    function Rectangle( normal::Union{Val{:x}, Val{:y}, Val{:z}},
                        ::Type{T},
                        l::Union{T, <:AbstractInterval{T}},
                        w::Union{T, <:AbstractInterval{T}},
                        loc::T) where {T}
        new{typeof(normal),T,typeof(l),typeof(w)}(l,w,loc)
    end
end

# Convenience functions
const RectangleX{T,TL,TW} = Rectangle{Val{:x},T,TL,TW}
const RectangleY{T,TL,TW}  = Rectangle{Val{:y},T,TL,TW}
const RectangleZ{T,TL,TW} = Rectangle{Val{:z},T,TL,TW}

RectangleX(args...) = Rectangle(Val(:x), args...)
RectangleY(args...)  = Rectangle(Val(:y), args...)
RectangleZ(args...) = Rectangle(Val(:z), args...)

#Constructors
RectangleX(b::Box{T}, x::Real) where {T} = RectangleX(T, b.y, b.z, T(x))
RectangleY(b::Box{T}, y::Real) where {T} = RectangleY(T, b.x, b.z, T(y))
RectangleZ(b::Box{T}, z::Real) where {T} = RectangleZ(T, b.x, b.y, T(z))

function Rectangle(;normal = Val(:z), lMin = -1, lMax = 1, wMin = -1, wMax = 1, loc = 0)
    T = float(promote_type(typeof.((lMin, lMax, wMin, wMax, loc))...))
    l = lMax == -lMin ? T(lMax) : T(lMin)..T(lMax)
    w = wMax == -wMin ? T(wMax) : T(wMin)..T(wMax)
    Rectangle(normal, T, l, w, T(loc))
end

Rectangle(normal, lMin, lMax, wMin, wMax, loc) = Rectangle(; normal = normal, lMin = lMin, lMax = lMax, wMin = wMin, wMax = wMax, loc = loc)

function Rectangle(normal::Union{Val{:x}, Val{:y}, Val{:z}}, l::L, w::W, loc::C) where {L<:Real, W<:Real, C<:Real}
    T = float(promote_type(L,W,C))
    Rectangle(normal, T, T(l/2), T(w/2), T(loc))
end

@inline in(p::AbstractCoordinatePoint, r::RectangleX) = begin
    _in_y(p, r.l) && _in_z(p, r.w) && _isapprox_x(p, r.loc)
end

@inline in(p::AbstractCoordinatePoint, r::RectangleY) = begin
    _in_x(p, r.l) && _in_z(p, r.w) && _isapprox_y(p, r.loc)
end

@inline in(p::AbstractCoordinatePoint, r::RectangleZ) = begin
    _in_x(p, r.l) && _in_y(p, r.w) && _eq_z(p, r.loc)
end    

get_l_limits(r::Rectangle) = (_left_linear_interval(r.l), _right_linear_interval(r.l))
get_w_limits(r::Rectangle) = (_left_linear_interval(r.w), _right_linear_interval(r.w))

function sample(r::RectangleX{T}, Nsamps::NTuple{3,Int} = (2,2,2)) where {T}
    lMin::T, lMax::T = get_l_limits(r)
    wMin::T, wMax::T = get_w_limits(r)
    samples = [
        CartesianPoint{T}(r.loc,y,z)
        for y in (Nsamps[2] ≤ 1 ? lMin : range(lMin, lMax, length = Nsamps[2]))
        for z in (Nsamps[3] ≤ 1 ? wMin : range(wMin, wMax, length = Nsamps[3]))
    ]
end

function sample(r::RectangleY{T}, Nsamps::NTuple{3,Int} = (2,2,2)) where {T}
    lMin::T, lMax::T = get_l_limits(r)
    wMin::T, wMax::T = get_w_limits(r)
    samples = [
        CartesianPoint{T}(x,r.loc,z)
        for x in (Nsamps[1] ≤ 1 ? lMin : range(lMin, lMax, length = Nsamps[1]))
        for z in (Nsamps[3] ≤ 1 ? wMin : range(wMin, wMax, length = Nsamps[3]))
    ]
end

function sample(r::RectangleZ{T}, Nsamps::NTuple{3,Int} = (2,2,2)) where {T}
    lMin::T, lMax::T = get_l_limits(r)
    wMin::T, wMax::T = get_w_limits(r)
    samples = [
        CartesianPoint{T}(x,y,r.loc)
        for x in (Nsamps[1] ≤ 1 ? lMin : range(lMin, lMax, length = Nsamps[1]))
        for y in (Nsamps[2] ≤ 1 ? wMin : range(wMin, wMax, length = Nsamps[2]))
    ]
end

function get_vertices(r::RectangleX{T}) where {T}
    lMin::T, lMax::T = get_l_limits(r)
    wMin::T, wMax::T = get_w_limits(r)
    (CartesianPoint{T}(r.loc, lMin, wMin),
    CartesianPoint{T}(r.loc, lMax, wMin),
    CartesianPoint{T}(r.loc, lMin, wMax),
    CartesianPoint{T}(r.loc, lMax, wMax))
end

function get_vertices(r::RectangleY{T}) where {T}
    lMin::T, lMax::T = get_l_limits(r)
    wMin::T, wMax::T = get_w_limits(r)
    (CartesianPoint{T}(lMin, r.loc, wMin),
    CartesianPoint{T}(lMax, r.loc, wMin),
    CartesianPoint{T}(lMin, r.loc, wMax),
    CartesianPoint{T}(lMax, r.loc, wMax))
end

function get_vertices(r::RectangleZ{T}) where {T}
    lMin::T, lMax::T = get_l_limits(r)
    wMin::T, wMax::T = get_w_limits(r)
    (CartesianPoint{T}(lMin, wMin, r.loc),
    CartesianPoint{T}(lMax, wMin, r.loc),
    CartesianPoint{T}(lMin, wMax, r.loc),
    CartesianPoint{T}(lMax, wMax, r.loc))
end