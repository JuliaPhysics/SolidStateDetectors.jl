struct Rectangle{T,TN,TL,TW} <: AbstractSurfacePrimitive{T}
    normal::TN
    l::TL
    w::TW
    loc::T
    function Rectangle( ::Type{T},
                        normal::Union{Val{:x}, Val{:y}, Val{:z}},
                        l::Union{T, <:AbstractInterval{T}},
                        w::Union{T, <:AbstractInterval{T}},
                        loc::T) where {T}
        new{T,typeof(normal),typeof(l),typeof(w)}(normal,l,w,loc)
    end
end

#Constructors
Rectangle(b::Box{T}, normal::Val{:x}, x::Real) where {T} = Rectangle(T, normal, b.y, b.z, T(x))
Rectangle(b::Box{T}, normal::Val{:y}, y::Real) where {T} = Rectangle(T, normal, b.x, b.z, T(y))
Rectangle(b::Box{T}, normal::Val{:z}, z::Real) where {T} = Rectangle(T, normal, b.x, b.y, T(z))

function Rectangle(;normal = Val(:z), lMin = -1, lMax = 1, wMin = -1, wMax = 1, loc = 0)
    T = float(promote_type(typeof.((lMin, lMax, wMin, wMax, loc))...))
    l = lMax == -lMin ? T(lMax) : T(lMin)..T(lMax)
    w = wMax == -wMin ? T(wMax) : T(wMin)..T(wMax)
    Rectangle(T, normal, l, w, T(loc))
end

Rectangle(normal, lMin, lMax, wMin, wMax, loc) = Rectangle(; normal = normal, lMin = lMin, lMax = lMax, wMin = wMin, wMax = wMax, loc = loc)

function Rectangle(normal::Union{Val{:x}, Val{:y}, Val{:z}}, l::L, w::W, loc::C) where {L<:Real, W<:Real, C<:Real}
    T = float(promote_type(L,W,C))
    Rectangle(T, normal, T(l/2), T(w/2), T(loc))
end

@inline in(p::AbstractCoordinatePoint, r::Rectangle{<:Any, Val{:x}}) = begin
    _in_y(p, r.l) && _in_z(p, r.w) && _isapprox_x(p, r.loc)
end

@inline in(p::AbstractCoordinatePoint, r::Rectangle{<:Any, Val{:y}}) = begin
    _in_x(p, r.l) && _in_z(p, r.w) && _isapprox_y(p, r.loc)
end

@inline in(p::AbstractCoordinatePoint, r::Rectangle{<:Any, Val{:z}}) = begin
    _in_x(p, r.l) && _in_y(p, r.w) && _eq_z(p, r.loc)
end    

get_l_limits(r::Rectangle) = (_left_linear_interval(r.l), _right_linear_interval(r.l))
get_w_limits(r::Rectangle) = (_left_linear_interval(r.w), _right_linear_interval(r.w))

function sample(r::Rectangle{T, Val{:x}}, Nsamps::NTuple{3,Int} = (2,2,2)) where {T}
    lMin::T, lMax::T = get_l_limits(r)
    wMin::T, wMax::T = get_w_limits(r)
    samples = [
        CartesianPoint{T}(r.loc,y,z)
        for y in (Nsamps[2] ≤ 1 ? lMin : range(lMin, lMax, length = Nsamps[2]))
        for z in (Nsamps[3] ≤ 1 ? wMin : range(wMin, wMax, length = Nsamps[3]))
    ]
end

function sample(r::Rectangle{T, Val{:y}}, Nsamps::NTuple{3,Int} = (2,2,2)) where {T}
    lMin::T, lMax::T = get_l_limits(r)
    wMin::T, wMax::T = get_w_limits(r)
    samples = [
        CartesianPoint{T}(x,r.loc,z)
        for x in (Nsamps[1] ≤ 1 ? lMin : range(lMin, lMax, length = Nsamps[1]))
        for z in (Nsamps[3] ≤ 1 ? wMin : range(wMin, wMax, length = Nsamps[3]))
    ]
end

function sample(r::Rectangle{T, Val{:z}}, Nsamps::NTuple{3,Int} = (2,2,2)) where {T}
    lMin::T, lMax::T = get_l_limits(r)
    wMin::T, wMax::T = get_w_limits(r)
    samples = [
        CartesianPoint{T}(x,y,r.loc)
        for x in (Nsamps[1] ≤ 1 ? lMin : range(lMin, lMax, length = Nsamps[1]))
        for y in (Nsamps[2] ≤ 1 ? wMin : range(wMin, wMax, length = Nsamps[2]))
    ]
end


function get_r(r::Rectangle{T, Val{:x}}, φ::T) where {T}
    atol = geom_atol_zero(T)
    x::T = r.loc
    yMin::T, yMax::T = get_l_limits(r)
    sφ::T, cφ::T = sincos(φ)
    if x * cφ < 0 return () end
    tanφ::T = tan(φ)
    yMin - atol ≤ x * tanφ ≤ yMax + atol ? (x / cφ,) : ()
end

function get_r(r::Rectangle{T, Val{:y}}, φ::T) where {T}
    atol = geom_atol_zero(T)
    y::T = r.loc
    xMin::T, xMax::T = get_l_limits(r)
    sφ::T, cφ::T = sincos(φ)
    if y * sφ < 0 return () end
    cotφ::T = cot(φ)
    xMin - atol ≤ y * cotφ ≤ xMax + atol ? (y / sφ,) : ()
end

function sample(r::Union{Rectangle{T, Val{:x}}, Rectangle{T, Val{:y}}}, g::CylindricalTicksTuple{T}) where {T}
    samples = [
        CylindricalPoint{T}(R, φ, z)
        for φ in g.φ
        for R in get_r(r,φ)
        for z in _get_ticks(g.z, get_w_limits(r)...)
    ]
end

function get_r_ticks(r::Rectangle{T, Val{:z}}, g::CylindricalTicksTuple{T}, φ::T)::Vector{T} where {T}
    xMin::T, xMax::T = get_l_limits(r)
    yMin::T, yMax::T = get_w_limits(r)
    sφ::T, cφ::T = sincos(φ)
    tanφ::T = tan(φ)
    cotφ::T = inv(tanφ)
    R = sort!(filter!(R -> R > 0,
        [
            (xMin ≤ yMin * cotφ ≤ xMax ? yMin / sφ : T(NaN)),
            (xMin ≤ yMax * cotφ ≤ xMax ? yMax / sφ : T(NaN)),
            (yMin ≤ xMin * tanφ ≤ yMax ? xMin / cφ : T(NaN)),
            (yMin ≤ xMax * tanφ ≤ yMax ? xMax / cφ : T(NaN))
        ]
    ))
    ticks = if length(R) == 2
        _get_ticks(g.r, R...)
    elseif length(R) == 1
        _get_ticks(g.r, T(0), R[1])
    else
        []
    end
end

function sample(r::Rectangle{T, Val{:z}}, g::CylindricalTicksTuple{T}) where {T}
    samples = [
        CylindricalPoint{T}(R, φ, z)
        for z in _get_ticks(g.z, r.loc, r.loc)
        for φ in g.φ
        for R in get_r_ticks(r, g, φ)
    ]
end

function sample(r::Rectangle{T, Val{:x}}, g::CartesianTicksTuple{T}) where {T}
    samples = [
        CartesianPoint{T}(r.loc,y,z)
        for y in _get_ticks(g.y, get_l_limits(r)...)
        for z in _get_ticks(g.z, get_w_limits(r)...)
    ]
end

function sample(r::Rectangle{T, Val{:y}}, g::CartesianTicksTuple{T}) where {T}
    samples = [
        CartesianPoint{T}(x,r.loc,z)
        for x in _get_ticks(g.x, get_l_limits(r)...)
        for z in _get_ticks(g.z, get_w_limits(r)...)
    ]
end

function sample(r::Rectangle{T, Val{:z}}, g::CartesianTicksTuple{T}) where {T}
    samples = [
        CartesianPoint{T}(x,y,r.loc)
        for x in _get_ticks(g.x, get_l_limits(r)...)
        for y in _get_ticks(g.y, get_w_limits(r)...)
    ]
end

function get_vertices(r::Rectangle{T, Val{:x}}) where {T}
    lMin::T, lMax::T = get_l_limits(r)
    wMin::T, wMax::T = get_w_limits(r)
    (CartesianPoint{T}(r.loc, lMin, wMin),
    CartesianPoint{T}(r.loc, lMax, wMin),
    CartesianPoint{T}(r.loc, lMin, wMax),
    CartesianPoint{T}(r.loc, lMax, wMax))
end

function get_vertices(r::Rectangle{T, Val{:y}}) where {T}
    lMin::T, lMax::T = get_l_limits(r)
    wMin::T, wMax::T = get_w_limits(r)
    (CartesianPoint{T}(lMin, r.loc, wMin),
    CartesianPoint{T}(lMax, r.loc, wMin),
    CartesianPoint{T}(lMin, r.loc, wMax),
    CartesianPoint{T}(lMax, r.loc, wMax))
end

function get_vertices(r::Rectangle{T, Val{:z}}) where {T}
    lMin::T, lMax::T = get_l_limits(r)
    wMin::T, wMax::T = get_w_limits(r)
    (CartesianPoint{T}(lMin, wMin, r.loc),
    CartesianPoint{T}(lMax, wMin, r.loc),
    CartesianPoint{T}(lMin, wMax, r.loc),
    CartesianPoint{T}(lMax, wMax, r.loc))
end