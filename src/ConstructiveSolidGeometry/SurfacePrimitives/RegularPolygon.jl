struct RegularPolygon{N,T,TR} <: AbstractSurfacePrimitive{T}
    r::TR
    z::T
    function RegularPolygon( N::Integer,
                   ::Type{T},
                   r::Union{T, <:AbstractInterval{T}},
                   z::T) where {T}
        @assert typeof(N) <: Integer && N >= 3 "Failed constructing RegularPolygon with $(N) points in the polygon base (needs at least 3 points)."
        new{N,T,typeof(r)}(r, z)
    end
end

# Convenience functions
const RegularTriangle{T,TR} = RegularPolygon{3,T,TR}
const Square{T,TR}  = RegularPolygon{4,T,TR}
const RegularPentagon{T,TR} = RegularPolygon{5,T,TR}
const RegularHexagon{T,TR}  = RegularPolygon{6,T,TR}

RegularTriangle(args...) = RegularPolygon(3, args...)
Square(args...)  = RegularPolygon(4, args...)
RegularPentagon(args...) = RegularPolygon(5, args...)
RegularHexagon(args...)  = RegularPolygon(6, args...)

#Constructors
RegularPolygon(p::RegularPrism{N,T}; z = 0) where {N,T} = RegularPolygon(N, T, p.r, T(z))

function RegularPolygon(N::Integer; rInner = 0, rOuter = 1, z = 0)
    T = float(promote_type(typeof.((rInner, rOuter, z))...))
    r = rInner == 0 ? T(rOuter) : T(rInner)..T(rOuter)
    RegularPolygon(N, T, r, T(z))
end

RegularPolygon(N::Integer, rInner, rOuter, z) = RegularPolygon(N; rInner = rInner, rOuter = rOuter, z = z)

function RegularPolygon(N::Integer, rOuter::R, z::Z) where {R<:Real, Z<:Real}
    T = float(promote_type(R,Z))
    RegularPolygon( N, T, T(rOuter), T(z))
end


@inline in(p::CylindricalPoint, rp::RegularPolygon{N, T, <:Real}) where {N,T} = begin
    _eq_z(p, rp.z) && p.r * cos(T(π/N) - mod(p.φ, T(2π/N))) / cos(T(π/N)) <= rp.r
end

@inline in(p::CylindricalPoint, rp::RegularPolygon{N, T, <:AbstractInterval}) where {N,T} = begin
    _eq_z(p, rp.z) && p.r * cos(T(π/N) - mod(p.φ, T(2π/N))) / cos(T(π/N)) in rp.r
end

@inline in(p::CartesianPoint, rp::RegularPolygon) = in(CylindricalPoint(p), rp)

# Special case: CartesianPoint in RegularHexagon: use analytical formulas
@inline in(p::CartesianPoint, rp::RegularHexagon{T, <:Real}) where {T} =
    _eq_z(p, rp.z) && _in_x(p, rp.r * sqrt(T(3))/2) && _in_y(p, rp.r - abs(p.x)/sqrt(T(3)))

@inline in(p::CartesianPoint, rp::RegularHexagon{T, <:AbstractInterval{T}}) where {T} =
    _eq_z(p, rp.z) && abs(p.y) <= rp.r.right * sqrt(T(3))/2 &&
    (
        abs(p.y) >= rp.r.left * sqrt(T(3))/2 && _in_x(p, rp.r.right * T(0.5)) ||
        abs(p.x) in (rp.r.left - abs(p.y) /sqrt(T(3)))..(rp.r.right - abs(p.y)/sqrt(T(3)))
    )

get_r_limits(rp::RegularPolygon) = (_left_radial_interval(rp.r), _right_radial_interval(rp.r))

function sample(rp::RegularPolygon{N,T}, Nsamps::NTuple{3,Int} = (2,N+1,2))::Vector{CylindricalPoint{T}} where {N,T}
    rMin::T, rMax::T = get_r_limits(rp)
    samples = [
        CylindricalPoint{T}(r*cos(π/N)/cos(π/N - mod(φ, 2π/N)),φ,rp.z)
        for φ in (Nsamps[2] ≤ 1 ? 0 : range(0, 2π, length = Nsamps[2])) #ignore the N given in Nsamps and sample with N
        for r in (Nsamps[1] ≤ 1 ? rMax : range(rMin, rMax, length = Nsamps[1]))
        ]
end
